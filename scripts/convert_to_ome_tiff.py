#!/usr/bin/env python3
"""
convert_to_ome_tiff.py
~~~~~~~~~~~~~~~~~~~~~~
Convert microscopy images to Avivator-compatible pyramidal OME-TIFF
using QuPath CLI (Bio-Formats for reading, OME-TIFF for writing).

Supports multi-series formats (SCN, CZI, LIF, NDPI, VSI, ...).
Produces tiled pyramidal OME-TIFF files with offsets.json for
Viv / Avivator viewers.

Requirements
------------
- QuPath built with ``./gradlew installDist`` (or an installed QuPath on PATH)
- Python 3.9+
- tifffile  (``pip install tifffile``)

Usage
-----
    # Auto-discover and convert all series
    python convert_to_ome_tiff.py image.scn -o output/

    # Convert specific series
    python convert_to_ome_tiff.py image.scn -o output/ --series 0 2

    # Fluorescence only (skip brightfield/RGB series)
    python convert_to_ome_tiff.py image.scn -o output/ --fluorescence-only

    # Custom QuPath install path
    python convert_to_ome_tiff.py image.scn -o output/ \\
        --qupath /path/to/QuPath/bin/QuPath

    # QuPath from gradle build
    python convert_to_ome_tiff.py image.scn -o output/ \\
        --qupath ./build/install/QuPath/bin/QuPath \\
        --java-home /path/to/jdk25

Bugs fixed vs. original oncosuite script
-----------------------------------------
1. Uses valid compression (DEFAULT/LZW/ZLIB) — ZSTD is NOT supported by
   QuPath's Bio-Formats writer and silently fails with exit code 0.
2. Always checks output file existence after conversion (QuPath may return
   exit code 0 even when the conversion fails).
3. Uses gentler series filtering — the original excluded valid series via
   over-aggressive area/pixel-size heuristics.
4. Prints QuPath stdout/stderr when conversion fails for diagnosis.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from xml.etree import ElementTree as ET

try:
    import tifffile
except ImportError:
    sys.exit("tifffile is required: pip install tifffile")

# ---------------------------------------------------------------------------
# Valid QuPath OME-TIFF compression types (from OMEPyramidWriter.CompressionType)
# ---------------------------------------------------------------------------
VALID_COMPRESSIONS = ("UNCOMPRESSED", "DEFAULT", "JPEG", "J2K", "J2K_LOSSY", "LZW", "ZLIB")


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def _parse_ome_pixels(ome_xml: str | None):
    """Parse <Pixels> element from OME-XML, return dict or None."""
    if not ome_xml:
        return None
    try:
        root = ET.fromstring(ome_xml)
        for el in root.iter():
            if el.tag.endswith("Pixels"):
                return el
    except ET.ParseError:
        pass
    return None


def read_ome_metadata(path: str, downsample_factor: int = 1) -> dict:
    """Read key metadata from an OME-TIFF file."""
    meta: dict = {
        "width": 0, "height": 0, "n_channels": 1,
        "pixel_type": "", "pixel_size_um": None,
        "is_rgb": False, "area": 0,
    }
    try:
        with tifffile.TiffFile(path) as tif:
            page = tif.pages[0]
            w = int(page.imagewidth * downsample_factor)
            h = int(page.imagelength * downsample_factor)
            meta["width"] = w
            meta["height"] = h
            meta["area"] = w * h

            ome_xml = getattr(tif, "ome_metadata", None)
            pixels = _parse_ome_pixels(ome_xml)
            if pixels is not None:
                meta["n_channels"] = int(pixels.get("SizeC", "1"))
                ptype = (pixels.get("Type") or "").upper()
                meta["pixel_type"] = ptype
                meta["is_rgb"] = meta["n_channels"] == 3 and "UINT8" in ptype

                sx = pixels.get("PhysicalSizeX")
                sy = pixels.get("PhysicalSizeY")
                if sx and sy:
                    px = (float(sx) + float(sy)) / 2.0
                    meta["pixel_size_um"] = px / downsample_factor if downsample_factor > 1 else px
                elif sx:
                    meta["pixel_size_um"] = float(sx) / downsample_factor if downsample_factor > 1 else float(sx)
    except Exception as exc:
        print(f"  Warning: metadata read failed for {path}: {exc}")
    return meta


def get_mpp_from_ometiff(path: str) -> float | None:
    """Extract microns-per-pixel from OME-TIFF OME-XML."""
    try:
        with tifffile.TiffFile(path) as tif:
            pixels = _parse_ome_pixels(getattr(tif, "ome_metadata", None))
            if pixels is None:
                return None
            sx = pixels.get("PhysicalSizeX")
            sy = pixels.get("PhysicalSizeY")
            if sx and sy:
                return (float(sx) + float(sy)) / 2.0
            return float(sx or sy) if (sx or sy) else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Avivator compatibility check
# ---------------------------------------------------------------------------

def check_avivator_compat(path: str) -> tuple[bool, str]:
    """
    Validate OME-TIFF for Avivator/Viv compatibility.

    Rejects 12-bit, YCbCr, float/signed pixel data across ALL pyramid levels.
    """
    try:
        with tifffile.TiffFile(path) as t:
            for idx, page in enumerate(t.pages):
                bps = page.bitspersample
                if bps not in (8, 16, 32):
                    return False, f"page {idx}: BitsPerSample={bps}"

                pi = page.tags.get("PhotometricInterpretation")
                if pi and pi.value == 6:
                    return False, f"page {idx}: YCbCr unsupported"

                sf = page.tags.get("SampleFormat")
                if sf is not None:
                    vals = sf.value if hasattr(sf.value, "__iter__") and not isinstance(sf.value, str) else (sf.value,)
                    for v in vals:
                        if v not in (1, 4):
                            return False, f"page {idx}: SampleFormat={v}"
        return True, "ok"
    except Exception as exc:
        return False, str(exc)


# ---------------------------------------------------------------------------
# Offsets.json generation for Avivator
# ---------------------------------------------------------------------------

def generate_offsets_json(ome_tiff_path: str) -> str | None:
    """
    Generate offsets.json containing IFD byte offsets for Avivator.

    Avivator (Viv) can use these offsets to seek directly to any
    tile/level in the OME-TIFF without parsing the full file.

    The output is a list of lists: one inner list per resolution level,
    each containing the IFD offsets for the pages at that level.
    """
    try:
        with tifffile.TiffFile(ome_tiff_path) as tif:
            if not tif.series:
                return None

            series = tif.series[0]
            offsets_by_level: list[list[int]] = []
            for level in series.levels:
                level_offsets = []
                for page in level.pages:
                    if hasattr(page, "offset"):
                        level_offsets.append(page.offset)
                offsets_by_level.append(level_offsets)

        out_path = ome_tiff_path.replace(".ome.tif", ".offsets.json")
        with open(out_path, "w") as f:
            json.dump(offsets_by_level, f)
        return out_path
    except Exception as exc:
        print(f"  Warning: offsets generation failed: {exc}")
        return None


# ---------------------------------------------------------------------------
# Series discovery
# ---------------------------------------------------------------------------

def discover_series(
    qupath_cmd: list[str],
    input_path: str,
    work_dir: str,
    env: dict | None = None,
    max_series: int = 20,
) -> list[tuple[int, dict]]:
    """
    Discover valid series in a multi-series image by probing with
    ``convert-ome --downsample=100``.

    Returns list of (series_index, metadata_dict).
    """
    input_abs = os.path.abspath(input_path)
    series_info: list[tuple[int, dict]] = []

    for s in range(max_series):
        probe_path = os.path.join(work_dir, f".probe_s{s}.ome.tif")
        args = qupath_cmd + [
            "convert-ome",
            input_abs,
            probe_path,
            f"--series={s}",
            "--downsample=100",
            "-c", "LZW",
            "--overwrite",
        ]
        try:
            result = subprocess.run(
                args, capture_output=True, text=True, env=env, timeout=180,
            )
            if result.returncode != 0:
                break
            if not os.path.isfile(probe_path):
                stderr = (result.stderr or "")[:200]
                stdout = (result.stdout or "")[:200]
                print(f"  Series {s}: probe returned 0 but no file (stderr={stderr!r})")
                break
            meta = read_ome_metadata(probe_path, downsample_factor=100)
            series_info.append((s, meta))
        except subprocess.TimeoutExpired:
            print(f"  Timeout probing series {s}, stopping discovery")
            break
        except Exception as exc:
            print(f"  Error probing series {s}: {exc}")
            break
        finally:
            try:
                os.remove(probe_path)
            except OSError:
                pass

    return series_info


# ---------------------------------------------------------------------------
# Series filtering
# ---------------------------------------------------------------------------

def filter_series(
    series_info: list[tuple[int, dict]],
    fluorescence_only: bool = False,
    min_area_frac: float = 0.01,
) -> list[tuple[int, dict]]:
    """
    Exclude macro / label / thumbnail series that are very small relative
    to the largest series.  Optionally exclude RGB/brightfield.

    The 1 % default threshold is much gentler than the original script's
    5 % + pixel-size heuristic, which incorrectly excluded valid series
    from multi-series SCN files.
    """
    if not series_info:
        return []

    max_area = max(m["area"] for _, m in series_info)
    threshold = max_area * min_area_frac

    kept, excluded = [], []
    for idx, meta in series_info:
        if meta["area"] < threshold:
            excluded.append((idx, "too small (macro/label/thumbnail)"))
            continue
        if fluorescence_only and meta.get("is_rgb"):
            excluded.append((idx, "RGB/brightfield (fluorescence-only mode)"))
            continue
        kept.append((idx, meta))

    for idx, reason in excluded:
        print(f"  Excluding series {idx}: {reason}")

    return kept


# ---------------------------------------------------------------------------
# Single-series conversion
# ---------------------------------------------------------------------------

def convert_one_series(
    qupath_cmd: list[str],
    input_path: str,
    output_path: str,
    series_idx: int,
    compression: str = "DEFAULT",
    tile_size: int = 512,
    pyramid_scale: float = 4.0,
    env: dict | None = None,
    timeout: int | None = None,
) -> tuple[bool, str | None]:
    """
    Convert one series to pyramidal OME-TIFF via QuPath CLI.

    Returns ``(success, error_message | None)``.
    """
    input_abs = os.path.abspath(input_path)
    args = qupath_cmd + [
        "convert-ome",
        input_abs,
        output_path,
        f"--series={series_idx}",
        "--big-tiff=true",
        "-c", compression,
        "--tile-width", str(tile_size),
        "--tile-height", str(tile_size),
        "--overwrite",
    ]
    if pyramid_scale > 1:
        args += ["--pyramid-scale", str(pyramid_scale)]

    try:
        result = subprocess.run(
            args, capture_output=True, text=True, env=env, timeout=timeout,
        )

        # QuPath may return 0 even on failure (e.g. invalid compression).
        # Always check that the output file was actually created.
        if result.returncode != 0:
            snippet = (result.stderr or result.stdout or "")[:500]
            return False, f"exit code {result.returncode}: {snippet}"

        if not os.path.isfile(output_path):
            stdout = (result.stdout or "")[:500]
            stderr = (result.stderr or "")[:500]
            return False, (
                f"QuPath returned 0 but no output file was created.\n"
                f"  stdout: {stdout}\n  stderr: {stderr}"
            )

        return True, None
    except subprocess.TimeoutExpired:
        return False, "conversion timed out"
    except Exception as exc:
        return False, str(exc)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Convert microscopy images to Avivator-compatible "
            "pyramidal OME-TIFF using QuPath CLI."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s image.scn -o output/
  %(prog)s image.scn -o output/ --series 0 2
  %(prog)s image.scn -o output/ --fluorescence-only
  %(prog)s image.scn -o output/ --compression LZW --tile-size 256
        """,
    )
    parser.add_argument("input", help="Input image file")
    parser.add_argument(
        "-o", "--output-dir", default=".",
        help="Output directory (default: current dir)",
    )
    parser.add_argument(
        "--qupath", default="QuPath",
        help="QuPath executable (default: QuPath on PATH)",
    )
    parser.add_argument(
        "--java-home", default=None,
        help="JAVA_HOME for QuPath (e.g. path to JDK 25)",
    )
    parser.add_argument(
        "--series", type=int, nargs="*", default=None,
        help="Series indices to convert (default: auto-discover)",
    )
    parser.add_argument(
        "--fluorescence-only", action="store_true",
        help="Skip brightfield/RGB series",
    )
    parser.add_argument(
        "--compression", default="DEFAULT", choices=VALID_COMPRESSIONS,
        help="TIFF compression (default: DEFAULT — JPEG for RGB, LZW/ZLIB otherwise)",
    )
    parser.add_argument("--tile-size", type=int, default=512, help="Tile size (default: 512)")
    parser.add_argument(
        "--pyramid-scale", type=float, default=4.0,
        help="Pyramid downscale factor between levels (default: 4). "
             "Set to 0 to use the source image's native pyramid levels only.",
    )
    parser.add_argument("--max-series", type=int, default=20, help="Max series to probe (default: 20)")
    parser.add_argument("--parallel", type=int, default=4, help="Parallel conversions (default: 4)")
    parser.add_argument("--timeout", type=int, default=None, help="Timeout per series in seconds")
    parser.add_argument("--skip-offsets", action="store_true", help="Skip offsets.json generation")
    parser.add_argument(
        "--java-opts", default="-Xmx16G -Xms4G",
        help="JAVA_OPTS for QuPath (default: -Xmx16G -Xms4G)",
    )
    args = parser.parse_args()

    # --- validate input ---
    if not os.path.isfile(args.input):
        sys.exit(f"Error: file not found: {args.input}")

    file_size_mb = os.path.getsize(args.input) / (1024 * 1024)
    stem = Path(args.input).stem

    # --- resolve QuPath command ---
    qupath_cmd = args.qupath.strip().split()

    # --- build env for subprocess ---
    env = os.environ.copy()
    if args.java_home:
        env["JAVA_HOME"] = args.java_home
    existing = env.get("JAVA_OPTS", "")
    env["JAVA_OPTS"] = f"{args.java_opts} {existing}".strip()

    # --- create dirs ---
    os.makedirs(args.output_dir, exist_ok=True)
    work_dir = os.path.join(args.output_dir, ".convert_work")
    os.makedirs(work_dir, exist_ok=True)

    print(f"Input:  {args.input}  ({file_size_mb:.1f} MB)")
    print(f"Output: {args.output_dir}")
    print(f"QuPath: {' '.join(qupath_cmd)}")
    print(f"Compression: {args.compression}")

    # ------------------------------------------------------------------
    # Step 1 — discover or use explicit series
    # ------------------------------------------------------------------
    if args.series is not None:
        print(f"\nUsing specified series: {args.series}")
        series_to_convert: list[tuple[int, dict]] = [(s, {}) for s in args.series]
    else:
        print(f"\nDiscovering series (max {args.max_series})...")
        all_series = discover_series(
            qupath_cmd, args.input, work_dir, env=env, max_series=args.max_series,
        )
        if not all_series:
            sys.exit("Error: no readable series found in image.")

        print(f"  Found {len(all_series)} series:")
        for idx, meta in all_series:
            kind = "RGB/brightfield" if meta.get("is_rgb") else "fluorescence"
            px = meta.get("pixel_size_um")
            px_s = f", {px:.4f} um/px" if px else ""
            print(f"    Series {idx}: {meta['width']}x{meta['height']}, "
                  f"{meta['n_channels']}ch, {kind}{px_s}")

        series_to_convert = filter_series(
            all_series, fluorescence_only=args.fluorescence_only,
        )
        if not series_to_convert:
            sys.exit("Error: no series remaining after filtering.")

    # ------------------------------------------------------------------
    # Step 2 — convert each series to pyramidal OME-TIFF
    # ------------------------------------------------------------------
    n = len(series_to_convert)
    print(f"\nConverting {n} series (compression={args.compression}, tile={args.tile_size})...")

    results: list[dict] = []

    def _do_convert(series_idx: int) -> dict:
        out_name = f"{stem}.series{series_idx}.ome.tif"
        out_path = os.path.join(args.output_dir, out_name)

        t0 = time.time()
        ok, err = convert_one_series(
            qupath_cmd, args.input, out_path, series_idx,
            compression=args.compression, tile_size=args.tile_size,
            pyramid_scale=args.pyramid_scale,
            env=env, timeout=args.timeout,
        )
        elapsed = time.time() - t0

        if not ok:
            return {"series": series_idx, "ok": False, "error": err, "elapsed": elapsed}

        compat, reason = check_avivator_compat(out_path)
        if not compat:
            return {
                "series": series_idx, "ok": False,
                "error": f"Avivator-incompatible: {reason}",
                "path": out_path, "elapsed": elapsed,
            }

        offsets_path = None
        if not args.skip_offsets:
            offsets_path = generate_offsets_json(out_path)

        mpp = get_mpp_from_ometiff(out_path)
        size_mb = os.path.getsize(out_path) / (1024 * 1024)

        return {
            "series": series_idx, "ok": True, "path": out_path,
            "offsets": offsets_path, "mpp": mpp, "size_mb": size_mb,
            "elapsed": elapsed,
        }

    workers = min(args.parallel, n)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(_do_convert, idx): idx
            for idx, _ in series_to_convert
        }
        for fut in as_completed(futures):
            r = fut.result()
            results.append(r)
            s = r["series"]
            if r["ok"]:
                mpp_s = f", {r['mpp']:.4f} um/px" if r.get("mpp") else ""
                print(f"  Series {s}: OK  ({r['size_mb']:.1f} MB, {r['elapsed']:.1f}s{mpp_s})")
                print(f"    -> {r['path']}")
                if r.get("offsets"):
                    print(f"    -> {r['offsets']}")
            else:
                print(f"  Series {s}: FAILED  ({r['elapsed']:.1f}s)")
                print(f"    {r['error']}")

    # --- clean up work dir ---
    shutil.rmtree(work_dir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    ok_results = [r for r in results if r["ok"]]
    fail_results = [r for r in results if not r["ok"]]

    print(f"\n{'=' * 60}")
    print(f"Done: {len(ok_results)} succeeded, {len(fail_results)} failed")
    if ok_results:
        print("\nSuccessful:")
        for r in sorted(ok_results, key=lambda x: x["series"]):
            print(f"  Series {r['series']}: {r['path']}")
    if fail_results:
        print("\nFailed:")
        for r in sorted(fail_results, key=lambda x: x["series"]):
            print(f"  Series {r['series']}: {r['error']}")

    if not ok_results:
        sys.exit(1)


if __name__ == "__main__":
    main()
