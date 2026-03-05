# AGENTS.md

## Cursor Cloud specific instructions

**Product**: QuPath — open-source desktop application (JavaFX) for bioimage analysis and quantitative digital pathology. Single monolithic Gradle multi-module build; no backend services, databases, or Docker required.

### Build, Test, Run

Standard commands via Gradle wrapper (see also `build.gradle.kts` and CI at `.github/workflows/gradle.yml`):

- **Build + tests**: `./gradlew build`
- **Tests only**: `./gradlew test`
- **Run GUI**: `./gradlew run` (requires a display; `DISPLAY=:1` is set in the cloud VM with Xvfb)
- **Lint**: No separate linter; Javadoc warnings surface during `./gradlew build`

### Key gotchas

- **JDK 25 required**: The project uses Java toolchains (`libs.versions.toml` → `jdk = "25"`). The Gradle `foojay-resolver-convention` plugin auto-downloads JDK 25 (Temurin) on first build. The system JDK (21) is not used for compilation — Gradle handles toolchain resolution automatically.
- **First build is slow** (~4-5 min) due to downloading the Gradle distribution, JDK 25 toolchain, and all Maven dependencies. Subsequent builds use the Gradle daemon and caches.
- **JavaFX warnings**: `Gdk-WARNING: XSetErrorHandler() called with a GDK error trap pushed` is harmless under Xvfb.
- **Native access warnings**: JavaCPP and JavaFX emit `restricted method` warnings on JDK 25. These are cosmetic and do not affect functionality.
- **Script editor**: Accessible via `Automate > Script editor` in the GUI. Supports Groovy scripting for testing features interactively.
