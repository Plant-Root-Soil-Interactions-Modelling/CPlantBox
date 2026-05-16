# Native dependency policy

CPlantBox links SuiteSparse and SUNDIALS only when the optional PiafMunch code is enabled.

```cmake
PIAFMUNCH=ON
```

`PIAFMUNCH` is currently enabled by default to preserve the documented native install workflow.

## Providers

The native dependency provider is selected independently for SuiteSparse and SUNDIALS:

```bash
-DCPB_SUITESPARSE_PROVIDER=bundled|system
-DCPB_SUNDIALS_PROVIDER=bundled|system
```

The default is:

```bash
-DCPB_SUITESPARSE_PROVIDER=bundled
-DCPB_SUNDIALS_PROVIDER=bundled
```

This keeps the existing `cmake . && make install` workflow working without requiring users to install SuiteSparse or SUNDIALS manually.

Here, "native default" means the behavior a developer or user gets when they run the documented CMake build without provider flags:

```bash
cmake .
make install
```

The long-term preferred default is source-built dependencies: build pinned SuiteSparse/SUNDIALS source inputs, link against those outputs, and keep `system` as an explicit override for developers and distribution packagers. The current bundled/prebuilt path is transitional compatibility for the existing Linux build.

## Bundled provider

The bundled provider uses the prebuilt archives and headers in the source tree:

```text
src/external/suitsparse/include/
src/external/suitsparse/lib/*.a
src/external/sundials/include/
src/external/sundials/lib/*.a
```

The currently bundled archives are Linux ELF archives. They are suitable for the current Linux native baseline CI, but they are not a general solution for macOS, Windows, or all wheel platforms.

Detected bundled versions:

| Dependency | Version |
| --- | --- |
| SUNDIALS | 4.0.2 |
| SuiteSparse | 5.3.0 |
| KLU | 1.3.9 |
| AMD | 2.4.6 |
| BTF | 1.2.6 |
| COLAMD | 2.9.6 |

## System provider

The system provider uses headers and libraries discovered by CMake. It verifies the required headers and libraries for the PiafMunch/SUNDIALS/KLU integration and fails during configure if any are missing:

```bash
-DCPB_SUITESPARSE_PROVIDER=system
-DCPB_SUNDIALS_PROVIDER=system
```

Use `CMAKE_PREFIX_PATH`, compiler environment variables, or standard system install locations to make the dependencies discoverable.

The system provider is opt-in. Release wheel builds must not accidentally depend on arbitrary libraries from a developer machine or CI image.

The system provider is intended for SUNDIALS/SuiteSparse versions compatible with the current PiafMunch integration. Arbitrary newer system packages may not be API-compatible.

## Wheel policy

Platform wheels should be built with explicit and audited native dependency choices.

For now:

- native/developer default: bundled provider
- wheel jobs: explicit provider selection per platform
- no accidental system native dependencies in release wheels
- no runtime dependency on `libpython`
- no unresolved native symbols in installed shared libraries

Before release-wheel publishing, wheel builds should use pinned source-built native dependencies rather than the checked-in prebuilt archives. That should be done explicitly in the wheel/dependency CI PRs, not implicitly in local builds.

## Intended provider lifecycle

The intended end state is:

```text
source = default for native builds and release wheels
system = explicit opt-in override for developers and distro packagers
prebuilt/bundled archives = removed
```

The current CMake provider branches should be treated as follows:

- bundled/prebuilt branch: compatibility path only, used to preserve the existing Linux `cmake . && make install` behavior until source builds are ready
- system branch: keep as an explicit opt-in escape hatch
- future source branch: preferred path and eventual default

It is safe to delete the bundled/prebuilt branch from `src/CMakeLists.txt` and remove the checked-in archives only after all of these are true:

1. SuiteSparse and SUNDIALS are built from pinned source inputs in CI.
2. The source provider works for the documented native build with `PIAFMUNCH=ON`.
3. Wheel-oriented builds use the source provider and no longer reference `src/external/*/lib/*.a`.
4. Linux native CI and wheel smoke tests pass with the source provider.
5. Native dependency audits show no missing libraries, no unresolved symbols, and no unexpected `libpython` dependency.
6. Third-party license/provenance metadata for bundled source-built artifacts is complete.
7. Documentation has been updated so users know how the source provider is obtained and built.

When those conditions are met, the cleanup should remove:

```text
src/external/suitsparse/lib/*.a
src/external/sundials/lib/*.a
```

and the corresponding CMake bundled/prebuilt lookup code. Do not remove the system-provider code unless the project explicitly decides to stop supporting external dependency installs.

## Link-time checks

Linux shared-library builds use a no-undefined linker check so missing native symbols fail during build instead of during `import plantbox`.

The Ubuntu baseline CI also audits the installed extension and `libCPlantBox` with:

```bash
ldd -r
readelf -d
```

The audit fails on:

- `not found`
- unresolved symbols
- `NEEDED` entries for `libpythonX.Y`
- dynamic `NEEDED` entries for bundled SuiteSparse/SUNDIALS libraries

## License and provenance status

SUNDIALS license metadata is present in the source tree:

```text
src/external/sundials/include/sundials/LICENSE
src/external/sundials/include/sundials/NOTICE
```

The vendored SuiteSparse subset currently does not include complete license/provenance files, although its headers reference upstream license files.

Do not publish release wheels that bundle SuiteSparse/SUNDIALS artifacts until third-party license and provenance metadata is complete.
