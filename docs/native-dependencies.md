# Native dependency policy

CPlantBox links SuiteSparse and SUNDIALS when `PIAFMUNCH=ON`.

`PIAFMUNCH` remains enabled by default to preserve the existing documented build flow:

```bash
cmake .
make install
```

## Providers

SuiteSparse and SUNDIALS providers are selected independently:

```bash
-DCPB_SUITESPARSE_PROVIDER=bundled|source|system
-DCPB_SUNDIALS_PROVIDER=bundled|source|system
```

Defaults remain compatible with the current source tree:

```bash
-DCPB_SUITESPARSE_PROVIDER=bundled
-DCPB_SUNDIALS_PROVIDER=bundled
```

Wheel builds override these defaults and use the `source` provider explicitly.

## Source provider

The source provider consumes a local prefix containing headers and static libraries built from pinned source inputs:

```bash
scripts/deps/source-native-deps.sh /tmp/cplantbox-native-deps
cmake . \
  -DCPB_SUITESPARSE_PROVIDER=source \
  -DCPB_SUNDIALS_PROVIDER=source \
  -DCPB_NATIVE_DEPS_PREFIX=/tmp/cplantbox-native-deps
make install
```

macOS uses the same source recipe, with deployment-target flags supplied by the wheel build before invoking the script.

Current pins:

| Dependency | Version | Source input |
| --- | --- | --- |
| SUNDIALS | 4.0.2 | Release tarball with SHA256 verification |
| SuiteSparse | 5.3.0 | Git tag `v5.3.0` verified at commit `e927f7a3fc82339755482e553df37d932ff30083` |

The source provider is the preferred path for wheel builds because it avoids relying on checked-in prebuilt archives or arbitrary system libraries.

## Bundled provider

The bundled provider uses the prebuilt archives and headers in the source tree:

```text
src/external/suitsparse/include/
src/external/suitsparse/lib/*.a
src/external/sundials/include/
src/external/sundials/lib/*.a
```

Those archives are kept for compatibility with the existing native build. They are Linux archives and are not used by the new macOS wheel path.

## System provider

The system provider uses headers and libraries discovered by CMake:

```bash
-DCPB_SUITESPARSE_PROVIDER=system
-DCPB_SUNDIALS_PROVIDER=system
```

Use `CMAKE_PREFIX_PATH`, compiler environment variables, or standard system install locations to make dependencies discoverable.

The system provider is explicit opt-in. Wheel builds should not accidentally depend on libraries from the CI image or a developer machine.

## Wheel policy

Wheel builds should have:

- explicit source-built native dependencies
- no runtime dependency on `libpython` / `Python.framework`
- no unexpected dynamic SuiteSparse/SUNDIALS dependencies
- no Homebrew, `/tmp`, user-home, or CI build-directory runtime paths

The wheel workflow enforces these with installed-wheel smoke tests and platform audit scripts.

## Publishing status

Wheel artifacts include the source-built SuiteSparse/SUNDIALS license and provenance files under `.dist-info/licenses/third_party/native-dependencies`. Publishing to TestPyPI/PyPI remains out of scope for this PR.
