# ibsimu: autotools → CMake migration plan

## What the current build actually does

The repo builds with GNU Autotools: `configure.ac` at the root, plus four
`Makefile.am` files (`/`, `src/`, `tests/`, `doc/`). `autoreconf -fi` has to
regenerate `configure`, `aclocal.m4`, libtool scaffolding (`compile`,
`test-driver`, `ltmain.sh`) and `src/config.h.in` from scratch every time,
because none of that generated output is committed — only the `.am`/`.ac`
sources are.

Three of the checks in `configure.ac` (`CHECK_ZLIB`, `CHECK_UMFPACK`,
`ACX_PTHREAD`) aren't defined anywhere in this repo. There's no `m4/`
directory. They only work because `aclocal` pulls them in from the
system-wide `autoconf-archive` package at build time — so the build is
silently dependent on whatever autoconf-archive version happens to be
installed on the machine you run `autoreconf` on. That's most of the
"complex and weird" feeling: the real macro definitions aren't visible in
the repo at all.

Dependencies checked in `configure.ac`, required unless noted:

| Dependency | How it's found today | Notes |
|---|---|---|
| pkg-config | `PKG_PROG_PKG_CONFIG` | hard requirement |
| zlib | `CHECK_ZLIB` (autoconf-archive macro) | required |
| libpng | manual `pkg-config --exists` | required |
| fontconfig | manual `pkg-config --exists` | required |
| freetype2 | manual `pkg-config --exists` | required |
| cairo >= 1.2.4 | `PKG_CHECK_MODULES` | required |
| GSL >= 1.12 | manual `pkg-config --exists` | required |
| pthreads | `ACX_PTHREAD` (autoconf-archive macro) | required |
| gtk+-3.0 | `PKG_CHECK_MODULES`, default "check" | optional (`GTK3`) |
| gtkglext-3.0 | `PKG_CHECK_MODULES`, default "check" | optional (`OPENGL`) |
| UMFPACK | `CHECK_UMFPACK` (autoconf-archive macro) | optional |
| csg | `PKG_CHECK_MODULES`, default "check" | optional |
| Intel MKL | hand-rolled `--with-mkl` block, no pkg-config, uses `$MKLROOT` | optional, on by default in your remote build |

There's also a hand-written `--disable-sigsegv_stack` flag (default on),
and a batch of feature checks (`siginfo_t`, `clockid_t`, `struct timespec`,
`getcwd`, `gettimeofday`, `clock_gettime`, `librt`, `strerror_r`) that feed
into a generated `src/config.h`. `src/Makefile.am` conditionally adds
source files for GTK3 / OpenGL / UMFPACK / CSG using Automake's
`if GTK3 ... endif` blocks tied to `AM_CONDITIONAL`.

## Why CMake is simpler here

- No `autoreconf`/`aclocal`/`libtool` regeneration step — `cmake -S . -B build` reads the `CMakeLists.txt` tree directly, nothing to bootstrap.
- `find_package(PkgConfig)` + `pkg_check_modules(... IMPORTED_TARGET ...)` replaces the manual `pkg-config --exists` + `LIBS=... CPPFLAGS=...` string concatenation for cairo, GSL, libpng, fontconfig, freetype2, gtk+-3.0, gtkglext-3.0, csg — each becomes one line plus `target_link_libraries(... PkgConfig::CAIRO)`.
- `find_package(ZLIB REQUIRED)` and `find_package(Threads REQUIRED)` are built into CMake itself, so `CHECK_ZLIB` and `ACX_PTHREAD` (the two macros silently borrowed from autoconf-archive) disappear entirely — no dependency on what's in `/usr/share/aclocal` on the build machine.
- Optional features become `option(IBSIMU_WITH_GTK3 "..." ON)` style flags instead of `--with-x`/`--without-x` configure flags — same idea, more standard, and visible in `cmake -LH`.
- Conditional sources are a plain `if(IBSIMU_WITH_GTK3) target_sources(...) endif()` instead of Automake's separate `if GTK3 ... endif` block syntax.
- `src/config.h` generation is a `configure_file(config.h.in config.h)` call, same concept as autoheader but no separate `autoheader` invocation.
- Out-of-source builds, multiple build directories (Debug/Release/Sanitizer), and IDE integration (CLion, VS Code CMake Tools) work for free.
- `make check` becomes `ctest`, with the same per-test pass/fail reporting but native parallelism (`ctest -j`).

## Proposed structure

```
CMakeLists.txt          top-level: project(), options, find_package/pkg_check_modules,
                         add_subdirectory(src [tests] [doc]), install + package export
src/CMakeLists.txt      library target, conditional sources, config.h generation
tests/CMakeLists.txt    one add_executable()+add_test() per test program
doc/CMakeLists.txt      optional Doxygen target via find_package(Doxygen)
cmake/                  small helper modules, only if needed (see UMFPACK below)
```

Top-level would expose options mirroring the current `--with-*` flags:
`IBSIMU_WITH_GTK3`, `IBSIMU_WITH_OPENGL`, `IBSIMU_WITH_UMFPACK`,
`IBSIMU_WITH_CSG`, `IBSIMU_WITH_MKL`, `IBSIMU_SIGSEGV_STACK` — each
default `ON` except MKL, matching current "check by default" behavior,
except MKL which you always pass explicitly.

Target remote build command becomes:

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=/opt/ibsimu/1.0.6 \
      -DIBSIMU_WITH_MKL=ON \
      -DCMAKE_CXX_FLAGS="-O3 -fopenmp"
cmake --build build -j
sudo cmake --install build
```

That's the direct replacement for `autoreconf -fi && ./configure ... && make -j && sudo make install`.

## Points that need a decision or a closer look before/while implementing

- **UMFPACK has no `.pc` file on most systems.** `CHECK_UMFPACK` today probably falls back to header/library search. CMake needs either a small `cmake/FindUMFPACK.cmake` (a `find_path`/`find_library` pair, maybe 15 lines) or `-DUMFPACK_ROOT=...` — this is the one spot that isn't a one-line `pkg_check_modules` swap.
- **MKL linking.** The current block hand-links `-lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp` against `$MKLROOT`. Recent oneAPI MKL ships its own `MKLConfig.cmake` under `$MKLROOT/lib/cmake/mkl`, so `find_package(MKL CONFIG)` is an option — more robust (matches threading layer automatically) but a bigger change in behavior. Plan is to start by porting the manual link line as-is, and only move to `find_package(MKL)` as a later cleanup if you want it.
- **gtkglext-3.0** is long unmaintained. Worth confirming it's still actually available/used on your remote box, or dropping `OPENGL` support in this fork rather than porting it.
- **Library versioning.** Autotools uses libtool version-info `0:1:0`. CMake's `VERSION`/`SOVERSION` properties aren't a direct translation of libtool's scheme. Since this is your personal fork (not meant to coexist with upstream `ibsimu` installs), simplest is to drop the `-1.0.6dev` suffix baked into the library/pkgconfig filenames and just version normally (`libibsimu.so.1.0.6`), unless something else on the remote machine links against the exact current `.so` name.
- **`src/config.h`** is currently listed directly in `Makefile.am`'s `SOURCES`, even though it's a generated file — that's an autotools quirk that should not carry over; under CMake it's purely a build-directory artifact, never committed.

## Rollout steps

1. Add `CMakeLists.txt` + `src/CMakeLists.txt` only. Build on the remote box, diff the resulting shared library's exported symbols against the autotools build, run your existing test binaries against it manually.
2. Add `tests/CMakeLists.txt`, wire into `ctest`, confirm the same tests pass as `make check` currently reports.
3. Add `doc/CMakeLists.txt` for the Doxygen target (optional, low risk).
4. Add `install()` rules + `install(EXPORT ...)`/package config so downstream CMake projects can `find_package(ibsimu)` instead of hand-linking against the `.pc` file. Keep generating a `.pc` file too, for any non-CMake consumers.
5. Once the CMake build is verified working end-to-end on the remote machine, delete the autotools files: `configure.ac`, all four `Makefile.am`, `compile`, `test-driver`, `ibsimu-1.0.6dev.pc.in`, and update `.gitignore` for `build/`.
6. Update `README`/`INSTALL` with the new command sequence.

Let me know if you want me to go ahead and write the actual `CMakeLists.txt` files next — I'd start with step 1 (library only) so you can verify it builds on the remote box before touching tests/docs or deleting anything.
