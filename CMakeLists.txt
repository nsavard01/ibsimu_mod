cmake_minimum_required(VERSION 3.16)

project(ibsimu
  VERSION 1.0.6
  DESCRIPTION "Ion Beam Simulation Package"
  LANGUAGES CXX
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# -----------------------------------------------------------------------
# Options -- mirror the old ./configure --with-*/--disable-* flags.
# Optional deps default ON ("check") like autotools did; if pkg-config
# can't find one, it's silently disabled below, same as before.
# MKL has no auto-detection (there's nothing to probe without MKLROOT
# being meaningful), so it stays opt-in, matching your remote build
# always passing --with-mkl explicitly.
# -----------------------------------------------------------------------

option(IBSIMU_WITH_GTK3     "Build the interactive Gtk3 plotter" ON)
option(IBSIMU_WITH_OPENGL   "Build the OpenGL/GtkGLExt 3D renderer (needs IBSIMU_WITH_GTK3)" ON)
option(IBSIMU_WITH_UMFPACK  "Build the UMFPACK direct solver" ON)
option(IBSIMU_WITH_CSG      "Build CSG solid support" ON)
option(IBSIMU_WITH_MKL      "Link against Intel MKL for BLAS" OFF)
option(IBSIMU_SIGSEGV_STACK "Print a stack trace on SIGSEGV" ON)

if(IBSIMU_WITH_OPENGL AND NOT IBSIMU_WITH_GTK3)
  message(WARNING "IBSIMU_WITH_OPENGL requires IBSIMU_WITH_GTK3 -- disabling IBSIMU_WITH_OPENGL")
  set(IBSIMU_WITH_OPENGL OFF)
endif()

# -----------------------------------------------------------------------
# Required dependencies
# -----------------------------------------------------------------------

find_package(PkgConfig REQUIRED)
find_package(ZLIB REQUIRED)          # replaces the autoconf-archive CHECK_ZLIB macro
find_package(Threads REQUIRED)       # replaces the autoconf-archive ACX_PTHREAD macro

pkg_check_modules(CAIRO      REQUIRED IMPORTED_TARGET cairo>=1.2.4)
pkg_check_modules(GSL        REQUIRED IMPORTED_TARGET gsl>=1.12)
pkg_check_modules(LIBPNG     REQUIRED IMPORTED_TARGET libpng)
pkg_check_modules(FONTCONFIG REQUIRED IMPORTED_TARGET fontconfig)
pkg_check_modules(FREETYPE   REQUIRED IMPORTED_TARGET freetype2)

# -----------------------------------------------------------------------
# Optional dependencies -- each disables its feature if not found,
# same "check" behavior configure.ac had by default.
# -----------------------------------------------------------------------

if(IBSIMU_WITH_GTK3)
  pkg_check_modules(GTK3 IMPORTED_TARGET gtk+-3.0)
  if(NOT GTK3_FOUND)
    message(STATUS "gtk+-3.0 not found -- disabling IBSIMU_WITH_GTK3 (and IBSIMU_WITH_OPENGL)")
    set(IBSIMU_WITH_GTK3 OFF)
    set(IBSIMU_WITH_OPENGL OFF)
  endif()
endif()

if(IBSIMU_WITH_OPENGL)
  pkg_check_modules(GTKGLEXT IMPORTED_TARGET gtkglext-3.0)
  if(NOT GTKGLEXT_FOUND)
    message(STATUS "gtkglext-3.0 not found -- disabling IBSIMU_WITH_OPENGL")
    set(IBSIMU_WITH_OPENGL OFF)
  endif()
endif()

if(IBSIMU_WITH_UMFPACK)
  find_package(UMFPACK)
  if(NOT UMFPACK_FOUND)
    message(STATUS "UMFPACK not found -- disabling IBSIMU_WITH_UMFPACK (set -DUMFPACK_ROOT=... if it's installed in a non-standard place)")
    set(IBSIMU_WITH_UMFPACK OFF)
  endif()
endif()

if(IBSIMU_WITH_CSG)
  pkg_check_modules(CSG IMPORTED_TARGET csg)
  if(NOT CSG_FOUND)
    message(STATUS "csg pkg-config module not found -- disabling IBSIMU_WITH_CSG")
    set(IBSIMU_WITH_CSG OFF)
  endif()
endif()

# -----------------------------------------------------------------------
# Intel MKL -- no pkg-config file, mirrors the hand-rolled --with-mkl
# block in configure.ac (same MKLROOT env var, same library set,
# threaded gnu_thread layer -- keep -fopenmp in CMAKE_CXX_FLAGS).
# -----------------------------------------------------------------------

if(IBSIMU_WITH_MKL)
  if(DEFINED ENV{MKLROOT})
    set(MKL_ROOT_DIR "$ENV{MKLROOT}")
  else()
    set(MKL_ROOT_DIR "/opt/intel/oneapi/mkl/latest")
  endif()

  if(NOT EXISTS "${MKL_ROOT_DIR}/include/mkl.h")
    message(FATAL_ERROR "IBSIMU_WITH_MKL is ON but mkl.h wasn't found under ${MKL_ROOT_DIR}. Set the MKLROOT environment variable or -DMKL_ROOT_DIR=...")
  endif()

  add_library(MKL::MKL INTERFACE IMPORTED)
  target_include_directories(MKL::MKL INTERFACE "${MKL_ROOT_DIR}/include")
  target_link_directories(MKL::MKL INTERFACE "${MKL_ROOT_DIR}/lib/intel64")
  target_link_options(MKL::MKL INTERFACE "-Wl,-rpath,${MKL_ROOT_DIR}/lib/intel64")
  target_link_libraries(MKL::MKL INTERFACE
    mkl_intel_lp64 mkl_gnu_thread mkl_core gomp dl ${CMAKE_THREAD_LIBS_INIT} m
  )
  message(STATUS "Intel MKL: ${MKL_ROOT_DIR}")
endif()

# -----------------------------------------------------------------------
# Feature checks -- replace the AC_CHECK_HEADERS/AC_CHECK_FUNCS/
# AC_CHECK_TYPES/AC_FUNC_STRERROR_R block in configure.ac. On any
# current Linux/glibc these all resolve true; kept as real checks
# rather than hardcoded so the build still behaves on other platforms.
# -----------------------------------------------------------------------

include(CheckIncludeFileCXX)
include(CheckSymbolExists)
include(CheckCXXSourceCompiles)

check_include_file_cxx(unistd.h HAVE_UNISTD_H)
check_symbol_exists(getcwd unistd.h HAVE_GETCWD)
check_include_file_cxx(sys/time.h HAVE_SYS_TIME_H)
check_symbol_exists(gettimeofday sys/time.h HAVE_GETTIMEOFDAY)
check_symbol_exists(clock_gettime time.h HAVE_CLOCK_GETTIME)

check_cxx_source_compiles("
  #include <signal.h>
  int main() { siginfo_t s; (void)s; return 0; }
" HAVE_SIGINFO_T)

check_cxx_source_compiles("
  #include <time.h>
  int main() { clockid_t c; (void)c; return 0; }
" HAVE_CLOCKID_T)

check_cxx_source_compiles("
  #include <time.h>
  int main() { struct timespec t; (void)t; return 0; }
" HAVE_STRUCT_TIMESPEC)

find_library(LIBRT_LIBRARY rt)
if(LIBRT_LIBRARY)
  set(HAVE_LIBRT 1)
endif()

check_cxx_source_compiles("
  #include <string.h>
  int main() { char buf[64]; char *p = strerror_r(0, buf, sizeof(buf)); (void)p; return 0; }
" STRERROR_R_CHAR_P)
if(STRERROR_R_CHAR_P)
  set(HAVE_STRERROR_R 1)
else()
  check_symbol_exists(strerror_r string.h HAVE_STRERROR_R)
endif()

# -----------------------------------------------------------------------
# Library
# -----------------------------------------------------------------------

add_subdirectory(src)

# -----------------------------------------------------------------------
# pkg-config file (kept for non-CMake consumers of the library)
# -----------------------------------------------------------------------

set(IBSIMU_LIB_NAME ibsimu)
set(IBSIMU_PC_INCLUDE_SUBDIR ibsimu)

set(IBSIMU_PC_REQUIRES_OPTIONAL "")
if(IBSIMU_WITH_GTK3)
  string(APPEND IBSIMU_PC_REQUIRES_OPTIONAL " gtk+-3.0")
endif()
if(IBSIMU_WITH_OPENGL)
  string(APPEND IBSIMU_PC_REQUIRES_OPTIONAL " gtkglext-3.0")
endif()
if(IBSIMU_WITH_CSG)
  string(APPEND IBSIMU_PC_REQUIRES_OPTIONAL " csg")
endif()

# NOTE: pkg_check_modules()'s *_LIBRARIES variables are bare library
# names (no -l, no CMake list flattening applied), and *_LDFLAGS /
# *_LDFLAGS_OTHER are themselves CMake ;-separated lists -- pasting
# either straight into a text file via configure_file() leaves literal
# semicolons in the output. Collect everything as a real CMake list
# first, then flatten it to a space-separated string right before
# writing the .pc file.
set(_ibsimu_pc_libs_list
  ${FREETYPE_LDFLAGS}
  ${FONTCONFIG_LDFLAGS}
  ${LIBPNG_LDFLAGS}
  -lz
  ${CMAKE_THREAD_LIBS_INIT}
)
set(_ibsimu_pc_cflags_list
  ${PTHREAD_CFLAGS}
)
if(IBSIMU_WITH_UMFPACK)
  list(APPEND _ibsimu_pc_libs_list -lumfpack)
  list(APPEND _ibsimu_pc_cflags_list "-I${UMFPACK_INCLUDE_DIRS}")
endif()

list(REMOVE_DUPLICATES _ibsimu_pc_libs_list)
string(REPLACE ";" " " IBSIMU_PC_LIBS_EXTRA "${_ibsimu_pc_libs_list}")
string(REPLACE ";" " " IBSIMU_PC_CFLAGS_EXTRA "${_ibsimu_pc_cflags_list}")

configure_file(
  "${CMAKE_SOURCE_DIR}/ibsimu.pc.in"
  "${CMAKE_BINARY_DIR}/ibsimu.pc"
  @ONLY
)
install(FILES "${CMAKE_BINARY_DIR}/ibsimu.pc"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

# -----------------------------------------------------------------------
# CMake package config, so downstream CMake projects can do
#   find_package(ibsimu REQUIRED)
#   target_link_libraries(your_target PRIVATE ibsimu::ibsimu)
# instead of hand-linking against the .pc file.
# -----------------------------------------------------------------------

install(EXPORT ibsimuTargets
  FILE ibsimuTargets.cmake
  NAMESPACE ibsimu::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ibsimu
)

configure_package_config_file(
  "${CMAKE_SOURCE_DIR}/cmake/ibsimuConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/ibsimuConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ibsimu
)
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/ibsimuConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)
install(FILES
  "${CMAKE_BINARY_DIR}/ibsimuConfig.cmake"
  "${CMAKE_BINARY_DIR}/ibsimuConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ibsimu
)

# -----------------------------------------------------------------------
# tests/ and doc/ CMakeLists.txt are a follow-up once this library
# build is verified against the autotools one on the remote box --
# see CMAKE_MIGRATION_PLAN.md.
# -----------------------------------------------------------------------

message(STATUS "")
message(STATUS "ibsimu ${PROJECT_VERSION} configuration")
message(STATUS "  Gtk3 plotter     : ${IBSIMU_WITH_GTK3}")
message(STATUS "  OpenGL renderer  : ${IBSIMU_WITH_OPENGL}")
message(STATUS "  UMFPACK solver   : ${IBSIMU_WITH_UMFPACK}")
message(STATUS "  CSG solid        : ${IBSIMU_WITH_CSG}")
message(STATUS "  Intel MKL        : ${IBSIMU_WITH_MKL}")
message(STATUS "  Install prefix   : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
