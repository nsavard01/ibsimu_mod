# src/CMakeLists.txt -- the ibsimu library target.
#
# Source discovery: core sources are picked up automatically with
# file(GLOB ... CONFIGURE_DEPENDS) so adding/removing/renaming a .cpp in
# src/ is picked up on the next build with no CMakeLists.txt edit --
# CONFIGURE_DEPENDS makes CMake re-glob and reconfigure automatically
# whenever the directory listing changes (CMake >= 3.12, all modern
# generators). The four optional feature groups (GTK3/OpenGL/UMFPACK/CSG)
# stay as explicit lists below since they're small, rarely change, and
# whether they compile at all depends on a build option -- mixing that
# into a glob would be more fragile, not less.

file(GLOB IBSIMU_ALL_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB IBSIMU_ALL_HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

set(IBSIMU_GTK3_SOURCES
  gtkfielddiagdialog.cpp
  gtkfielddiagexportdialog.cpp
  gtkfielddiagwindow.cpp
  gtkframewindow.cpp
  gtkgeom3dwindow.cpp
  gtkgeomwindow.cpp
  gtkhardcopy.cpp
  gtkparticlediagdialog.cpp
  gtkparticlediagexportdialog.cpp
  gtkparticlediagwindow.cpp
  gtkplotter.cpp
  gtkpreferences.cpp
  gtkwindow.cpp
  icons.cpp
  softwarerenderer.cpp
)
set(IBSIMU_OPENGL_SOURCES glrenderer.cpp)
set(IBSIMU_UMFPACK_SOURCES epot_umfpacksolver.cpp)
set(IBSIMU_CSG_SOURCES csgobject_solid.cpp)

set(IBSIMU_OPTIONAL_SOURCES
  ${IBSIMU_GTK3_SOURCES}
  ${IBSIMU_OPENGL_SOURCES}
  ${IBSIMU_UMFPACK_SOURCES}
  ${IBSIMU_CSG_SOURCES}
)

set(IBSIMU_CORE_SOURCES ${IBSIMU_ALL_SOURCES})
foreach(_opt_src ${IBSIMU_OPTIONAL_SOURCES})
  list(REMOVE_ITEM IBSIMU_CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${_opt_src}")
endforeach()

# Public headers: everything except icon_*.hpp, which are embedded XPM
# icon data for the Gtk3 UI, not part of the library's API.
set(IBSIMU_PUBLIC_HEADERS ${IBSIMU_ALL_HEADERS})
list(FILTER IBSIMU_PUBLIC_HEADERS EXCLUDE REGEX "/icon_[^/]*\\.hpp$")

# ---------------------------------------------------------------------
# Generated config.h (replaces autoheader's src/config.h)
# ---------------------------------------------------------------------

set(GTK3 ${IBSIMU_WITH_GTK3})
set(OPENGL ${IBSIMU_WITH_OPENGL})
set(USE_MKL ${IBSIMU_WITH_MKL})
set(CSG ${IBSIMU_WITH_CSG})
set(SIGSEGV_STACK ${IBSIMU_SIGSEGV_STACK})

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  @ONLY
)

# ---------------------------------------------------------------------
# Library target
# ---------------------------------------------------------------------

add_library(ibsimu ${IBSIMU_CORE_SOURCES})
add_library(ibsimu::ibsimu ALIAS ibsimu)

target_compile_features(ibsimu PUBLIC cxx_std_14)

set_target_properties(ibsimu PROPERTIES
  VERSION   ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  OUTPUT_NAME ibsimu
)

# config.h lives in the build dir; sources #include "config.h" (quoted),
# so it just needs to be on the include path -- no source tree pollution.
target_include_directories(ibsimu
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/ibsimu>
)

target_compile_options(ibsimu PRIVATE -Wall)

target_link_libraries(ibsimu
  PUBLIC
    ZLIB::ZLIB
    Threads::Threads
    PkgConfig::CAIRO
    PkgConfig::GSL
    PkgConfig::LIBPNG
    PkgConfig::FONTCONFIG
    PkgConfig::FREETYPE
)

if(IBSIMU_WITH_GTK3)
  target_sources(ibsimu PRIVATE ${IBSIMU_GTK3_SOURCES})
  target_link_libraries(ibsimu PUBLIC PkgConfig::GTK3)
endif()

if(IBSIMU_WITH_OPENGL)
  target_sources(ibsimu PRIVATE ${IBSIMU_OPENGL_SOURCES})
  target_link_libraries(ibsimu PUBLIC PkgConfig::GTKGLEXT)
endif()

if(IBSIMU_WITH_UMFPACK)
  target_sources(ibsimu PRIVATE ${IBSIMU_UMFPACK_SOURCES})
  target_link_libraries(ibsimu PUBLIC UMFPACK::UMFPACK)
endif()

if(IBSIMU_WITH_CSG)
  target_sources(ibsimu PRIVATE ${IBSIMU_CSG_SOURCES})
  target_link_libraries(ibsimu PUBLIC PkgConfig::CSG)
endif()

if(IBSIMU_WITH_MKL)
  target_link_libraries(ibsimu PUBLIC MKL::MKL)
endif()

if(HAVE_LIBRT)
  target_link_libraries(ibsimu PUBLIC rt)
endif()

# ---------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------

install(TARGETS ibsimu
  EXPORT ibsimuTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES ${IBSIMU_PUBLIC_HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ibsimu
)

# config.h is generated (see configure_file() above), not part of the
# glob -- install it alongside the public headers same as autotools did
# (it was listed in src/Makefile.am's include_HEADERS). error.hpp and
# others #include "config.h" directly, so consumers building against the
# installed headers need this file present, not just available at
# library-build time.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/config.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ibsimu
)
