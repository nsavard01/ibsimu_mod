#[=======================================================================[.rst:
FindUMFPACK
-----------

UMFPACK (part of SuiteSparse) generally ships without a pkg-config file,
so it can't be found with pkg_check_modules() the way the rest of
ibsimu's dependencies are. This module does a plain header/library search
instead.

Result variables
^^^^^^^^^^^^^^^^^

``UMFPACK_FOUND``
  true if umfpack.h and libumfpack were both found
``UMFPACK_INCLUDE_DIRS``
``UMFPACK_LIBRARIES``

Imported target
^^^^^^^^^^^^^^^^

``UMFPACK::UMFPACK``

Hints
^^^^^

Set ``UMFPACK_ROOT`` (cmake -DUMFPACK_ROOT=... or environment variable)
if UMFPACK is installed somewhere find_path/find_library won't look by
default, e.g. a SuiteSparse build under /opt or /usr/local.
#]=======================================================================]

find_path(UMFPACK_INCLUDE_DIR
  NAMES umfpack.h
  HINTS ${UMFPACK_ROOT} ENV UMFPACK_ROOT
  PATH_SUFFIXES include include/suitesparse suitesparse
)

find_library(UMFPACK_LIBRARY
  NAMES umfpack
  HINTS ${UMFPACK_ROOT} ENV UMFPACK_ROOT
  PATH_SUFFIXES lib lib64
)

# UMFPACK also needs the rest of SuiteSparse's support libraries at link
# time on most distros. These are optional finds -- if a distro's
# libumfpack.so already carries these as transitive deps (common when
# installed from a package manager), the extra target_link_libraries
# entries are harmless no-ops.
find_library(AMD_LIBRARY NAMES amd)
find_library(SUITESPARSECONFIG_LIBRARY NAMES suitesparseconfig)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK
  REQUIRED_VARS UMFPACK_LIBRARY UMFPACK_INCLUDE_DIR
)

if(UMFPACK_FOUND AND NOT TARGET UMFPACK::UMFPACK)
  add_library(UMFPACK::UMFPACK UNKNOWN IMPORTED)
  set_target_properties(UMFPACK::UMFPACK PROPERTIES
    IMPORTED_LOCATION "${UMFPACK_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIR}"
  )
  foreach(_extra_lib AMD_LIBRARY SUITESPARSECONFIG_LIBRARY)
    if(${_extra_lib})
      set_property(TARGET UMFPACK::UMFPACK APPEND PROPERTY
        INTERFACE_LINK_LIBRARIES "${${_extra_lib}}")
    endif()
  endforeach()
endif()

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY AMD_LIBRARY SUITESPARSECONFIG_LIBRARY)
set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})
set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY})
