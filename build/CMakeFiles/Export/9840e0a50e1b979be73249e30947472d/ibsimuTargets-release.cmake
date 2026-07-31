#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ibsimu::ibsimu" for configuration "Release"
set_property(TARGET ibsimu::ibsimu APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(ibsimu::ibsimu PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libibsimu.a"
  )

list(APPEND _cmake_import_check_targets ibsimu::ibsimu )
list(APPEND _cmake_import_check_files_for_ibsimu::ibsimu "${_IMPORT_PREFIX}/lib64/libibsimu.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
