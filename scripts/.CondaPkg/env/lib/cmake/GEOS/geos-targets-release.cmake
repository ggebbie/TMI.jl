#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GEOS::geos" for configuration "Release"
set_property(TARGET GEOS::geos APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GEOS::geos PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgeos.3.12.2.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libgeos.3.12.2.dylib"
  )

list(APPEND _cmake_import_check_targets GEOS::geos )
list(APPEND _cmake_import_check_files_for_GEOS::geos "${_IMPORT_PREFIX}/lib/libgeos.3.12.2.dylib" )

# Import target "GEOS::geos_c" for configuration "Release"
set_property(TARGET GEOS::geos_c APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GEOS::geos_c PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "GEOS::geos"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libgeos_c.1.18.2.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libgeos_c.1.dylib"
  )

list(APPEND _cmake_import_check_targets GEOS::geos_c )
list(APPEND _cmake_import_check_files_for_GEOS::geos_c "${_IMPORT_PREFIX}/lib/libgeos_c.1.18.2.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
