#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "PROJ::proj" for configuration "Release"
set_property(TARGET PROJ::proj APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(PROJ::proj PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libproj.25.9.2.0.dylib"
  IMPORTED_SONAME_RELEASE "@rpath/libproj.25.dylib"
  )

list(APPEND _cmake_import_check_targets PROJ::proj )
list(APPEND _cmake_import_check_files_for_PROJ::proj "${_IMPORT_PREFIX}/lib/libproj.25.9.2.0.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
