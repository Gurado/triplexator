
cmake_minimum_required (VERSION 2.6)
project (triplexator)

# The version number.
set (triplexator_VERSION_MAJOR 1)
set (triplexator_VERSION_MINOR 3)
set (triplexator_VERSION_PATCH 2)

# does this system provide the required dependencies?
include (${CMAKE_ROOT}/Modules/CheckFunctionExists.cmake)
find_package (Boost COMPONENTS iostreams)
find_package (OpenMP)

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_subdirectory (src) 

# ---------------------------------------------------------------------------
# build a CPack driven installer package
INCLUDE(InstallRequiredSystemLibraries)
SET(CPACK_CMAKE_GENERATOR "Unix Makefiles")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "triplexator: Finding nucleic acid triple helices")
SET(CPACK_PACKAGE_VENDOR "Fabian Buske")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENCE.txt")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
SET(CPACK_PACKAGE_VERSION_MAJOR "${triplexator_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${triplexator_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${triplexator_VERSION_PATCH}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "triplexator ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
IF(WIN32 AND NOT UNIX)
  # There is a bug in NSI that does not handle full unix paths properly. Make
  # sure there is at least one set of four (4) backlasshes.
#  SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\InstallIcon.bmp")
  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\triplexator.exe")
  SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} triplexator")
  SET(CPACK_NSIS_HELP_LINK "http:\\\\www.bioinformatics.org.au/triplexator")
  SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\www.bioinformatics.org.au/triplexator")
  SET(CPACK_NSIS_CONTACT "f.buske@uq.edu.au")
  SET(CPACK_NSIS_MODIFY_PATH ON)
ELSE(NOT WIN32 AND UNIX)
  SET(CPACK_STRIP_FILES "bin/triplexator")
  SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(WIN32 AND NOT UNIX)
SET(CPACK_PACKAGE_EXECUTABLES "triplexator" "triplexator")
INCLUDE(CPack)


