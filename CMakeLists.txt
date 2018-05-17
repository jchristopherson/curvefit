# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.0)
project(curvefit C CXX Fortran)

# Define version information
set(CURVEFIT_MAJOR_VERSION 1)
set(CURVEFIT_MINOR_VERSION 2)
set(CURVEFIT_PATCH_VERSION 0)
set(CURVEFIT_VERSION ${CURVEFIT_MAJOR_VERSION}.${CURVEFIT_MINOR_VERSION}.${CURVEFIT_PATCH_VERSION})

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

# By default, shared library
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Locate Dependencies
find_package(ferror 1.3.0)
find_package(linalg 1.5.0)
find_package(nonlin 1.3.0)

# Define dependency directories
set(NONLIN_DIRECTORY src/external/nonlin)
set(LINALG_DIRECTORY ${NONLIN_DIRECTORY}/src/external/linalg)
set(FERROR_DIRECTORY ${LINALG_DIRECTORY}/src/external/ferror)

# If FERROR is not installed on the system, build the default implementation
if (${ferror_FOUND})
    message(STATUS "An acceptable version of FERROR (v" ${ferror_VERSION} ") was found, and will be utilized.")
    set(FERROR_LIBRARIES ferror)
    get_target_property(ferror_LibLocation ferror LOCATION)
endif()

# If LINALG is not installed on the system, build the default implementation
if (${linalg_FOUND})
    message(STATUS "An acceptable version of LINALG (v" ${linalg_VERSION} ") was found, and will be utilized.")
    set(LINALG_LIBRARIES linalg)
    get_target_property(linalg_LibLocation linalg LOCATION)
endif()

# If NONLIN is not installed on the system, build the default implementation
if (${nonlin_FOUND})
    message(STATUS "An acceptable version of NONLIN (v" ${nonlin_VERSION} ") was found, and will be utilized.")
    set(NONLIN_LIBRARIES nonlin)
    get_target_property(nonlin_LibLocation nonlin LOCATION)
else()
    message(STATUS "NONLIN not found.  The default implementation will be used.")
    add_subdirectory(${NONLIN_DIRECTORY})
    set(nonlin_INCLUDE_DIRS ${NONLIN_DIRECTORY}/include)
    set(NONLIN_LIBRARIES nonlin)
    set(nonlin_LibLocation ${nonlin_BINARY_DIR})
    include_directories(${LINALG_DIRECTORY}/include)
    include_directories(${FERROR_DIRECTORY}/include)
endif()

# Include the dependency module files
include_directories(${nonlin_INCLUDE_DIRS} ${linalg_INCLUDE_DIRS} ${ferror_INCLUDE_DIRS})

# Export all symbols on Windows when building shared libraries
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

# Locate the module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/include)

# Define output directories, if undefined
if (NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  message(STATUS "CURVEFIT library output directories undefined.  Using default directories.")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    # Debug Build
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
    endif ()
  elseif (CMAKE_BUILD_TYPE MATCHES Release)
    # Release Build
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Release)
    endif ()
  else ()
    # Default Condition
    if (BUILD_SHARED_LIBS)
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
    else()
      set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/Debug)
    endif ()
  endif ()
endif ()

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set(CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -Wall -Wno-c-binding-type")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Locate the source directory
add_subdirectory(src)

# ------------------------------------------------------------------------------
# EXAMPLES
# ------------------------------------------------------------------------------
option(BUILD_CURVEFIT_EXAMPLES "Build CURVEFIT examples?" OFF)
if (BUILD_CURVEFIT_EXAMPLES)
  # Inform the user we're building the examples
  message(STATUS "Building CURVEFIT examples.")

  # Build the examples
  add_subdirectory(examples)
endif()

# ------------------------------------------------------------------------------
# TESTS
# ------------------------------------------------------------------------------
option(BUILD_CURVEFIT_TESTS "Build CURVEFIT tests?" OFF)
if (BUILD_CURVEFIT_TESTS)
    # Inform the user we're building the tests
    message(STATUS "Building CURVEFIT tests.")

    # Build the tests
    add_subdirectory(tests)
endif()
