# Define the CURVEFIT source files
set(curvefit_sources curvefit_c_binding.f90 curvefit_calibration.f90
    curvefit_core.f90 curvefit_interp.f90 curvefit_regression.f90
    curvefit_statistics.f90)

# Build the library
add_library(curvefit ${curvefit_sources})
target_link_libraries(curvefit
  ${NONLIN_LIBRARIES}
  ${LINALG_LIBRARIES}
  ${FERROR_LIBRARIES}
)

# ------------------------------------------------------------------------------
# INSTALLATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Define target information
set_property(TARGET curvefit PROPERTY VERSION ${CURVEFIT_VERSION})
set_property(TARGET curvefit PROPERTY SOVERSION ${CURVEFIT_MAJOR_VERSION})
set_property(TARGET curvefit PROPERTY INTERFACE_curvefit_MAJOR_VERSION ${CURVEFIT_MAJOR_VERSION})
set_property(TARGET curvefit APPEND PROPERTY COMPATIBLE_INTERFACE_STRING curvefit_MAJOR_VERSION)

# Locate the "include" directories
set(curvefit_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

# Define the installation instructions
install(TARGETS curvefit EXPORT curvefitTargets
    RUNTIME DESTINATION curvefit/bin
    LIBRARY DESTINATION curvefit/lib
    ARCHIVE DESTINATION curvefit/lib
    CONFIGURATIONS Release
    INCLUDES DESTINATION curvefit/include
)
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}/curvefit
    COMPONENT Devel
)

# Include the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION curvefit/doc)
install(FILES ${PROJECT_SOURCE_DIR}/doc/refman.pdf DESTINATION curvefit/doc)

# Define the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/curvefitConfigVersion.cmake"
    VERSION ${CURVEFIT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT curvefitTargets
    FILE "${CMAKE_BINARY_DIR}/curvefitTargets.cmake"
)

# Define the configuration file
configure_file(
    "${PROJECT_SOURCE_DIR}/curvefitConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/curvefitConfig.cmake"
    COPYONLY
)

set(ConfigPackageLocation lib/cmake/curvefit)
install(
    EXPORT curvefitTargets
    FILE curvefitTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)
install(
    FILES
        "${CMAKE_BINARY_DIR}/curvefitConfig.cmake"
        "${CMAKE_BINARY_DIR}/curvefitConfigVersion.cmake"
    DESTINATION
        ${ConfigPackageLocation}
    COMPONENT
        Devel
)
