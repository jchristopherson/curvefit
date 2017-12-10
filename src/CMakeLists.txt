# Locate the dependencies
add_subdirectory(external/nonlin)

# Define the CURVEFIT source files
set(curvefit_sources curvefit_c_binding.f90 curvefit_calibration.f90
    curvefit_core.f90 curvefit_interp.f90 curvefit_regression.f90
    curvefit_statistics.f90)

# Build the library
add_library(curvefit ${curvefit_sources})
target_link_libraries(curvefit nonlin)
