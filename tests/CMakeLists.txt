# Define the source files
set(curvefit_test_sources
    curvefit_test.f90
    curvefit_test_interp.f90
    curvefit_test_statistics.f90
    ${PROJECT_SOURCE_DIR}/src/external/nonlin/src/external/linalg/tests/test_core.f90
    curvefit_test_regression.f90
    curvefit_test_calibration.f90
)

# Build the library
add_executable(curvefit_test ${curvefit_test_sources})
target_link_libraries(curvefit_test curvefit ${LINALG_LIBRARIES})

# Run the tests
add_custom_command(
    OUTPUT curvefit_tests
    DEPENDS ${curvefit_test_sources}
    COMMAND curvefit_test
)
add_custom_target(run_curvefit_tests ALL DEPENDS curvefit_tests)

# Ensure NONLIN, LINALG, FERROR, and BLAS libraries are copied to the runtime directory
get_target_property(ferror_LibLocation ferror LOCATION)
add_custom_command(TARGET curvefit_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${ferror_LibLocation} $<TARGET_FILE_DIR:curvefit_test>
)

get_target_property(linalg_LibLocation linalg LOCATION)
add_custom_command(TARGET curvefit_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${linalg_LibLocation} $<TARGET_FILE_DIR:curvefit_test>
)

get_target_property(blas_LibLocation blas LOCATION)
add_custom_command(TARGET curvefit_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${blas_LibLocation} $<TARGET_FILE_DIR:curvefit_test>
)

get_target_property(nonlin_LibLocation nonlin LOCATION)
add_custom_command(TARGET curvefit_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${nonlin_LibLocation} $<TARGET_FILE_DIR:curvefit_test>
)
