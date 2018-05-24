# Link in dependencies necessary to generate plots
find_package(fplot)
find_package(ferror 1.3.0)
find_package(collections 1.1.0)

# Interpolation Example
add_executable(curvefit_interp_example curvefit_interp_example.f90)
target_link_libraries(curvefit_interp_example curvefit)

# LOWESS Example
add_executable(curvefit_lowess_example curvefit_lowess_example.f90)
target_link_libraries(curvefit_lowess_example curvefit fplot)

# Nonlinear Regression Example
add_executable(curvefit_nlreg_example curvefit_nlreg_example.f90)
target_link_libraries(curvefit_nlreg_example curvefit fplot)

# Calibration Example
add_executable(curvefit_cal_example curvefit_cal_example.f90)
target_link_libraries(curvefit_cal_example curvefit)

# Linear Interpolation Example
add_executable(curvefit_linear_interp_example curvefit_linear_interp_example.f90)
target_link_libraries(curvefit_linear_interp_example curvefit fplot)

# Polynomial Interpolation Example
add_executable(curvefit_poly_interp_example curvefit_poly_interp_example.f90)
target_link_libraries(curvefit_poly_interp_example curvefit fplot)

# Spline Interpolation Example
add_executable(curvefit_spline_interp_example curvefit_spline_interp_example.f90)
target_link_libraries(curvefit_spline_interp_example curvefit fplot)

# Moving Average Example
add_executable(curvefit_moving_avg_example curvefit_moving_avg_example.f90)
target_link_libraries(curvefit_moving_avg_example curvefit fplot)

# Linear Least-Squares Scalar Example
add_executable(curvefit_lls_scalar_example curvefit_lls_scalar_example.f90)
target_link_libraries(curvefit_lls_scalar_example curvefit fplot)

# ------------------------------------------------------------------------------
if (WIN32)
    # Ensure NONLIN, LINALG, FERROR, and BLAS libraries are copied to the runtime directory
    get_target_property(ferror_LibLocation ferror LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ferror_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )

    get_target_property(linalg_LibLocation linalg LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${linalg_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )

    get_target_property(blas_LibLocation blas LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${blas_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )

    get_target_property(nonlin_LibLocation nonlin LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${nonlin_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )

    # Copy shared libraries necessary to generate plots to the runtime directory
    get_target_property(ferror_LibLocation ferror LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ferror_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )

    get_target_property(fplot_LibLocation fplot LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${fplot_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )

    get_target_property(collections_LibLocation collections LOCATION)
    add_custom_command(TARGET curvefit_interp_example POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${collections_LibLocation} $<TARGET_FILE_DIR:curvefit_interp_example>
    )
endif()
