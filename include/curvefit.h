// curvefit.h
#ifndef CURVEFIT_H_DEFINED
#define CURVEFIT_H_DEFINED

#include <stdbool.h>
#include "ferror.h"
#include "nonlin.h"

/** @brief An error flag denoting an improperly sized array. */
#define CF_ARRAY_SIZE_ERROR NL_ARRAY_SIZE_ERROR
/** @brief An error denoting that there is insufficient memory available. */
#define CF_OUT_OF_MEMORY_ERROR NL_OUT_OF_MEMORY_ERROR
/** @brief An error denoting that no data has been defined. */
#define CF_NO_DATA_DEFINED_ERROR 300
/** @brief An error flag denoting an invalid input. */
#define CF_INVALID_INPUT_ERROR NL_INVALID_INPUT_ERROR
/** @brief An error flag denoting a non-monotonic array was given when a
!! monotonic array was expected. */
#define CF_NONMONOTONIC_ARRAY_ERROR 301
/** @brief An error resulting from an invalid operation. */
#define CF_INVALID_OPERATION_ERROR NL_INVALID_OPERATION_ERROR
/** @brief An error resulting from a lack of convergence. */
#define CF_CONVERGENCE_ERROR NL_CONVERGENCE_ERROR
/** @brief An error indicating the user-requested tolerance is too small to be
!! practical for the problem at hand. */
#define CF_TOLERANCE_TOO_SMALL_ERROR NL_TOLERANCE_TOO_SMALL_ERROR

/** Indicates that the spline is quadratic over the interval under
!! consideration (beginning or ending interval).  This is equivalent to
!! allowing a "natural" boundary condition at either the initial or final
!! point. */
#define SPLINE_QUADRATIC_OVER_INTERVAL 1000
/** Indicates a known first derivative at either the beginning or ending 
!! point. */
#define SPLINE_KNOWN_FIRST_DERIVATIVE 1001
/** Indicates a known second derivative at either the beginning or ending 
!! point. */
#define SPLINE_KNOWN_SECOND_DERIVATIVE 1002
/** Indicates a continuous third derivative at either the beginning or ending 
!! point. */
#define SPLINE_CONTINUOUS_THIRD_DERIVATIVE 1003

/** @brief Describes a routine for finding the coefficients of a function
!! of one variable.
!!
!! @param[in] x The independent variable.
!! @param[in] n The number of coefficients in @p c.
!! @param[in] c An array of function coefficients.
!!
!! @result The value of the function at @p x.
 */
typedef double (*reg_fcn)(double x, int n, const double *c);


/** @brief A type encapsulating the Fortran linear_interp type. */
typedef struct {
    /** @brief A pointer to the Fortran linear_interp object. */
    void *ptr;
    /** @brief The size of the Fortran linear_interp object, in bytes. */
    int n;
} linear_interp;

/** @brief A type encapsulating the Fortran polynomial_interp type. */
typedef struct {
    /** @brief A pointer to the Fortran polynomial_interp object. */
    void *ptr;
    /** @brief The size of the Fortran polynomial_interp object, in bytes. */
    int n;
} polynomial_interp;

/** @brief A type encapsulating the Fortran spline_interp type. */
typedef struct {
    /** @brief A pointer to the Fortran spline_interp object. */
    void *ptr;
    /** @brief The size of the Fortran spline_interp object, in bytes. */
    int n;
} spline_interp;


#ifdef __cplusplus
extern "C" {
#endif

/** @brief Tests to see if an array is montonically increasing or decreasing.
!!
!! @param[in] n The number of elements in the array.
!! @param[in] x The array to test.
!!
!! @return Returns true if @p x is monotonic; else, false.
 */
bool is_monotonic(int n, const double *x);

/** @brief Initializes a new linear_inter object.
!!
!! @param[out] obj The linear_inter object to initialize.
!! @param[in] n The number of data points.
!! @param[in] x An N-element array containing the x-components of each data
!!  point.  This array must be monotonic (ascending or descending only).
!! @param[in] y An N-element array containing the y-components of each data
!!  point.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
!!      available.
!!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically 
!!      increasing or decreasing.
 */
void alloc_linear_interp(linear_interp *obj, int n, const double *x, 
                         const double *y, errorhandler *err);

/** @brief Frees resources held by a linear_inter object.
!!
!! @param[in,out] obj The linear_inter object.
 */
void free_linear_interp(linear_interp *obj);

/** @brief Performs a linear interpolation to determine the points @p y that
!! for the requested indendent variable values in @p x.
!!
!! @param[in] obj The linear_interp object.
!! @param[in] n The number of points to interpolate.
!! @param[in] x An N-element array containing the values of the independent
!!  variable at which to interpolate.
!! @param[out] y An N-element array where the interpolated values can be
!!  written.
 */
void linear_interpolate(const linear_interp *obj, int n, const double *x, 
                        double *y);

/** @brief Gets the number of points used by the interpolation object.
!!
!! @param[in] obj The linear_interp object.
!!
!! @return The number of points.
 */
int linear_interp_get_point_count(const linear_interp *obj);

/** @brief Gets a copy of the data points stored by the interpolation object.
!!
!! @param[in] obj The linear_interp object.
!! @param[in] n The size of the buffer arrays.
!! @param[out] x An N-element array where the x-coordinate data will be 
!!  written.
!! @param[out] y An N-element array where the y-coordinate data will be 
!!  written.
!!
!! @par Remarks
!! If @p n is different than the actual number of points that exist, the 
!! lesser of the two values will be utilized.  The interpolation object
!! can be queried to determine the quantity of stored points.
 */
void linear_interp_get_points(const linear_interp *obj, int n, double *x, 
                              double *y);

/** @brief Initializes a new polynomial_interp object.
!!
!! @param[out] obj The polynomial_interp object to initialize.
!! @param[in] n The number of data points.
!! @param[in] x An N-element array containing the x-components of each data
!!  point.  This array must be monotonic (ascending or descending only).
!! @param[in] y An N-element array containing the y-components of each data
!!  point.
!! @param[in] order The order of the interpolating polynomial.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
!!      available.
!!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically 
!!      increasing or decreasing.
!!  - CF_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
 */
void alloc_polynomial_interp(polynomial_interp *obj, int n, const double *x, 
                             const double *y, int order, errorhandler *err);

/** @brief Frees resources held by a polynomial_interp object.
!!
!! @param[in,out] obj The polynomial_interp object.
 */
void free_polynomial_interp(polynomial_interp *obj);

/** @brief Performs a polynomial interpolation to determine the points @p y 
!! that for the requested indendent variable values in @p x.
!!
!! @param[in] obj The polynomial_interp object.
!! @param[in] n The number of points to interpolate.
!! @param[in] x An N-element array containing the values of the independent
!!  variable at which to interpolate.
!! @param[out] y An N-element array where the interpolated values can be
!!  written.
 */
void polynomial_interpolate(const polynomial_interp *obj, int n, 
                            const double *x, double *y);

/** @brief Gets the number of points used by the interpolation object.
!!
!! @param[in] obj The polynomial_interp object.
!!
!! @return The number of points.
 */
int polynomial_interp_get_point_count(const polynomial_interp *obj);

/** @brief Gets a copy of the data points stored by the interpolation object.
!!
!! @param[in] obj The polynomial_interp object.
!! @param[in] n The size of the buffer arrays.
!! @param[out] x An N-element array where the x-coordinate data will be 
!!  written.
!! @param[out] y An N-element array where the y-coordinate data will be 
!!  written.
!!
!! @par Remarks
!! If @p n is different than the actual number of points that exist, the 
!! lesser of the two values will be utilized.  The interpolation object
!! can be queried to determine the quantity of stored points.
 */
void polynomial_interp_get_points(const polynomial_interp *obj, int n, 
                                  double *x, double *y);

/** @brief Initializes a new spline_interp object.
!!
!! @param[out] obj The spline_interp object to initialize.
!! @param[in] n The number of data points.
!! @param[in] x An N-element array containing the x-components of each data
!!  point.  This array must be monotonic (ascending or descending only).
!! @param[in] y An N-element array containing the y-components of each data
!!  point.
!! @param[in] ibcbeg An input that defines the nature of the 
!!  boundary condition at the beginning of the spline.  If an invalid
!!  parameter is used, the code defaults to SPLINE_QUADRATIC_OVER_INTERVAL.
!!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
!!      initial interval.  No value is required for @p ybcbeg.
!!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
!!      initial point is provided in @p ybcbeg.
!!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
!!      initial point is provided in @p ybcbeg.
!!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
!!      at x(2).  No value is required for @p ybcbeg.
!! @param[in] ybcbeg If needed, the value of the initial point boundary
!!  condition.  If not needed, this parameter is ignored.
!! @param[in] ibcend An input that defines the nature of the 
!!  boundary condition at the end of the spline.  If an invalid
!!  parameter is used, the code defaults to SPLINE_QUADRATIC_OVER_INTERVAL.
!!  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
!!      final interval.  No value is required for @p ybcend.
!!  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
!!      initial point is provided in @p ybcend.
!!  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
!!      initial point is provided in @p ybcend.
!!  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
!!      at x(n-1).  No value is required for @p ybcend.
!! @param[in] ybcend If needed, the value of the final point boundary
!!  condition.  If not needed, this parameter is ignored.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
!!      available.
!!  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically 
!!      increasing or decreasing.
 */
void alloc_spline_interp(spline_interp *obj, int n, const double *x, 
                         const double *y, int ibcbeg, double ybcbeg, int ibcend,
                         double ybcend, errorhandler *err);

/** @brief Frees resources held by a spline_interp object.
!!
!! @param[in,out] obj The spline_interp object.
 */
void free_spline_interp(spline_interp *obj);

/** @brief Performs a spline interpolation to determine the points @p y 
!! that for the requested indendent variable values in @p x.
!!
!! @param[in] obj The spline_interp object.
!! @param[in] n The number of points to interpolate.
!! @param[in] x An N-element array containing the values of the independent
!!  variable at which to interpolate.
!! @param[out] y An N-element array where the interpolated values can be
!!  written.
 */
void spline_interpolate(const spline_interp *obj, int n, const double *x, 
                        double *y);

/** @brief Gets the number of points used by the interpolation object.
!!
!! @param[in] obj The spline_interp object.
!!
!! @return The number of points.
 */
int spline_interp_get_point_count(const spline_interp *obj);

/** @brief Gets a copy of the data points stored by the interpolation object.
!!
!! @param[in] obj The spline_interp object.
!! @param[in] n The size of the buffer arrays.
!! @param[out] x An N-element array where the x-coordinate data will be 
!!  written.
!! @param[out] y An N-element array where the y-coordinate data will be 
!!  written.
!!
!! @par Remarks
!! If @p n is different than the actual number of points that exist, the 
!! lesser of the two values will be utilized.  The interpolation object
!! can be queried to determine the quantity of stored points.
 */
void spline_interp_get_points(const spline_interp *obj, int n, double *x, 
                              double *y);

/** @brief Computes the interpolated first derivative.
!!
!! @param[in] obj The spline_interp object.
!! @param[in] n The number of points to interpolate.
!! @param[in] x An N-element array containing the values of the independent
!!  variable at which to interpolate.
!! @param[out] y An N-element array where the interpolated values can be
!!  written.
 */
void spline_interp_diff1(const spline_interp *obj, int n, double *x, double *y);

/** @brief Computes the interpolated second derivative.
!!
!! @param[in] obj The spline_interp object.
!! @param[in] n The number of points to interpolate.
!! @param[in] x An N-element array containing the values of the independent
!!  variable at which to interpolate.
!! @param[out] y An N-element array where the interpolated values can be
!!  written.
 */
void spline_interp_diff2(const spline_interp *obj, int n, double *x, double *y);

#ifdef __cplusplus
}
#endif
#endif // CURVEFIT_H_DEFINED
