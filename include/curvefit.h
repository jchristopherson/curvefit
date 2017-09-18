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
 * monotonic array was expected. */
#define CF_NONMONOTONIC_ARRAY_ERROR 301
/** @brief An error resulting from an invalid operation. */
#define CF_INVALID_OPERATION_ERROR NL_INVALID_OPERATION_ERROR
/** @brief An error resulting from a lack of convergence. */
#define CF_CONVERGENCE_ERROR NL_CONVERGENCE_ERROR
/** @brief An error indicating the user-requested tolerance is too small to be
 * practical for the problem at hand. */
#define CF_TOLERANCE_TOO_SMALL_ERROR NL_TOLERANCE_TOO_SMALL_ERROR

/** Indicates that the spline is quadratic over the interval under
 * consideration (beginning or ending interval).  This is equivalent to
 * allowing a "natural" boundary condition at either the initial or final
 * point. */
#define SPLINE_QUADRATIC_OVER_INTERVAL 1000
/** Indicates a known first derivative at either the beginning or ending 
 * point. */
#define SPLINE_KNOWN_FIRST_DERIVATIVE 1001
/** Indicates a known second derivative at either the beginning or ending 
 * point. */
#define SPLINE_KNOWN_SECOND_DERIVATIVE 1002
/** Indicates a continuous third derivative at either the beginning or ending 
 * point. */
#define SPLINE_CONTINUOUS_THIRD_DERIVATIVE 1003

/** @brief Describes a routine for finding the coefficients of a function
 * of one variable.
 *
 * @param x The independent variable.
 * @param n The number of coefficients in @p c.
 * @param c An array of function coefficients.
 *
 * @result The value of the function at @p x.
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

/** @brief A type encapsulating the Fortran lowess_smoothing type. */
typedef struct {
    /** @brief A pointer to the Fortran lowess_smoothing object. */
    void *ptr;
    /** @brief The size of the Fortran lowess_smoothing object, in bytes. */
    int n;
} lowess_smoothing;

/** @brief A type encapsulating the Fortran nonlinear_regression type. */
typedef struct {
    /** @brief A pointer to the Fortran nonlinear_regression object. */
    void *ptr;
    /** @brief The size of the Fortran nonlinear_regression object, in bytes. */
    int n;
} nonlinear_regression;


#ifdef __cplusplus
extern "C" {
#endif

/** @brief Tests to see if an array is montonically increasing or decreasing.
 *
 * @param n The number of elements in the array.
 * @param x The array to test.
 *
 * @return Returns true if @p x is monotonic; else, false.
 */
bool is_monotonic(int n, const double *x);

/** @brief Initializes a new linear_inter object.
 *
 * @param obj The linear_inter object to initialize.
 * @param n The number of data points.
 * @param x An N-element array containing the x-components of each data
 *  point.  This array must be monotonic (ascending or descending only).
 * @param y An N-element array containing the y-components of each data
 *  point.
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 *  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically 
 *      increasing or decreasing.
 */
void alloc_linear_interp(linear_interp *obj, int n, const double *x, 
                         const double *y, errorhandler *err);

/** @brief Frees resources held by a linear_inter object.
 *
 * @param obj The linear_inter object.
 */
void free_linear_interp(linear_interp *obj);

/** @brief Performs a linear interpolation to determine the points @p y that
 * for the requested indendent variable values in @p x.
 *
 * @param obj The linear_interp object.
 * @param n The number of points to interpolate.
 * @param x An N-element array containing the values of the independent
 *  variable at which to interpolate.
 * @param y An N-element array where the interpolated values can be
 *  written.
 */
void linear_interpolate(const linear_interp *obj, int n, const double *x, 
                        double *y);

/** @brief Gets the number of points used by the interpolation object.
 *
 * @param obj The linear_interp object.
 *
 * @return The number of points.
 */
int linear_interp_get_point_count(const linear_interp *obj);

/** @brief Gets a copy of the data points stored by the interpolation object.
 *
 * @param obj The linear_interp object.
 * @param n The size of the buffer arrays.
 * @param x An N-element array where the x-coordinate data will be 
 *  written.
 * @param y An N-element array where the y-coordinate data will be 
 *  written.
 *
 * @par Remarks
 * If @p n is different than the actual number of points that exist, the 
 * lesser of the two values will be utilized.  The interpolation object
 * can be queried to determine the quantity of stored points.
 */
void linear_interp_get_points(const linear_interp *obj, int n, double *x, 
                              double *y);

/** @brief Initializes a new polynomial_interp object.
 *
 * @param obj The polynomial_interp object to initialize.
 * @param n The number of data points.
 * @param x An N-element array containing the x-components of each data
 *  point.  This array must be monotonic (ascending or descending only).
 * @param y An N-element array containing the y-components of each data
 *  point.
 * @param order The order of the interpolating polynomial.
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 *  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically 
 *      increasing or decreasing.
 *  - CF_INVALID_INPUT_ERROR: Occurs if @p order is less than 1.
 */
void alloc_polynomial_interp(polynomial_interp *obj, int n, const double *x, 
                             const double *y, int order, errorhandler *err);

/** @brief Frees resources held by a polynomial_interp object.
 *
 * @param obj The polynomial_interp object.
 */
void free_polynomial_interp(polynomial_interp *obj);

/** @brief Performs a polynomial interpolation to determine the points @p y 
 * that for the requested indendent variable values in @p x.
 *
 * @param obj The polynomial_interp object.
 * @param n The number of points to interpolate.
 * @param x An N-element array containing the values of the independent
 *  variable at which to interpolate.
 * @param y An N-element array where the interpolated values can be
 *  written.
 */
void polynomial_interpolate(const polynomial_interp *obj, int n, 
                            const double *x, double *y);

/** @brief Gets the number of points used by the interpolation object.
 *
 * @param obj The polynomial_interp object.
 *
 * @return The number of points.
 */
int polynomial_interp_get_point_count(const polynomial_interp *obj);

/** @brief Gets a copy of the data points stored by the interpolation object.
 *
 * @param obj The polynomial_interp object.
 * @param n The size of the buffer arrays.
 * @param x An N-element array where the x-coordinate data will be 
 *  written.
 * @param y An N-element array where the y-coordinate data will be 
 *  written.
 *
 * @par Remarks
 * If @p n is different than the actual number of points that exist, the 
 * lesser of the two values will be utilized.  The interpolation object
 * can be queried to determine the quantity of stored points.
 */
void polynomial_interp_get_points(const polynomial_interp *obj, int n, 
                                  double *x, double *y);

/** @brief Initializes a new spline_interp object.
 *
 * @param obj The spline_interp object to initialize.
 * @param n The number of data points.
 * @param x An N-element array containing the x-components of each data
 *  point.  This array must be monotonic (ascending or descending only).
 * @param y An N-element array containing the y-components of each data
 *  point.
 * @param ibcbeg An input that defines the nature of the 
 *  boundary condition at the beginning of the spline.  If an invalid
 *  parameter is used, the code defaults to SPLINE_QUADRATIC_OVER_INTERVAL.
 *  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
 *      initial interval.  No value is required for @p ybcbeg.
 *  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
 *      initial point is provided in @p ybcbeg.
 *  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
 *      initial point is provided in @p ybcbeg.
 *  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
 *      at x(2).  No value is required for @p ybcbeg.
 * @param ybcbeg If needed, the value of the initial point boundary
 *  condition.  If not needed, this parameter is ignored.
 * @param ibcend An input that defines the nature of the 
 *  boundary condition at the end of the spline.  If an invalid
 *  parameter is used, the code defaults to SPLINE_QUADRATIC_OVER_INTERVAL.
 *  - SPLINE_QUADRATIC_OVER_INTERVAL: The spline is quadratic over its
 *      final interval.  No value is required for @p ybcend.
 *  - SPLINE_KNOWN_FIRST_DERIVATIVE: The spline's first derivative at its
 *      initial point is provided in @p ybcend.
 *  - SPLINE_KNOWN_SECOND_DERIVATIVE: The spline's second derivative at its
 *      initial point is provided in @p ybcend.
 *  - SPLINE_CONTINUOUS_THIRD_DERIVATIVE: The third derivative is continuous
 *      at x(n-1).  No value is required for @p ybcend.
 * @param ybcend If needed, the value of the final point boundary
 *  condition.  If not needed, this parameter is ignored.
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 *  - CF_NONMONOTONIC_ARRAY_ERROR: Occurs if @p x is not monotonically 
 *      increasing or decreasing.
 */
void alloc_spline_interp(spline_interp *obj, int n, const double *x, 
                         const double *y, int ibcbeg, double ybcbeg, int ibcend,
                         double ybcend, errorhandler *err);

/** @brief Frees resources held by a spline_interp object.
 *
 * @param obj The spline_interp object.
 */
void free_spline_interp(spline_interp *obj);

/** @brief Performs a spline interpolation to determine the points @p y 
 * that for the requested indendent variable values in @p x.
 *
 * @param obj The spline_interp object.
 * @param n The number of points to interpolate.
 * @param x An N-element array containing the values of the independent
 *  variable at which to interpolate.
 * @param y An N-element array where the interpolated values can be
 *  written.
 */
void spline_interpolate(const spline_interp *obj, int n, const double *x, 
                        double *y);

/** @brief Gets the number of points used by the interpolation object.
 *
 * @param obj The spline_interp object.
 *
 * @return The number of points.
 */
int spline_interp_get_point_count(const spline_interp *obj);

/** @brief Gets a copy of the data points stored by the interpolation object.
 *
 * @param obj The spline_interp object.
 * @param n The size of the buffer arrays.
 * @param x An N-element array where the x-coordinate data will be 
 *  written.
 * @param y An N-element array where the y-coordinate data will be 
 *  written.
 *
 * @par Remarks
 * If @p n is different than the actual number of points that exist, the 
 * lesser of the two values will be utilized.  The interpolation object
 * can be queried to determine the quantity of stored points.
 */
void spline_interp_get_points(const spline_interp *obj, int n, double *x, 
                              double *y);

/** @brief Computes the interpolated first derivative.
 *
 * @param obj The spline_interp object.
 * @param n The number of points to interpolate.
 * @param x An N-element array containing the values of the independent
 *  variable at which to interpolate.
 * @param y An N-element array where the interpolated values can be
 *  written.
 */
void spline_interp_diff1(const spline_interp *obj, int n, double *x, double *y);

/** @brief Computes the interpolated second derivative.
 *
 * @param obj The spline_interp object.
 * @param n The number of points to interpolate.
 * @param x An N-element array containing the values of the independent
 *  variable at which to interpolate.
 * @param y An N-element array where the interpolated values can be
 *  written.
 */
void spline_interp_diff2(const spline_interp *obj, int n, double *x, double *y);

/** @brief Computes the mean of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @return The mean of @p x.
 */
double mean(int n, const double *x);

/** @brief Computes the median of a data set.
 *
 * @param n The number of data points.
 * @param x The data set whose median is to be found.  Ideally, the
 *  data set should be monotonically increasing; however, if it is not, it
 *  may be sorted by the routine, dependent upon the value of @p srt.  On
 *  output, the array contents are unchanged; however, they may be sorted
 *  into ascending order (dependent upon the value of @p srt).
 * @param srt A logical flag determining if @p x should be sorted.
 *
 * @return The median of @p x.
 */
double median(int n, double *x, bool srt);

/** @brief Computes the sample variance of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @return The variance of @p x.
 *
 * @par Remarks
 * To avoid overflow-type issues, Welford's algorithm is employed.  A simple
 * illustration of this algorithm can be found 
 * [here](https://www.johndcook.com/blog/standard_deviation/).
 */
double variance(int n, const double *x);

/** @brief Computes the covariance matrix of N data sets of M observations.
 *
 * @param m The number of observations.
 * @param n The number of data sets.
 * @param x The M-by-N matrix of data.
 * @param c The N-by-N matrix where the resulting covariance matrix
 *  will be written.
 * @param err The errorhandler object.  If no error handling is
 *  desired, simply pass NULL, and errors will be dealt with by the default
 *  internal error handler.  Possible errors that may be encountered are as
 *  follows.
 *  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
 *      available.
 */
void covariance(int m, int n, const double *x, double *c, errorhandler *err);

/** @brief Computes the corrected standard deviation of a data set.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 *
 * @return The standard deviation of @p x.
 */
double standard_deviation(int n, const double *x);

/** @brief Computes the confidence interval based upon a standard normal 
 * distribution.
 *
 * @param n The number of data points.
 * @param x An N-element array containing the data set.
 * @param alpha The confidence level.  This value must lie between
 * zero and one such that: 0 < alpha < 1.
 *
 * @return The confidence interval as the deviation from the mean.
 *
 * @par Remarks
 * The confidence interval, assuming a standard normal distribution, is
 * as follows: mu +/- z * s / sqrt(n), where mu = the mean, and s = the
 * standard deviation.  This routine computes the z * s / sqrt(n) portion 
 * leaving the computation of the mean to the user.
 */
double confidence_interval(int n, const double *x, double alpha);

/** @brief Applies a moving average to smooth a data set.
!!
!! @param[in] n The number of data points.
!! @param[in,out] x On input, the signal to smooth.  On output, the smoothed
!!  signal.
!! @param[in] npts The size of the averaging window.  This value must be
!!  at least 2, but no more than the number of elements in @p x.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_INVALID_INPUT_ERROR: Occurs if @p npts is less than 2, or greater
!!      than the length of @p x.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
!!      available.
 */
void moving_average(int n, double *x, int npts, errorhandler *err);

/** @brief Employs a least squares fit to determine the coefficient A in the
!! linear system: Y = A * X.
!!
!! @param[in] n The number of data points.
!! @param[in] x An N-element array containing the independent variable data.
!! @param[in,out] y An N-element array containing the dependent variable
!!  data corresponding to @p x.  On output, the contents of this array are
!!  overwritten as it is used for storage purposes by the algorithm.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if insufficient memory is available.
!!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are different sizes.
!!
!! @return The scalar coefficient A.
 */
double least_squares_fit_1var(int n, const double *x, double *y, 
                              errorhandler *err);

/** @brief Employs a least squares fit to determine the coefficient A in the
!! linear system: Y = A * X.
!!
!! @param[in] m The number of dependent variables.
!! @param[in] n The number of independent variables.
!! @param[in,out] x An N-by-NPTS matrix containing the P data points of the
!!  N independent variables.
!! @param[in] y An M-by-NPTS matrix containing the P data points of the M
!!  dependent variables.
!! @param[out] a The M-by-N matrix where the resulting coefficient matrix A
!!  will be written.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_ARRAY_SIZE_ERROR: Occurs if any of the matrix dimensions are not
!!      compatiable.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
!!      available.
!!
!! @par Remarks
!! The algorithm attempts to compute the coefficient matrix A as follows.
!! Y * X**T = A * X * X**T
!! Y * X**T * INV(X * X**T) = A
!! This does require that X * X**T does not result in a singular matrix.  To
!! handle the situation where X * X**T is singular, the Moore-Penrose
!! pseudo-inverse, computed by means of singular value decomposition, is
!! utilized to still arrive at a solution that, at minimum, has a minimum
!! Euclidean norm of its residual.
!! Let: PINV(X) = X**T * INV(X * X**T),
!! Then: A = Y * PINV(X)
 */
void least_squares_fit_nvar(int m, int n, int npts, double *x, const double *y,
                            double *a, errorhandler *err);

/** @brief Initializes a new c_lowess_smoothing object.
!!
!! @param[out] obj The c_lowess_smoothing object.
!! @param[in] n The number of data points.
!! @param[in] x An N-element array containing the x-coordinate data.  
!!  Ideally, the data set should be monotonically increasing; however, if 
!!  it is not, it may be sorted by the routine, dependent upon the value 
!!  of @p srt.
!! @param[in] y An N-element array containing the y-coordinate data.
!! @param[in] srt A logical flag determining if @p x should be sorted.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_ARRAY_SIZE_ERROR: Occurs if @p x and @p y are not the same size.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory
!!      available.
 */
void alloc_lowess(lowess_smoothing *obj, int n, const double *x, 
                  const double *y, bool srt, errorhandler *err);

/** @brief Frees resources held by a c_lowess_smoothing object.
!!
!! @param[in,out] obj The c_lowess_smoothing object.
 */
void free_lowess(lowess_smoothing *obj);

/** @brief Performs the actual smoothing operation.
!!
!! @param[in,out] obj The c_lowess_smoothing object.
!! @param[in] f Specifies the amount of smoothing.  More specifically, this
!! value is the fraction of points used to compute each value.  As this 
!! value increases, the output becomes smoother.  Choosing a value in the
!! range of 0.2 to 0.8 usually results in a good fit.  As such, a reasonable
!! starting point, in the absence of better information, is a value of 0.5.
!! @param[in] n The size of the buffer @p y.  Ideally, this parameter is
!!  equal to the number of points stored in @p obj; however, the routine
!!  will only traverse the minimum of the this parameter or the number of
!!  points stored in @p obj.
!! @param[out] y An N-element array to which the smoothed data will be
!!  written.
!! @param[in,out] err The errorhandler object.  If no error handling is
!!  desired, simply pass NULL, and errors will be dealt with by the default
!!  internal error handler.  Possible errors that may be encountered are as
!!  follows.
!!  - CF_NO_DATA_DEFINED_ERROR: Occurs if no data has been defined.
!!  - CF_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
!!      available.
 */
void lowess_smooth(lowess_smoothing *obj, double f, int n, double *y, 
                   errorhandler *err);

/** @brief Gets the number of points used by the lowess_smoothing object.
!!
!! @param[in] obj The c_lowess_smoothing object.
!!
!! @return The number of points.
 */
int lowess_get_point_count(const lowess_smoothing *obj);

/** @brief Gets a copy of the data points stored by the lowess_smoothing
!! object.
!!
!! @param[in] obj The c_lowess_smoothing object.
!! @param[in] n The size of the buffer arrays.
!! @param[out] x An N-element array where the x-coordinate data will be 
!!  written.
!! @param[out] y An N-element array where the y-coordinate data will be 
!!  written.
!!
!! @par Remarks
!! If @p n is different than the actual number of points that exist, the 
!! lesser of the two values will be utilized.  The lowess_smoothing object
!! can be queried to determine the quantity of stored points.
 */
void lowess_get_points(const lowess_smoothing *obj, int n, double *x, 
                       double *y);

/** @brief Gets the residuals from each data point.
!!
!! @param[in] this The c_lowess_smoothing object.
!! @param[in] n The number of elements available in the buffer array @p x.
!! @param[out] x An N-element array where the residual data should be 
!!  written.
 */
void lowess_get_residuals(const lowess_smoothing *obj, int n, double *x);

#ifdef __cplusplus
}
#endif
#endif // CURVEFIT_H_DEFINED
