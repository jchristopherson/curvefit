// curvefit.h
#ifndef CURVEFIT_H_DEFINED
#define CURVEFIT_H_DEFINED

#include <stdbool.h>
#include "ferror.h"
#include "nonlin.h"

/** @brief An error flag denoting an improperly sized array. */
#define CF_ARRAY_SIZE_ERROR = NL_ARRAY_SIZE_ERROR;
/** @brief An error denoting that there is insufficient memory available. */
#define CF_OUT_OF_MEMORY_ERROR = NL_OUT_OF_MEMORY_ERROR;
/** @brief An error denoting that no data has been defined. */
#define CF_NO_DATA_DEFINED_ERROR = 300;
/** @brief An error flag denoting an invalid input. */
#define CF_INVALID_INPUT_ERROR = NL_INVALID_INPUT_ERROR;
/** @brief An error flag denoting a non-monotonic array was given when a
!! monotonic array was expected. */
#define CF_NONMONOTONIC_ARRAY_ERROR = 301;
/** @brief An error resulting from an invalid operation. */
#define CF_INVALID_OPERATION_ERROR = NL_INVALID_OPERATION_ERROR;
/** @brief An error resulting from a lack of convergence. */
#define CF_CONVERGENCE_ERROR = NL_CONVERGENCE_ERROR;
/** @brief An error indicating the user-requested tolerance is too small to be
!! practical for the problem at hand. */
#define CF_TOLERANCE_TOO_SMALL_ERROR = NL_TOLERANCE_TOO_SMALL_ERROR;

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
!! @param[in] obj The c_linear_interp object.
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

#ifdef __cplusplus
}
#endif
#endif // CURVEFIT_H_DEFINED
