// curvefit_c_test.c

#include <stdio.h>
#include "curvefit.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "c_test_core.h"

bool test_mean();
bool test_variance();
bool test_confidence_interval();
bool test_median();
bool test_covariance();

int main() {
    // Local Variables
    bool rst, overall;

    // Initialization
    overall = true;

    // Statistics Tests
    rst = test_mean();
    if (!rst) overall = false;

    rst = test_variance();
    if (!rst) overall = false;

    rst = test_confidence_interval();
    if (!rst) overall = false;

    rst = test_median();
    if (!rst) overall = false;

    rst = test_covariance();
    if (!rst) overall = false;

    // End
    if (overall)
        printf("C API CURVEFIT TEST STATUS: PASS\n");
    else
        printf("C API CURVEFIT TEST STATUS: FAILED\n");
    return 0;
}


/* ************************************************************************** */
bool test_mean() {
    // Constants
    const double ans = 10.0;
    const double tol = 1.0e-8;
    const int npts = 5;

    // Local Variables
    bool rst;
    double mu, x[5] = {9.0, 10.0, 11.0, 7.0, 13.0};

    // Process
    rst = true;
    mu = mean(npts, x);
    if (fabs(mu - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a mean value of %f, but computed a value of %f.\n",
            ans, mu);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_variance() {
    // Constants
    const double ans = 5.0;
    const double tol = 1.0e-8;
    const int npts = 5;

    // Local Variables
    bool rst;
    double v, x[5] = {9.0, 10.0, 11.0, 7.0, 13.0};

    // Process
    rst = true;
    v = variance(npts, x);
    if (fabs(v - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a variance of %f, but computed a value of %f.\n",
            ans, v);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_confidence_interval() {
    // Constants
    const int n = 1000;
    const double alpha = 0.95;
    const double z = 1.96;
    const double tol = 1.0e-3;  // Due to the precision of Z

    // Local Variables
    bool rst;
    double x[n], c, ans;
    int i;
    time_t t;

    // Initialization
    rst = true;
    srand((unsigned int)time(&t));
    for (i = 0; i < n; ++i) x[i] = ((double)rand()) / RAND_MAX;
    ans = z * standard_deviation(n, x) / sqrt((double)n);

    // Process
    c = confidence_interval(n, x, alpha);
    if (fabs(c - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a confidence interval of %f, but computed a value of %f\n", 
            ans, c);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_median() {
    // Constants
    const int n = 7;
    const double tol = 1.0e-8;
    const double ans = 1.0;

    // Local Variables
    bool rst;
    double m, x[7] = {2.3, 3.0, 1.0, -4.0, -7.0, 1.0, 3.0};

    // Process
    rst = true;
    m = median(n, x, true);
    if (fabs(m - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a median of %f, but computed a value of %f\n", 
            ans, m);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_covariance() {
    // Constants
    const int m = 3;
    const int n = 4;
    const double tol = 1.0e-8;

    // Local Variables
    bool rst;
    double c[n*n], 
        x[12] = {5.0, 1.0, 4.0, 0.0, -5.0, 9.0, 3.0, 7.0, 8.0, 7.0, 3.0, 1.0e1},
        ans[16] = {4.33333333333333, 8.83333333333333, -3.0,
            5.66666666666667, 8.83333333333333, 50.3333333333333,
            6.5, 24.1666666666667, -3.0, 6.5, 7.0, 1.0,
            5.66666666666667, 24.1666666666667, 1.0, 12.3333333333333};
    
    // Process
    rst = true;
    covariance(m, n, x, c, NULL);
    if (!is_dbl_mtx_equal(n, n, c, ans, tol)) {
        rst = false;
        printf("Test Failed: Covariance matrix of 4 data sets.\n");
    }

    // End
    return rst;
}

/* ************************************************************************** */

/* ************************************************************************** */