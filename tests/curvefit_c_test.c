// curvefit_c_test.c

#include <stdio.h>
#include "curvefit.h"
#include "linalg.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "c_test_core.h"


bool test_mean();
bool test_variance();
bool test_confidence_interval();
bool test_median();
bool test_covariance();
bool test_seb();
bool test_nonlin();
bool test_term_nonlin();
bool test_hysteresis();
bool test_return_to_zero();
bool test_repeatability();
bool test_crosstalk();

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

    // Calibration Tests
    rst = test_seb();
    if (!rst) overall = false;

    rst = test_nonlin();
    if (!rst) overall = false;

    rst = test_term_nonlin();
    if (!rst) overall = false;

    rst = test_hysteresis();
    if (!rst) overall = false;

    rst = test_return_to_zero();
    if (!rst) overall = false;

    rst = test_repeatability();
    if (!rst) overall = false;

    rst = test_crosstalk();
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
bool test_seb() {
    // Constants
    const double fullscale = 5.0e2;
    const double sebout_ans = 2.8001;
    const double sebpfs_ans = 0.0104;
    const double tol = 1.0e-4;
    const int n = 11;

    // Local Variables
    bool rst;
    seb_results s;
    double spfs,
        applied[] = {0.0, 1.0e2, 2.0e2, 3.0e2, 4.0e2, 5.0e2, 4.0e2, 3.0e2, 
            2.0e2, 1.0e2, 0.0},
        output[] = {0.0, 0.55983, 1.11975, 1.67982, 2.24005, 2.80039, 
            2.24023, 1.68021, 1.12026, 0.56021, 0.00006};
    
    // Test
    rst = true;
    seb(n, applied, output, fullscale, &s, NULL);
    spfs = 1.0e2 * s.seb / fullscale;
    if (fabs(spfs - sebpfs_ans) > tol) {
        rst = false;
        printf("Test Failed: Expected an SEB of: %f, but computed: %f.\n",
            sebpfs_ans, spfs);
    }
    if (fabs(s.output - sebout_ans) > tol) {
        rst = false;
        printf("Test Failed: Expected an SEB output of: %f, but computed: %f.\n",
            sebout_ans, s.output);
    }

    // End
    return rst;
}


/* ************************************************************************** */
bool test_nonlin() {
    // Constants
    const double ans = -0.062122611;
    const double slope = 178.5558182;
    const double tol = 1.0e-4;
    const int n = 11;

    // Local Variables
    bool rst;
    int i;
    double x, measured[n],
        applied[] = {0.0, 1.0e2, 2.0e2, 3.0e2, 4.0e2, 5.0e2, 4.0e2, 3.0e2, 
            2.0e2, 1.0e2, 0.0},
        output[] = {0.0, 0.55983, 1.11975, 1.67982, 2.24005, 2.80039, 
            2.24023, 1.68021, 1.12026, 0.56021, 0.00006};
    
    // Apply the calibration gain
    for (i = 0; i < n; ++i) measured[i] = slope * output[i];

    // Test
    x = nonlinearity(n, applied, measured);
    rst = true;
    if (fabs(x - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a nonlinearity error of: %f, but computed: %f.\n",
            ans, x);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_term_nonlin() {
    // Constants
    const double ans = -0.073918276;
    const double slope = 178.5558182;
    const double tol = 1.0e-4;
    const int n = 11;

    // Local Variables
    bool rst;
    int i;
    double x, measured[n],
        applied[] = {0.0, 1.0e2, 2.0e2, 3.0e2, 4.0e2, 5.0e2, 4.0e2, 3.0e2, 
            2.0e2, 1.0e2, 0.0},
        output[] = {0.0, 0.55983, 1.11975, 1.67982, 2.24005, 2.80039, 
            2.24023, 1.68021, 1.12026, 0.56021, 0.00006};
    
    // Apply the calibration gain
    for (i = 0; i < n; ++i) measured[i] = slope * output[i];

    // Test
    x = terminal_nonlinearity(n, applied, measured);
    rst = true;
    if (fabs(x - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a terminal nonlinearity error of: %f, but computed: %f.\n",
            ans, x);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_hysteresis() {
    // Constants
    const double ans = 0.091063467;
    const double slope = 178.5558182;
    const double tol = 1.0e-4;
    const int n = 11;

    // Local Variables
    bool rst;
    int i;
    double x, measured[n],
        applied[] = {0.0, 1.0e2, 2.0e2, 3.0e2, 4.0e2, 5.0e2, 4.0e2, 3.0e2, 
            2.0e2, 1.0e2, 0.0},
        output[] = {0.0, 0.55983, 1.11975, 1.67982, 2.24005, 2.80039, 
            2.24023, 1.68021, 1.12026, 0.56021, 0.00006};
    
    // Apply the calibration gain
    for (i = 0; i < n; ++i) measured[i] = slope * output[i];

    // Test
    x = hysteresis(n, applied, measured);
    rst = true;
    if (fabs(x - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a hysteresis error of: %f, but computed: %f.\n",
            ans, x);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_return_to_zero() {
    // Constants
    const double ans = 0.010713349;
    const double slope = 178.5558182;
    const double tol = 1.0e-4;
    const int n = 11;

    // Local Variables
    bool rst;
    int i;
    double x, measured[n],
        applied[] = {0.0, 1.0e2, 2.0e2, 3.0e2, 4.0e2, 5.0e2, 4.0e2, 3.0e2, 
            2.0e2, 1.0e2, 0.0},
        output[] = {0.0, 0.55983, 1.11975, 1.67982, 2.24005, 2.80039, 
            2.24023, 1.68021, 1.12026, 0.56021, 0.00006};
    
    // Apply the calibration gain
    for (i = 0; i < n; ++i) measured[i] = slope * output[i];

    // Test
    x = return_to_zero(n, applied, measured, 1e-4);
    rst = true;
    if (fabs(x - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a return to zero error of: %f, but computed: %f.\n",
            ans, x);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_repeatability() {
    // Constants
    const int npts = 15;
    const int ntests = 3;
    const double slope = 2736.031907;
    const double ans = 0.30096351;
    const double tol = 1.0e-4;

    // Local Variables
    bool rst;
    int i;
    double x, measured[npts * ntests],
        applied[] = {0.0, 100.0, 300.0, 600.0, 900.0,
            1200.0, 1500.0, 1800.0, 2100.0, 2400.0, 2700.0,
            3000.0, 3400.0, 1500.0, 0.0, 0.0, 100.0, 300.0,
            600.0, 900.0, 1200.0, 1500.0, 1800.0, 2100.0,
            2400.0, 2700.0, 3000.0, 3400.0, 1500.0, 0.0, 0.0,
            100.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0,
            2100.0, 2400.0, 2700.0, 3000.0, 3400.0, 1500.0,
            0.0},
        output[] = {0.0, 0.03654, 0.10961, 0.2192, 0.32883,
            0.43845, 0.5481, 0.65775, 0.76743, 0.87711, 0.9868,
            1.09646, 1.2427, 0.5481, 0.0, 0.0, 0.03655, 0.10965,
            0.21925, 0.32888, 0.43852, 0.54818, 0.65784, 0.76752,
            0.8772, 0.98689, 1.09656, 1.24281, 0.54814, 0.00001,
            0.0, 0.03655, 0.10961, 0.21921, 0.32884, 0.43847,
            0.54813, 0.65779, 0.76746, 0.87715, 0.98682, 1.09648,
            1.24274, 0.54813, 0.00001};
    
    // Apply the calibration
    for (i = 0; i < npts * ntests; ++i) measured[i] = slope * output[i];

    // Test
    x = repeatability(npts, ntests, applied, measured);
    rst = true;
    if (fabs(x - ans) > tol) {
        rst = false;
        printf("Test Failed: Expected a repeatability error of: %f, but computed: %f.\n",
            ans, x);
    }

    // End
    return rst;
}

/* ************************************************************************** */
bool test_crosstalk() {
    // Constants
    const int npts = 34;
    const int ndof = 2;
    const double tol = 1.0e-4;

    // Local Variables
    bool rst;
    int i, j, indices[] = {1, 17, 18, 34};
    double xerr[npts*ndof], xmeas[npts*ndof], xint[ndof*npts], 
        xmeast[ndof*npts], c[ndof*ndof], ans[ndof*ndof], xt[ndof*ndof],
        xin[] = {0.0, 3000.0, 6000.0, 7500.0, 9000.0, 12000.0, 
            15000.0, 7500.0, 0.0, 0.0, -3000.0, -6000.0, -7500.0, -9000.0, 
            -12000.0, -15000.0, -7500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
            0.0, 0.0, 0.0, 67.7908728, 135.5817456, 203.3726184, 271.1634912, 
            338.954364, 203.3726184, 0.0, 0.0, -67.7908728, -135.5817456, 
            -203.3726184, -271.1634912, -338.954364, -203.3726184, 0.0},
        xout[] = {0.0, 0.38905, 0.77816, 0.97269, 1.16714, 1.556, 
            1.94484, 0.9726, -1.0e-5, 0.0, -0.388886, -0.77775, -0.97215, 
            -1.16654, -1.55533, -1.9441, -0.97171, 4.0e-5, 0.0, -0.00044, 
            -0.0013, -0.0024, -0.00382, -0.00528, -0.00257, 0.00015, 0.0, 
            0.00144, 0.00306, 0.00446, 0.00567, 0.00688, 0.00451, -2.0e-5, 
            0.0, 0.00122, 0.00259, 0.0029, 0.00314, 0.00338, 0.00356, 0.00477, 
            -1.0e-5, 0.0, 0.00021, 0.00051, 0.00069, 0.00088, 0.0013, 0.00178, 
            0.00058, 3.0e-5, 0.0, 0.27156, 0.54329, 0.81507, 1.08682, 1.35881, 
            0.81553, 1.0e-5, 0.0, -0.27145, -0.54312, -0.81493, -1.0868, 
            -1.35879, -0.81548, 0.0};
    
    // Transpose the input and output matrices
    for (j = 0; j < ndof; ++j) {
        for (i = 0; i < npts; ++i) {
            if (i != j) {
                xint[INDEX(j,i,ndof)] = xin[INDEX(i,j,npts)];
                xmeast[INDEX(j,i,ndof)] = xout[INDEX(i,j,npts)];
            }
        }
    }

    // Compute the calibration gains, and then apply the results
    least_squares_fit_nvar(ndof, ndof, npts, xmeast, xint, c, NULL);
    mtx_mult_(false, true, npts, ndof, ndof, 1.0, xout, npts, c, ndof, 0.0,
        xmeas);

    // Compute the measurement errors
    for (i = 0; i < ndof * npts; ++i) xerr[i] = xmeas[i] - xin[i];

    // Compute the crosstalk errors
    crosstalk(npts, ndof, xerr, indices, xt, NULL);

    // The solution is:
    ans[0] = 0.0;
    ans[1] = nonlinearity(17, &xin[INDEX(17,0,npts)], &xmeas[INDEX(17,0,npts)]);

    ans[2] = nonlinearity(17, &xin[INDEX(0,1,npts)], &xmeas[INDEX(0,1,npts)]);
    ans[3] = 0.0;

    // Test
    rst = true;
    if (!is_dbl_mtx_equal(ndof, ndof, xt, ans, tol)) {
        rst = false;
        printf("Test Failed: The crosstalk error test failed.\n");
    }

    // End
    return rst;
}

/* ************************************************************************** */