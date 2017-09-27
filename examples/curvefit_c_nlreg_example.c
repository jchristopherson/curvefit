// curvefit_c_nlreg_example.c

#include <stdio.h>
#include "curvefit.h"

#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
double cfcn(double x, int n, const double *c);

int main() {
    // Local Variables
    const int n = 21;
    const int ncoeff = 4;
    nonlinear_regression obj;
    double 
        cp[4] = {1.0, 1.0, 1.0, 1.0},
        xp[21] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 
            1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0},
        yp[21] = {1.216737514, 1.250032542, 1.305579195, 1.040182335,
            1.751867738, 1.109716707, 2.018141531, 1.992418729,
            1.807916923, 2.078806005, 2.698801324, 2.644662712,
            3.412756702, 4.406137221, 4.567156645, 4.999550779,
            5.652854194, 6.784320119, 8.307936836, 8.395126494,
            10.30252404};
    iteration_behavior ib;
    
    // Initialization
    alloc_nonlinear_regression(&obj, n, xp, yp, cfcn, ncoeff, NULL);

    // Solve
    nonlinreg_solve(&obj, ncoeff, cp, &ib, NULL);

    // Print the results
    printf("c0: %f\nc1: %f\nc2: %f\nc3: %f\n", cp[3], cp[2], cp[1], cp[0]);

    // End
    free_nonlinear_regression(&obj);
    return 0;
}

double cfcn(double x, int n, const double *c) {
    return c[0] * CUBE(x) + c[1] * SQR(x) + c[2] * x + c[3];
}
