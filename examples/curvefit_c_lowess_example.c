// curvefit_c_lowess_example.c

#include "curvefit.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main() {
    // Constants
    const int n = 100;
    const double maxX = 1.0;
    const double minX = 0.0;

    // Local Variables
    int i, id;
    double x[n], y[n], yr[n], ys[n], ys2[n], dx, r;
    lowess_smoothing obj;
    FILE *f;
    time_t t;

    // Initialization
    srand((unsigned int)time(&t));
    dx = (maxX - minX) / (n - 1.0);
    r = ((double)rand()) / RAND_MAX - 0.5;
    x[0] = minX;
    y[0] = 0.5 * sin(20.0 * x[0]) + cos(5.0 * x[0]) * exp(-0.1 * x[0]);
    yr[0] = y[0] + r;
    for (i = 1; i < n; ++i) {
        x[i] = x[i-1] + dx;
        y[i] = 0.5 * sin(20.0 * x[i]) + cos(5.0 * x[i]) * exp(-0.1 * x[i]);

        r = ((double)rand()) / RAND_MAX - 0.5;
        yr[i] = y[i] + r;
    }
    alloc_lowess(&obj, n, x, yr, true, NULL);

    // Compute the fit for both f = 0.2, and f = 0.8
    lowess_smooth(&obj, 0.2, n, ys, NULL);
    lowess_smooth(&obj, 0.8, n, ys2, NULL);

    // Write the results to file
    f = fopen("c_lowess.txt", "w");
    for (i = 0; i < n; ++i) {
        fprintf(f, "%f,%f,%f,%f,%f\n", x[i], y[i], yr[i], ys[i], ys2[i]);
    }
    fclose(f);

    // End
    return 0;
}