// curvefit_c_interp_example.c

#include "curvefit.h"
#include <stdio.h>

int main() {
    // Local Variables
    const int npts = 100;
    int i;
    double dx, xi[npts], yi[npts],
        x[9] = {-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0},
        y[9] = {0.0, 0.15, 1.12, 2.36, 2.36, 1.46, 0.49, 0.06, 0.0};
    linear_interp obj;
    FILE *f;
    
    // Initialization
    dx = 8.0 / (npts - 1.0);
    xi[0] = -4.0;
    for (i = 1; i < npts; ++i) xi[i] = xi[i-1] + dx;
    alloc_linear_interp(&obj, 9, x, y, NULL);

    // Interpolate
    linear_interpolate(&obj, npts, xi, yi);

    // Write the results to file
    f = fopen("c_linear_interp.txt", "w");
    for (i = 0; i < npts; ++i) {
        if (i < 9)
            fprintf(f, "%f,%f,%f,%f\n", x[i], y[i], xi[i], yi[i]);
        else
            fprintf(f, ",,%f,%f\n", xi[i], yi[i]);
    }
    fclose(f);

    // End
    free_linear_interp(&obj);
    return 0;
}