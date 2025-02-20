//
//  interp.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/12.
//

#include "interp.hpp"

#include <iostream>

double interp_poly(double z, int n, double * y, double * x) {
    double s = 0.0;
    for (int i=0; i<n; i++) {
        double prod = 1.0;
        for (int j=0; j<n; j++) {
            if (j != i) {
                prod *= (z - x[j]) / (x[i] - x[j]);
            }
        }
        s += y[i] * prod;
    }
    return s;
}

void get_mid_array_poly(double * mid, int size, double * y, double * x, int ninterp) {
    if (size < ninterp) {
        std::cerr << "array size should be larger than number of interpolation points!" << std::endl;
        exit(EXIT_FAILURE);
    }
    int ninterp_l = ninterp / 2;
    int ninterp_r = ninterp - ninterp_l;
    int start_idx;
    for (int i=0; i<size-1; i++) {
        if (i < ninterp_l) {
            start_idx = 0;
        } else if (i > size - ninterp_r) {
            start_idx = size - ninterp;
        } else {
            start_idx = i - ninterp_l;
        }
        double x_mid = 0.5 * (x[i] + x[i+1]);
        mid[i] = interp_poly(x_mid, ninterp, y+start_idx, x+start_idx);
    }
}
