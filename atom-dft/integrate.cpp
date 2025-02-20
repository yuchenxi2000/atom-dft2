//
//  integrate.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#include "integrate.hpp"

#include <iostream>

double simpson(int size, double * y, double * x) {
    if (size < 3) {
        std::cerr << "simpson rule requires at least 3 points, but " << size << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    double I;
    int idx_start;
    if (size % 2 == 0) {
        double c1 = (x[1]-x[0]) / (x[2]-x[0]);
        double c2 = (x[1]-x[0]) * (x[1]-x[0]) / (x[2]-x[1]) / (x[2]-x[0]);
        I = (x[1]-x[0])/6.0 * ((3.0 - c1) * y[0] + (3.0 + c2 + c1) * y[1] - c2 * y[2]);
        idx_start = 3;
    } else {
        I = 0.0;
        idx_start = 2;
    }
    for (int i=idx_start; i<size; i+=2) {
        double c1 = (x[i]-x[i-1]) / (x[i-1]-x[i-2]);
        double c2 = (x[i]-x[i-2]) * (x[i]-x[i-2]) / (x[i]-x[i-1]) / (x[i-1]-x[i-2]);
        double c3 = (x[i-1]-x[i-2]) / (x[i]-x[i-1]);
        I += (x[i]-x[i-2])/6.0 * ((2.0 - c1) * y[i-2] + c2 * y[i-1] + (2.0 - c3) * y[i]);
    }
    return I;
}

double simpson_uniform(int size, double * y, double step) {
    if (size < 3) {
        std::cerr << "simpson rule requires at least 3 points, but " << size << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    double I;
    int idx_start;
    if (size % 2 == 0) {
        I = step / 3.0 * (5.0/4.0 * y[0] + 2.0 * y[1] - 1.0/4.0 * y[2]);
        idx_start = 3;
    } else {
        I = 0.0;
        idx_start = 2;
    }
    for (int i=idx_start; i<size; i+=2) {
        I += step / 3.0 * (y[i-2] + 4.0 * y[i-1] + y[i]);
    }
    return I;
}

double trapz_7_uniform(int size, double * y, double step) {
    double s = (
                36799 * (y[0] + y[size-1])
                + 176648 * (y[1] + y[size-2])
                + 54851 * (y[2] + y[size-3])
                + 177984 * (y[3] + y[size-4])
                + 89437 * (y[4] + y[size-5])
                + 130936 * (y[5] + y[size-6])
                + 119585 * (y[6] + y[size-7])
                ) / 120960;
    for (int i=7; i<size-7; i++) {
        s += y[i];
    }
    return step * s;
}

void simpson_cumsum(int size, double * y, double * x, double * s) {
    if (size < 3) {
        std::cerr << "simpson rule requires at least 3 points, but " << size << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    double c1 = (x[1]-x[0]) / (x[2]-x[0]);
    double c2 = (x[1]-x[0]) * (x[1]-x[0]) / (x[2]-x[1]) / (x[2]-x[0]);
    s[1] = s[0] + (x[1]-x[0])/6.0 * ((3.0 - c1) * y[0] + (3.0 + c2 + c1) * y[1] - c2 * y[2]);
    for (int i=2; i<size; i++) {
        double c1 = (x[i]-x[i-1]) / (x[i-1]-x[i-2]);
        double c2 = (x[i]-x[i-2]) * (x[i]-x[i-2]) / (x[i]-x[i-1]) / (x[i-1]-x[i-2]);
        double c3 = (x[i-1]-x[i-2]) / (x[i]-x[i-1]);
        s[i] = s[i-2] + (x[i]-x[i-2])/6.0 * ((2.0 - c1) * y[i-2] + c2 * y[i-1] + (2.0 - c3) * y[i]);
    }
}
