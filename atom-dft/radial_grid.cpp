//
//  radial_grid.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#include "radial_grid.hpp"

#include <cmath>
#include <iostream>

void radial_grid::set_exp_grid(double rc, double r0) {
    this->r0 = r0;
    for (int i=0; i<N; i++) {
        u[i] = double(i) / double(N-1) * log(rc / r0);
        r[i] = r0 * exp(u[i]);
        drdu[i] = r[i];
    }
    u_step = u[1];
}

double radial_grid::get_r_by_u_exp_grid(double u) {
    return r0 * exp(u);
}
