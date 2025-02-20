//
//  poisson.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#include "poisson.hpp"

#include "radial_grid.hpp"
#include "ode.hpp"
#include "interp.hpp"
#include "integrate.hpp"
#include <cmath>

extern "C" void solve_poisson(double * VH, double * dVH, double * density_r, double r0, double rc, int N) {
    radial_grid grid(N);
    grid.set_exp_grid(rc, r0);
    
    ode_2d_func func = [density_r, grid] (int i, double y1, double y2, double * dy1, double * dy2) {
        *dy1 = grid.r[i] * y2;
        *dy2 = -2.0 * y2 - 4.0 * M_PI * density_r[i];
    };
    double density_r_mid[3];
    double r_mid[3];
    for (int i=0; i<3; i++) {
        density_r_mid[i] = interp_poly(grid.u[i] + 0.5 * grid.u_step, 4, density_r, grid.u);
        r_mid[i] = interp_poly(grid.u[i] + 0.5 * grid.u_step, 4, grid.r, grid.u);
    }
    ode_2d_func func_mid = [density_r_mid, r_mid] (int i, double y1, double y2, double * dy1, double * dy2) {
        *dy1 = r_mid[i] * y2;
        *dy2 = -2.0 * y2 - 4.0 * M_PI * density_r_mid[i];
    };
    // we use 4th-order runge-kutta for the first 4 points
    double * density_r2 = new double[N];
    for (int i=0; i<N; i++) {
        density_r2[i] = density_r[i] * grid.r[i];
    }
    VH[0] = 4.0 * M_PI * trapz_7_uniform(N, density_r2, grid.u_step);
    dVH[0] = 0.0;
    runge_kutta_4_outward(4, VH, dVH, func, func_mid, grid.u_step);
    // 4th-order explicit adam for the remaining points
    adam_4_outward(grid.N, VH, dVH, func, grid.u_step);
    delete [] density_r2;
}
