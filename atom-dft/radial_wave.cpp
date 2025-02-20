//
//  radial_wave.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//
//  Some of the eigenvalue finding code is rewritten from Fortran code in project https://github.com/certik/dftatom
//  which is licensed under MIT license. Copyright (c) 2006-2012 Ondřej Čertík, Jiří Vackář, John Pask
//

#include "radial_wave.hpp"

#include "ode.hpp"
#include "integrate.hpp"
#include "radial_grid.hpp"
#include <cmath>

extern "C" int solve_radial_from_zero(double * chi, double * dchi, int * num_nodes, double E, int l, double Z, double * Vext, double r0, double rc, int N) {
    radial_grid grid(N);
    grid.set_exp_grid(rc, r0);
    
    double * p2r = new double[N];
    for (int i=0; i<N; i++) {
        p2r[i] = 2.0 * ((E - Vext[i]) * grid.r[i] + Z - 0.5 * l * (l+1) / grid.r[i]);
    }
    
    ode_2d_linear_coeff coeff_func = [grid, p2r] (int i, double * mat, double * vec) {
        mat[0] = 0.0;
        mat[1] = grid.r[i];
        mat[2] = -p2r[i];
        mat[3] = 0.0;
        vec[0] = 0.0;
        vec[1] = 0.0;
    };
    
    for (int i=0; i<3; i++) {
        chi[i] = pow(grid.r[i], l+1);
        dchi[i] = (l + 1) * pow(grid.r[i], l);
    }
    implicit_adams_4_linear_outward(N, chi, dchi, coeff_func, grid.u_step);
    
    // count number of nodes
    *num_nodes = 0;
    for (int i=0; i<N-1; i++) {
        if (abs(chi[i]) > 1e10) {
            break;
        } else if (chi[i] == 0.0) {
            *num_nodes += 1;
        } else if ((chi[i] > 0.0 && chi[i+1] < 0.0) || (chi[i] < 0.0 && chi[i+1] > 0.0)) {
            *num_nodes += 1;
        }
    }
    
    delete [] p2r;
    
    return 0;
}

extern "C" int solve_radial_from_zero_and_inf(double * chi, double * dchi, double * dE, int * num_nodes, double E, int l, double Z, double * Vext, double r0, double rc, int N) {
    
    if (E >= 0.0) {
        return 1;
    }
    
    radial_grid grid(N);
    grid.set_exp_grid(rc, r0);
    
    double * p2r = new double[N];
    for (int i=0; i<N; i++) {
        p2r[i] = 2.0 * ((E - Vext[i]) * grid.r[i] + Z - 0.5 * l * (l+1) / grid.r[i]);
    }
    
    // find index of last zero point of p2r
    int idx_last_zero;
    
    for (idx_last_zero = grid.N-1; idx_last_zero>0; idx_last_zero--) {
        if (p2r[idx_last_zero] == 0.0) {
            break;
        } else if ((p2r[idx_last_zero] < 0.0 && p2r[idx_last_zero-1] > 0.0) || (p2r[idx_last_zero] > 0.0 && p2r[idx_last_zero-1] < 0.0)) {
            break;
        }
    }
    
    if (idx_last_zero == grid.N-1 || idx_last_zero < 10) {
        return 2;
    }
    
    // find index of inf
    int idx_inf;
    const double pseudo_phase_min = 36.0;
    for (idx_inf = idx_last_zero; idx_inf<N; idx_inf++) {
        double pseudo_phase = sqrt(-2.0 * E) * (grid.r[idx_inf] - grid.r[idx_last_zero]);
        if (pseudo_phase >= pseudo_phase_min) {
            break;
        }
    }
    
    if (idx_inf == grid.N-1) {
        return 3;
    }

    if (idx_inf - idx_last_zero < 10) {
        return 4;
    }

    if (grid.r[idx_last_zero] / grid.r[idx_inf] > 0.5) {
        return 5;
    }
    
    ode_2d_linear_coeff coeff_func = [grid, p2r] (int i, double * mat, double * vec) {
        mat[0] = 0.0;
        mat[1] = grid.r[i];
        mat[2] = -p2r[i];
        mat[3] = 0.0;
        vec[0] = 0.0;
        vec[1] = 0.0;
    };
    
    ode_2d_linear_coeff coeff_func_reverse = [grid, p2r, idx_last_zero] (int i, double * mat, double * vec) {
        mat[0] = 0.0;
        mat[1] = grid.r[i + idx_last_zero];
        mat[2] = -p2r[i + idx_last_zero];
        mat[3] = 0.0;
        vec[0] = 0.0;
        vec[1] = 0.0;
    };
    
    // solve from zero
    for (int i=0; i<3; i++) {
        chi[i] = pow(grid.r[i], l+1);
        dchi[i] = (l + 1) * pow(grid.r[i], l);
    }
    implicit_adams_4_linear_outward(idx_last_zero + 1, chi, dchi, coeff_func, grid.u_step);
    
    // backup
    double chi_last_zero = chi[idx_last_zero];
    double dchi_last_zero = dchi[idx_last_zero];
    
    // solve from inf
    for (int i=idx_inf-3; i<idx_inf; i++) {
        chi[i] = exp(-sqrt(-2.0 * E) * (grid.r[i] - grid.r[idx_last_zero]));
        dchi[i] = -sqrt(-2.0 * E) * chi[i];
    }
    implicit_adams_4_linear_inward(idx_inf - idx_last_zero + 1, chi, dchi, coeff_func_reverse, grid.u_step);
    
    // match two solutions
    double factor = chi_last_zero / chi[idx_last_zero];
    for (int i=idx_last_zero; i<idx_inf+1; i++) {
        chi[i] *= factor;
        dchi[i] *= factor;
    }
    
    // count number of nodes
    *num_nodes = 0;
    for (int i=0; i<idx_inf; i++) {
        if (abs(chi[i]) > 1e10) {
            break;
        } else if (chi[i] == 0.0) {
            *num_nodes += 1;
        } else if ((chi[i] > 0.0 && chi[i+1] < 0.0) || (chi[i] < 0.0 && chi[i+1] > 0.0)) {
            *num_nodes += 1;
        }
    }
    
    // calculate dE
    double * chi2r = new double[N];
    for (int i=0; i<N; i++) {
        chi2r[i] = chi[i] * chi[i] * grid.r[i];
    }
    double S = trapz_7_uniform(idx_inf + 1, chi2r, grid.u_step);
    *dE = chi_last_zero * (dchi_last_zero - dchi[idx_last_zero]) / (2.0 * S);
    
    delete [] p2r;
    delete [] chi2r;
    
    return 0;
}
