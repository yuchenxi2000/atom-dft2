//
//  ode.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/18.
//

#include "ode.hpp"

#include <iostream>
#include "linalg.hpp"

void runge_kutta_4_outward(int nsteps, double * y1, double * y2, ode_2d_func func, ode_2d_func func_mid, double h_step) {
    double k_chi[4];
    double k_dchi[4];
    for (int i=0; i<nsteps-1; i++) {
        func(i, y1[i], y2[i], k_chi, k_dchi);
        func_mid(i, y1[i] + 0.5*h_step*k_chi[0], y2[i] + 0.5*h_step*k_dchi[0], k_chi+1, k_dchi+1);
        func_mid(i, y1[i] + 0.5*h_step*k_chi[1], y2[i] + 0.5*h_step*k_dchi[1], k_chi+2, k_dchi+2);
        func(i+1, y1[i] + h_step*k_chi[2], y2[i] + h_step*k_dchi[2], k_chi+3, k_dchi+3);
        y1[i+1] = y1[i] + h_step/6.0 * (k_chi[0] + 2.0 * k_chi[1] + 2.0 * k_chi[2] + k_chi[3]);
        y2[i+1] = y2[i] + h_step/6.0 * (k_dchi[0] + 2.0 * k_dchi[1] + 2.0 * k_dchi[2] + k_dchi[3]);
    }
}

void runge_kutta_4_inward(int nsteps, double * chi, double * dchi, ode_2d_func func, ode_2d_func func_mid, double h_step) {
    double k_chi[4];
    double k_dchi[4];
    for (int i=nsteps-1; i>0; i--) {
        func(i, chi[i], dchi[i], k_chi, k_dchi);
        func_mid(i-1, chi[i] - 0.5*h_step*k_chi[0], dchi[i] - 0.5*h_step*k_dchi[0], k_chi+1, k_dchi+1);
        func_mid(i-1, chi[i] - 0.5*h_step*k_chi[1], dchi[i] - 0.5*h_step*k_dchi[1], k_chi+2, k_dchi+2);
        func(i-1, chi[i] - h_step*k_chi[2], dchi[i] - h_step*k_dchi[2], k_chi+3, k_dchi+3);
        chi[i-1] = chi[i] - h_step/6.0 * (k_chi[0] + 2.0 * k_chi[1] + 2.0 * k_chi[2] + k_chi[3]);
        dchi[i-1] = dchi[i] - h_step/6.0 * (k_dchi[0] + 2.0 * k_dchi[1] + 2.0 * k_dchi[2] + k_dchi[3]);
    }
}

void adam_4_outward(int nsteps, double * y1, double * y2, ode_2d_func func, double h_step) {
    if (nsteps < 4) {
        std::cerr << "4th-order explicit adams for linear 2d ode requires at least 4 points, but " << nsteps << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    for (int i=0; i<nsteps-4; i++) {
        double dy[5][2];
        for (int j=0; j<4; j++) {
            func(i+j, y1[i+j], y2[i+j], dy[j], dy[j]+1);
        }
        double y1p, y2p;
        y1p = y1[i+3] + h_step/24.0 * (55.0 * dy[3][0] - 59.0 * dy[2][0] + 37.0 * dy[1][0] - 9.0 * dy[0][0]);
        y2p = y2[i+3] + h_step/24.0 * (55.0 * dy[3][1] - 59.0 * dy[2][1] + 37.0 * dy[1][1] - 9.0 * dy[0][1]);
        func(i+4, y1p, y2p, dy[4], dy[4]+1);
        y1[i+4] = y1[i+3] + h_step/24.0 * (9.0 * dy[4][0] + 19.0 * dy[3][0] - 5.0 * dy[2][0] + dy[1][0]);
        y2[i+4] = y2[i+3] + h_step/24.0 * (9.0 * dy[4][1] + 19.0 * dy[3][1] - 5.0 * dy[2][1] + dy[1][1]);
    }
}

void adam_4_inward(int nsteps, double * y1, double * y2, ode_2d_func func, double h_step) {
    if (nsteps < 4) {
        std::cerr << "4th-order explicit adams for linear 2d ode requires at least 4 points, but " << nsteps << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    for (int i=nsteps-1; i>3; i--) {
        double dy[5][2];
        for (int j=0; j<4; j++) {
            func(i-j, y1[i-j], y2[i-j], dy[j], dy[j]+1);
        }
        double y1p, y2p;
        y1p = y1[i-3] - h_step/24.0 * (55.0 * dy[3][0] - 59.0 * dy[2][0] + 37.0 * dy[1][0] - 9.0 * dy[0][0]);
        y2p = y2[i-3] - h_step/24.0 * (55.0 * dy[3][1] - 59.0 * dy[2][1] + 37.0 * dy[1][1] - 9.0 * dy[0][1]);
        func(i-4, y1p, y2p, dy[4], dy[4]+1);
        y1[i-4] = y1[i-3] - h_step/24.0 * (9.0 * dy[4][0] + 19.0 * dy[3][0] - 5.0 * dy[2][0] + dy[1][0]);
        y2[i-4] = y2[i-3] - h_step/24.0 * (9.0 * dy[4][1] + 19.0 * dy[3][1] - 5.0 * dy[2][1] + dy[1][1]);
    }
}

void implicit_adams_4_linear_outward(int nsteps, double * y1, double * y2, ode_2d_linear_coeff coeff_func, double h_step) {
    if (nsteps < 3) {
        std::cerr << "4th-order implicit adams for linear 2d ode requires at least 3 points, but " << nsteps << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    double mat[4][4];
    double vec[4][2];
    for (int i=0; i<nsteps-3; i++) {
        for (int j=0; j<4; j++) {
            coeff_func(i+j, mat[j], vec[j]);
        }
        double z1, z2;
        z1 = y1[i+2] + h_step/24.0 * (
                                   19.0*(mat[2][0]*y1[i+2] + mat[2][1]*y2[i+2] + vec[2][0])
                                   -5.0*(mat[1][0]*y1[i+1] + mat[1][1]*y2[i+1] + vec[1][0])
                                   +    (mat[0][0]*y1[i] + mat[0][1]*y2[i] + vec[0][0])
                                   +9.0*(vec[3][0])
                                   );
        z2 = y2[i+2] + h_step/24.0 * (
                                   19.0*(mat[2][2]*y1[i+2] + mat[2][3]*y2[i+2] + vec[2][1])
                                   -5.0*(mat[1][2]*y1[i+1] + mat[1][3]*y2[i+1] + vec[1][1])
                                   +    (mat[0][2]*y1[i] + mat[0][3]*y2[i] + vec[0][1])
                                   +9.0*(vec[3][1])
                                   );
        for (int j=0; j<4; j++) {
            mat[3][j] *= -h_step * 9.0/24.0;
        }
        mat[3][0] += 1.0;
        mat[3][3] += 1.0;
        solve_2d_linear_equation(mat[3], z1, z2, y1+i+3, y2+i+3);
    }
}

void implicit_adams_4_linear_inward(int nsteps, double * y1, double * y2, ode_2d_linear_coeff coeff_func, double h_step) {
    if (nsteps < 3) {
        std::cerr << "4th-order implicit adams for linear 2d ode requires at least 3 points, but " << nsteps << " were given!" << std::endl;
        exit(EXIT_FAILURE);
    }
    double mat[4][4];
    double vec[4][2];
    for (int i=nsteps-1; i>=3; i--) {
        for (int j=0; j<4; j++) {
            coeff_func(i-j, mat[j], vec[j]);
        }
        double z1, z2;
        z1 = y1[i-2] - h_step/24.0 * (
                                   19.0*(mat[2][0]*y1[i-2] + mat[2][1]*y2[i-2] + vec[2][0])
                                   -5.0*(mat[1][0]*y1[i-1] + mat[1][1]*y2[i-1] + vec[1][0])
                                   +(mat[0][0]*y1[i] + mat[0][1]*y2[i] + vec[0][0])
                                   +9.0*(vec[3][0])
                                   );
        z2 = y2[i-2] - h_step/24.0 * (
                                   19.0*(mat[2][2]*y1[i-2] + mat[2][3]*y2[i-2] + vec[2][1])
                                   -5.0*(mat[1][2]*y1[i-1] + mat[1][3]*y2[i-1] + vec[1][1])
                                   +(mat[0][2]*y1[i] + mat[0][3]*y2[i] + vec[0][1])
                                   +9.0*(vec[3][1])
                                   );
        for (int j=0; j<4; j++) {
            mat[3][j] *= h_step * 9.0/24.0;
        }
        mat[3][0] += 1.0;
        mat[3][3] += 1.0;
        solve_2d_linear_equation(mat[3], z1, z2, y1+i-3, y2+i-3);
    }
}
