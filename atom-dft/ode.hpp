//
//  ode.hpp
//  atom-dft
//
//  Created by ycx on 2024/11/18.
//

#ifndef ode_hpp
#define ode_hpp

#include <functional>

typedef std::function<void (int, double, double, double *, double *)> ode_2d_func;
void runge_kutta_4_outward(int nsteps, double * chi, double * dchi, ode_2d_func func, ode_2d_func func_mid, double h_step);
void runge_kutta_4_inward(int nsteps, double * chi, double * dchi, ode_2d_func func, ode_2d_func func_mid, double h_step);

void adam_4_outward(int nsteps, double * y1, double * y2, ode_2d_func func, double h_step);
void adam_4_inward(int nsteps, double * y1, double * y2, ode_2d_func func, double h_step);

typedef std::function<void (int, double *, double *)> ode_2d_linear_coeff;
void implicit_adams_4_linear_outward(int nsteps, double * y1, double * y2, ode_2d_linear_coeff coeff_func, double h_step);
void implicit_adams_4_linear_inward(int nsteps, double * y1, double * y2, ode_2d_linear_coeff coeff_func, double h_step);

#endif /* ode_hpp */
