//
//  poisson.hpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#ifndef poisson_hpp
#define poisson_hpp

#include "radial_grid.hpp"

#ifdef __cplusplus
extern "C" {
#endif

void solve_poisson(double * VH, double * dVH, double * density_r, double r0, double rc, int N);

#ifdef __cplusplus
}
#endif

#endif /* poisson_hpp */
