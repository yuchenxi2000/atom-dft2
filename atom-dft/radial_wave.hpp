//
//  radial_wave.hpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#ifndef radial_wave_hpp
#define radial_wave_hpp

#ifdef __cplusplus
extern "C" {
#endif

int solve_radial_from_zero(double * chi, double * dchi, int * num_nodes, double E, int l, double Z, double * Vext, double r0, double rc, int N);
int solve_radial_from_zero_and_inf(double * chi, double * dchi, double * dE, int * num_nodes, double E, int l, double Z, double * Vext, double r0, double rc, int N);

#ifdef __cplusplus
}
#endif

#endif /* radial_wave_hpp */
