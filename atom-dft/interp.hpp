//
//  interp.hpp
//  atom-dft
//
//  Created by ycx on 2024/11/12.
//

#ifndef interp_hpp
#define interp_hpp

double interp_poly(double z, int n, double * y, double * x);

void get_mid_array_poly(double * mid, int size, double * y, double * x, int ninterp);

#endif /* interp_hpp */
