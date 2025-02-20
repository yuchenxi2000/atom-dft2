//
//  integrate.hpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#ifndef integrate_hpp
#define integrate_hpp

double simpson(int size, double * y, double * x);
double simpson_uniform(int size, double * y, double step);
void simpson_cumsum(int size, double * y, double * x, double * s);
double trapz_7_uniform(int size, double * y, double step);

#endif /* integrate_hpp */
