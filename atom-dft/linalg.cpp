//
//  linalg.cpp
//  atom-dft
//
//  Created by ycx on 2024/11/18.
//

#include "linalg.hpp"

void solve_2d_linear_equation(double A[4], double b1, double b2, double * x1, double * x2) {
    // A    = ( A0  A1 )
    //        ( A2  A3 )
    // Ainv = ( A3 -A1 )
    //        ( -A2 A0 ) / det(A)
    double detA = A[0]*A[3] - A[1]*A[2];
    *x1 = (A[3]*b1 - A[1]*b2) / detA;
    *x2 = (-A[2]*b1 + A[0]*b2) / detA;
}
