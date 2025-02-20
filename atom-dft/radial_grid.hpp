//
//  radial_grid.hpp
//  atom-dft
//
//  Created by ycx on 2024/11/4.
//

#ifndef radial_grid_hpp
#define radial_grid_hpp

#include <iostream>

class radial_grid {
    double r0;
public:
    int N;
    double u_step;
    double * r;
    double * u;
    double * drdu;
    
    void set_exp_grid(double rc, double r0);
    double get_r_by_u_exp_grid(double u);
    
    radial_grid() {
        this->N = 0;
        this->r = nullptr;
        this->u = nullptr;
        this->drdu = nullptr;
    }
    radial_grid(int N) {
        if (N <= 0) {
            std::cerr << "radial_grid.N must be positive! but I found N = " << N << std::endl;
            exit(EXIT_FAILURE);
        }
        this->N = N;
        this->r = new double[N];
        this->u = new double[N];
        this->drdu = new double[N];
    }
    radial_grid(const radial_grid & other) {
        this->N = other.N;
        this->u_step = other.u_step;
        if (other.N != 0) {
            this->r = new double[N];
            this->u = new double[N];
            this->drdu = new double[N];
            memcpy(this->r, other.r, sizeof(double)*N);
            memcpy(this->u, other.u, sizeof(double)*N);
            memcpy(this->drdu, other.drdu, sizeof(double)*N);
        } else {
            this->r = nullptr;
            this->u = nullptr;
            this->drdu = nullptr;
        }
    }
    radial_grid & operator=(const radial_grid & other) {
        if (this == &other) {
            return *this;
        }
        this->u_step = other.u_step;
        if (this->N != 0) {
            delete [] this->r;
            delete [] this->u;
            delete [] this->drdu;
        }
        this->N = other.N;
        if (this->N != 0) {
            this->r = new double[this->N];
            this->u = new double[this->N];
            this->drdu = new double[this->N];
        }
        if (other.N != 0) {
            memcpy(this->r, other.r, sizeof(double) * this->N);
            memcpy(this->u, other.u, sizeof(double) * this->N);
            memcpy(this->drdu, other.drdu, sizeof(double) * this->N);
        }
        return *this;
    }
    radial_grid(radial_grid && other) noexcept {
        this->u_step = other.u_step;
        this->N = other.N;
        this->r = other.r;
        this->u = other.u;
        this->drdu = other.drdu;
        other.N = 0;
        other.r = nullptr;
        other.u = nullptr;
        other.drdu = nullptr;
    }
    radial_grid & operator=(radial_grid && other) noexcept {
        if (this == &other) {
            return *this;
        }
        this->u_step = other.u_step;
        this->N = other.N;
        this->r = other.r;
        this->u = other.u;
        this->drdu = other.drdu;
        other.N = 0;
        other.r = nullptr;
        other.u = nullptr;
        other.drdu = nullptr;
        return *this;
    }
    ~radial_grid() {
        if (this->N != 0) {
            delete [] this->r;
            delete [] this->u;
            delete [] this->drdu;
        }
    }
};

#endif /* radial_grid_hpp */
