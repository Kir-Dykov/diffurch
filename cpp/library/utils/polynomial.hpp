#pragma once

// tool

#include "real.hpp"

using namespace std;

template <size_t order, size_t degree>
array<Real, degree - order + 1> derivative_coefficients() {
    array<Real, degree - order + 1> C;
    for (int i = 0; i <= degree - order; i++) {
        C[i] = 1.;
        for (Real j = i+1; j <= i + order; j++) {
            C[i] *= j;
        }
    }
    return C;   
}

template <size_t degree>
class Polynomial {
private:
    array<Real, degree+1> coefs;
public:   
    
    Polynomial(const array<Real, degree+1>& coefs) : coefs(coefs) {};
    Polynomial(std::initializer_list<Real> list) : coefs(list) {};
    
    
    
    template <int order=0>
    Real eval(Real theta) const { 
        if constexpr (order > degree) {
            return 0;
        } 
        const static array<Real, degree - order + 1> d_coefs = derivative_coefficients<order, degree>();
        Real res = coefs[degree]*d_coefs[degree-order];
        for (int j = degree - 1; j >= order; j--) {
            res *= theta;
            res += coefs[j]*d_coefs[j-order];
        }
        // Real res = 0r;
        // Real power = 1r;
        // for (int j = 0; j <= degree - order; j++) {
        //     res += power * coefs[j+order]*d_coefs[j];
        //     power *= theta;
        // }
        return res;
    }
    
};
