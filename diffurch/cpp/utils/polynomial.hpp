#pragma once


using namespace std;

template <size_t order, size_t degree>
array<double, degree - order + 1> derivative_coefficients() {
    array<double, degree - order + 1> C;
    for (int i = 0; i <= degree - order; i++) {
        C[i] = 1.;
        for (double j = i+1; j <= i + order; j++) {
            C[i] *= j;
        }
    }
    return C;   
}

template <size_t degree>
class Polynomial {
private:
    array<double, degree+1> coefs;
public:   
    
    Polynomial(const array<double, degree+1>& coefs) : coefs(coefs) {};
    Polynomial(std::initializer_list<double> list) : coefs(list) {};
    
    
    
    template <int order=0>
    double eval(double theta) const { 
        if constexpr (order > degree) {
            return 0;
        } 
        const static array<double, degree - order + 1> d_coefs = derivative_coefficients<order, degree>();
        double res = coefs[degree]*d_coefs[degree-order];
        for (int j = degree - 1; j >= order; j--) {
            res *= theta;
            res += coefs[j]*d_coefs[j-order];
        }
        // double res = 0r;
        // double power = 1r;
        // for (int j = 0; j <= degree - order; j++) {
        //     res += power * coefs[j+order]*d_coefs[j];
        //     power *= theta;
        // }
        return res;
    }
    
};
