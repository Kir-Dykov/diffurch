#pragma once

#include <vector>

using namespace std;

// template
// vector<double> derivative_coefs(const vector<double>& coefs)


template <size_t order, size_t degree>
array<double, degree - order + 1> derivative_coefficients() {
    // cout << "derivative_coefficients is called with order=" << order << endl;
    array<double, degree - order + 1> C;
    for (int i = 0; i <= degree - order; i++) {
        C[i] = 1.;
        for (double j = i+1; j <= i + order; j++) {
            C[i] *= j;
        }
    }
    // cout << C << endl;
    return C;   
}

// TODO: template for numerical type?
template <size_t degree>
class Polynomial {
private:
    array<double, degree+1> coefs;
public:
    
    template <int order=0>
    double eval(double theta) const { 
        if constexpr (order > degree) {
            return 0;
        } 
        
        const static array<double, degree - order + 1> d_coefs = derivative_coefficients<order, degree>();
        
        
        double res = coefs[degree]*d_coefs[degree-order];
        for (int j = degree - 1; j >= order; j--) {
            res *= theta;
            // cout << "j  = " << j << endl << flush;
            res += coefs[j]*d_coefs[j-order];
        }
        // double res = coefs[degree];
        // for (int j = degree - 1; j >= 0; j--) {
        //     res *= theta;
        //     res += coefs[j];
        // }
        return res;
    }
    
    Polynomial(const array<double, degree+1>& coefs) : coefs(coefs) {
        // cout << "Eval(0) = " << eval(0.) << " and Eval(1) = " << eval(1.) << endl;
    };
    
};



enum class RK_method {
    Euler,
    RK4,
    Zero,
    RK5_4_7FM
};


template <RK_method rk = RK_method::Euler>
struct RK {
    /*
    #   Euler   #
    */
    static const int s = 1;
    const vector<vector<double>> A = 
        {{}};
    const vector<double> C =
        {0};
    const vector<double> B =
        {1};
    const vector<vector<double>> BB =
        {{0., 1}}; // linear interpolation coefficients from low to high power
};

template <>
struct RK<RK_method::RK4> {
    /*
    #   RK4   #
    */
    static const int s = 4;
    const vector<vector<double>> A = 
        {{}, 
         {1./2.}, 
         {0.,     1./2.}, 
         {0.,     0.,      1.}};
    const vector<double> C =
        {0,
         1./2.,
         1./2.,
         1.};
    const vector<double> B =
        {1./6., 1./3., 1./3., 1./6.};
    const vector<Polynomial<1>> BB =
        {array{0., 1./6.}, array{0., 1./3.}, array{0., 1./3.}, array{0., 1./6.}}; // linear interpolation coefficients from low to high power
};

template <>
struct RK<RK_method::RK5_4_7FM> {
    /*
    #   RK4   #
    */
    static const int s = 9;
    const vector<vector<double>> A = 
    {
{},
{1./5.},
{3./40.          , 9./40.},
{44./45.         , -56./15.     , 32./9.},
{19372./6561.    , -25360./2187., 64448./6561.     , -212./729.},
{9017./3168.     , -355./33.    , 46732./5247.     , 49./176.          , -5103./18656.},
{35./384.        , 0.           , 500./1113.       , 125./192.         , -2187./6784.   , 11./84.        },
{6245./62208.    , 0.           , 8875./103032.    , -125./1728.       , 801./13568.    , -13519./368064.    , 11105./368064.  },
{632855./4478976., 0.           , 4146875./6391016., 5390625./14183424., -15975./108544., 8295925./220286304., -1779595./62938944., -805./4104.}
    };
    const vector<double> C =
        {0,
        1./5.,
        3./10.,
        4./5.,
        8./9.,
        1.,
        1.,
        1./6.,
        5./6};
    const vector<double> B =
        {5179./57600., 0., 7571./16695., 393./640., -92097./339200., 187./2100., 1./40., 0., 0.};
    const vector<Polynomial<5>> BB =
        {array{0., 1., -38039./7040., 125923./10560., -19683./1760., 3303./880.},
         array{0., 0., 0., 0., 0., 0., }, 
         array{0., 0., -12500./4081., 205000./12243., -90000./4081., 36000./4081.}, 
         array{0., 0., -3125./704., 25625./1056., -5625./176., 1125./88.},
         array{0., 0., 164025./74624., -448335./37312., 295245./18656., -59049./9328.},
         array{0., 0., -25./28., 205./42., -45./7., 18./7.},
         array{0., 0., -2./11., 73./55., -171./55., 108./55.},
         array{0., 0., 189./22., -1593./55., 3537./110., -648./55.},
         array{0., 0., 351./110., -999./55., 2943./110., -648./55.}};
};

template <>
struct RK<RK_method::Zero> {
    /*
    #   RK4   #
    */
    static const int s = 1;
    const vector<vector<double>> A = 
        {{}};
    const vector<double> C =
        {0};
    const vector<double> B =
        {0};
    const vector<Polynomial<0>> BB =
        {array{0.}};
};