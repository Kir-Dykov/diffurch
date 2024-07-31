#pragma once

#include <vector>

using namespace std;

enum RK_method {
    Euler,
    RK4,
    Zero
};


template <RK_method rk = Euler>
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
    const vector<vector<double>> B_polynomial_coefs =
        {{1}}; // linear interpolation coefficients from low to high power
};

template <>
struct RK<RK4> {
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
    const vector<vector<double>> B_polynomial_coefs =
        {{1./6.}, {1./3.}, {1./3.}, {1./6.}}; // linear interpolation coefficients from low to high power
};


template <>
struct RK<Zero> {
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
    const vector<vector<double>> B_polynomial_coefs =
        {{0}}; // linear interpolation coefficients from low to high power
};