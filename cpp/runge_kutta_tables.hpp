#pragma once

#include <vector>

using namespace std;

enum RK_method {
    Euler,
    RK4,
    Zero,
    RK5_4_7FM
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
struct RK<RK5_4_7FM> {
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
    const vector<vector<double>> B_polynomial_coefs =
        {{1., -38039./7040., 125923./10560., -19683./1760., 3303./880.},
         {0., 0., 0., 0., 0., }, 
         {0., -12500./4081., 205000./12243., -90000./4081., 36000./4081.}, 
         {0., -3125./704., 25625./1056., -5625./176., 1125./88.},
         {0., 164025./74624., -448335./37312., 295245./18656., -59049./9328.},
         {0., -25./28., 205./42., -45./7., 18./7.},
         {0., -2./11., 73./55., -171./55., 108./55.},
         {0., 189./22., -1593./55., 3537./110., -648./55.},
         {0., 351./110., -999./55., 2943./110., -648./55.}};
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