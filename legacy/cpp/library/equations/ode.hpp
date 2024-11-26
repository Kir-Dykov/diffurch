#pragma once

#pragma once

#include "../include.hpp"


template <size_t n>
using ODE = DifferentialEquation<n, ArgSpec<true,  0, 0, 0>,  ArgSpec<false, 0, 0, 0>, 0>;


struct ODE_lin_1 : public ODE<1> {
    PARAMS(ODE_lin_1, k)

    ODE_lin_1() {
        f[0] = VEC_LAMBDA((this), (x), (k*x));
        analytic_solutions.push_back(VEC_LAMBDA((this), (t), (exp(k*t))));
    };
        
};

// x'' + k^2 x = 0
struct ODE_harmonic_oscillator : public ODE<2> {
    PARAMS(ODE_harmonic_oscillator, k)
    
    ODE_harmonic_oscillator() {
        f[0] = [this](VecArg<2> X){
            const auto& [x, dot_x] = X; 
            return Vec<2>{dot_x, -k*k*x};
        };
        analytic_solutions.push_back([this](double t){return Vec<2>{sin(k*t), k*cos(k*t)};});        
    }
    
};

// x'' = k^2 x(1 - ln x - ln^2 x)         x = e^sin(k t), x(0) = 1, x'(0) = k  
struct ODE_ln_2 : public ODE<2> {
    PARAMS(ODE_ln_2, k)
    
    ODE_ln_2() {
        f[0] = [this](const Vec<2>& x){Real logx = log(x[0]); return Vec<2>{x[1], k*k*x[0]*(1 - logx - logx*logx)};};
        analytic_solutions.push_back([this](Real t){return Vec<2>{exp(sin(k*t)), k*cos(k*t)*exp(sin(k*t))};});
    }
    
};

    
// x'' = k^2 x(1 - ln x - ln^2 x)         x = e^sin(k t), x(0) = 1, x'(0) = k   
struct ODE_Lorenz : public ODE<3> {
    PARAMS(ODE_Lorenz, sigma, beta, rho)
    
    ODE_Lorenz() {
        f[0] = VEC_LAMBDA((this), 
                          (x,y,z), 
                          (
                              sigma*(y - x), 
                              rho*x - y - x*z, 
                              -beta*z + x*z
                          ));
    }
};

