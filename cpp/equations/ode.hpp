#pragma once

#pragma once

#include "../library/differential_equation.hpp"
#include "../library/vec.hpp"
#include "../library/macros.hpp"
#include "../library/param_names.hpp"

template <size_t n>
using ODE = DifferentialEquation<n, ArgSpec<true,  0, 0, 0>,  ArgSpec<false, 0, 0, 0>, 0>;


struct ODE_lin_1 : public ODE<1> {
    // x' = k x
    double k;
    using param_names = ParamNames<"k">;
    AUTO(param_tuple, tie(k));
    

    ODE_lin_1(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);    
        
        f = { VecMapC<1,1>([this](double x){return k*x;})  };
        analytic_solutions.push_back(VecMapC<1,1>([this](double t){return exp(k*t);}));
    };
        
};

struct ODE_harmonic_oscillator : public ODE<2> {
    // x'' + k^2 x = 0
    double k;
    using param_names = ParamNames<"k">;
    AUTO(param_tuple, tie(k));
    
    ODE_harmonic_oscillator(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);    
        
        f[0] = VecMapC<2,2>([this](const Vec<2>& X){
            const auto& [x, dot_x] = X; 
            return Vec<2>{dot_x, -k*k*x};
        });
        analytic_solutions.push_back(VecMap<1,2>([this](double t){return Vec<2>{sin(k*t), k*cos(k*t)};}));
        
    }
};
    
struct ODE_ln_2 : public ODE<2> {
    // x'' = k^2 x(1 - ln x - ln^2 x)         x = e^sin(k t), x(0) = 1, x'(0) = k   
    double k;
    using param_names = ParamNames<"k">;
    AUTO(param_tuple, tie(k));
    
    ODE_ln_2(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);    
        
        f = { VecMapC<2,2>([this](const Vec<2>& x){double logx = log(x[0]); return Vec<2>{x[1], k*k*x[0]*(1 - logx - logx*logx)};}) };
        
        analytic_solutions.push_back(VecMapC<1,2>([this](double t){return Vec<2>{exp(sin(k*t)), k*cos(k*t)*exp(sin(k*t))};}));
        
    }
    
};


    
struct ODE_Lorenz : public ODE<3> {
    // x'' = k^2 x(1 - ln x - ln^2 x)         x = e^sin(k t), x(0) = 1, x'(0) = k   
    double sigma, beta, rho;
    using param_names = ParamNames<"sigma", "beta", "rho">;
    AUTO(param_tuple, tie(sigma, beta, rho));
    
    ODE_Lorenz(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);    
        f[0] = VecMapC<3,3>([this](const Vec<3>& X){ 
            const auto& [x, y, z] = X;
            return Vec<3>{sigma*(y - x), 
                          rho*x - y - x*z, 
                          -beta*z + x*z}; });
    }
    
};

// Lorenz

// FUNC(t, exp(k*t))
// FUNC((x,y,z), (sigma*(y-x), rho*x - y - x*z, -beta*z+ x*y), (sigma*(dy-dx), rho*dx - dy - x*dz - z*dx, -beta*dz + y*dx + x*dy)) 
// x' = 