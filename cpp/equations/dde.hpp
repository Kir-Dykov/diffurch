#pragma once

#pragma once

#include "../library/differential_equation.hpp"
#include "../library/vec.hpp"

template <size_t n, size_t m>
using DDE = DifferentialEquation<n, ArgSpec<true,  m, 0, 0>,  ArgSpec<false, 0, 0, 0>, 0>;

// Continuous DDE
// x' = x/2 + e/2 x(t - 1)     has solution x = C e^t
// x' = e x(t - 1)             has solution x = C e^t
// x' = - x ln(x(t-pi/2))      has solution x = e^{C sin(t)}

struct DDE_lin_1 : public DDE<1, 1> {
    // x' = k (1- theta) x + k e^tau theta x(t - tau)
    // has solution x = e^{kt}
    double theta, k, tau;
    using param_names = ParamNames<"theta", "k", "tau">;
    AUTO(param_tuple, tie(theta, k, tau));
    
    DDE_lin_1(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   
        
        double A = k*(1. - theta);
        double B = k * exp(tau) * theta;
        f[0] = VecMapC<2,1>([A, B](const Vec<2>& X){
            const double& x     = X[0]; 
            const double& x_tau = X[1];
            return A*x +  B*x_tau;
        });
        f_spec = {{Delay(tau)}, {}};
    
        analytic_solutions.push_back(VecMap<1,1>([this](double t){return exp(k*t);}));
    }
    
    
    
};

// struct DDE_ln_1 : public DDE<1, 1> {
//     // x' = - x ln(x(t-pi/2))      has solution x = e^{C sin(t)}
//     DDE_ln_1() {
//         f = { VecMapC<2,1>([this](const Vec<2>& x){return - x[0] * log(x[1]);})  };
//         f_spec.delays = {Delay(M_PI * 0.5)};
//     }
    
//     VecMapC<1,1> analytic_solution {[this](double t){return exp(sin(t));}};
// };

