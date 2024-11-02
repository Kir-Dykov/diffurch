#pragma once

#include "../include.hpp"
#include "dde.hpp"

// template <size_t n, size_t m>
// using DDE = DifferentialEquation<n, ArgSpec<true,  m, 0, 0>,  ArgSpec<false, 0, 0, 0>, 0>;

// Continuous DDE
// x' = x/2 + e/2 x(t - 1)     has solution x = C e^t
// x' = e x(t - 1)             has solution x = C e^t
// x' = - x ln(x(t-pi/2))      has solution x = e^{C sin(t)}

struct DDE_relay_1 : public DifferentialEquation<1, ArgSpec<true,  0, 0, 0>,  ArgSpec<false, 1, 0, 0>, 1> {
    // x' = k x + alpha sign(x_tau), alpha < 0
    double alpha, k, tau;
    using param_names = ParamNames<"alpha", "k", "tau">;
    AUTO(param_tuple, tie(alpha, k, tau));
    
    DDE_relay_1(const tuple_add_const<decltype(param_tuple)>& params) {
        copy_tuple(params, param_tuple);   

        f[0] = VecMapC<1,1>([this](const Vec<1>& x){
            return k*x - alpha;
        });
        f[1] = VecMapC<1,1>([this](const Vec<1>& x){
            return k*x + alpha;
        });
        
        b = VecMapC<1,1>([](Vec<1> x_tau){return x_tau;});
        b_levels = {0.};
        b_spec = {{Delay(tau)}, {}};
    
        
        // analytic solution
        /*
            if initial function phi(t) is negative except for phi(0) = 0, then  
        */
        analytic_solutions.push_back(VecMap<1,1>([this](Real t){
            Real a, b;
            
            if (k == 0) {a = -tau; b = 3*tau; } 
            else { a = tau + log(-1+2*exp(-k*tau))/k; b = tau - log(-1+2*exp(-k*tau))/k; }
            
            t = fmod(t - a, b - a); // Map x to interval [a, b)
            t += (b-a)*(t < 0) + a;

            if (k == 0) {
                if (t <= tau) {
                    return -t*alpha;
                } else {
                    return (t - 2 * tau)*alpha;
                }
            } else {
                if (t <= tau) {
                    return (1- exp(k*t))*alpha/k;
                } else {
                    return -(1 + (1 - 2*exp(-k*tau))*exp(k*t))*alpha/k;
                }
            }
        }));
    }
};






struct ResurgentNeuron : public DDE<2, 2>
{
    PARAMS(ResurgentNeuron, alpha, beta, eta, xi, lambda, h);
    
    Real f_alpha(Real u) {
        return alpha*(1.r - u)/(alpha+u);
    }
    
    Real f_beta(Real u) {
        return beta*(1.r - u)/(beta+u);
    }
    
    Real g(Real u) {
        return xi*(u-eta)/(u+xi);
    }
    
    ResurgentNeuron() {
        f[0] = VEC_LAMBDA((this), (x, y, x_1, y_1, x_h, y_h), (
            f_alpha(exp(lambda*x_1)),
            f_beta( exp(lambda*y_h)) + g(exp(lambda*x))
        ));
        f_spec = {{Delay(1.r), Delay(h)}, {}};
    }
};
