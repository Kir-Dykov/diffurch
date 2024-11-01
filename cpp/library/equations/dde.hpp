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



// x' = k (1- theta) x + k e^tau theta x(t - tau)
// has solution x = e^{kt}
struct DDE_lin_1 : public DDE<1, 1> {
    PARAMS(DDE_lin_1, theta, k, tau);
    
    DDE_lin_1() {
        f[0] = VEC_LAMBDA((this), (x, x_tau), ((k*(1. - theta))*x + (k * exp(tau) * theta)*x_tau));
        f_spec = {{Delay(tau)}, {}};
        analytic_solutions.push_back( VEC_LAMBDA((this), (t), (exp(k*t))) );
    }
};


struct DDE_harmonic : public DDE<2, 1> {
    PARAMS(DDE_harmonic, theta, k, tau);
    
    DDE_harmonic() {
        f[0] = VEC_LAMBDA((this), (x, dot_x, x_tau, dot_x_tau), (
            dot_x,
            -(1-theta)*(k*k*x) - (theta)*(k * sin(k*tau) *dot_x_tau + k*k* cos(k*tau)*x_tau)
        ));
        f_spec = {{Delay(tau)}, {}};
        analytic_solutions.push_back(VEC_LAMBDA((this), (t), (sin(k*t), k*cos(k*t))));
    }
};

// x' = - x ln(x(t-pi/2))      has solution x = e^{C sin(t)}
struct DDE_ln_1 : public DDE<1, 1> {
    PARAMS(DDE_ln_1, C);
    
    DDE_ln_1() {
        f[0] = VEC_LAMBDA((this), (x, x_tau), (-x * log(x_tau)));
        f_spec.delays = {Delay(M_PI * 0.5)};
        analytic_solutions.push_back(VEC_LAMBDA((this), (t), (exp(C*sin(t)))));
    }    
};

