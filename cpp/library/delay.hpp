#pragma once

#include <vector>
#include <functional>
#include "real.hpp"

struct Delay {
    std::function<Real(Real)> tau;
    std::function<std::vector<Real>(Real)> propagate;
    
    Delay() {
        tau       = function<Real(Real)>([](Real){              throw("Delay is called but is uninitialized.");             return 0.; });
        propagate = function<std::vector<Real>(Real)>([](Real){ throw("Delay propagation is called but is uninitialized."); return vector<Real>(); });
    }
    
    Delay(const std::function<Real(Real)>& tau, const std::function<std::vector<Real>(Real)>& propagate_) : tau(tau), propagate(propagate) {};
    
    Delay(Real tau_) : tau([tau_](Real){return tau_;}), propagate([tau_](Real t){return std::vector<Real>{t + tau_};}) {};
    
    inline Real operator()(Real t) const {
        return tau(t);
    }
};


using namespace std;

template<bool non_delayed = true, size_t delays_n = 0, size_t neutral_delays_n = 0, size_t derivatives_n = 0>
struct ArgSpec {
    array<Delay, delays_n> delays;
    array<Delay, neutral_delays_n> neutral_delays;
    
    ArgSpec() {};
    
    ArgSpec(array<Delay, delays_n>                 delays,
            array<Delay, neutral_delays_n> neutral_delays) : delays(delays), neutral_delays(neutral_delays) {};
};


