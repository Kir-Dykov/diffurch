#pragma once

// related to discoque and integration and definition of ddes

#include <vector>
#include <functional>

#include "../utils/real.hpp"

using namespace std;

struct Delay {
    std::function<Real(Real)> tau;
    std::function<std::vector<Real>(Real)> propagate;
    
    Delay() {
        tau       = std::function<Real(Real)>([](Real){              throw("Delay is called but is uninitialized.");             return 0.; });
        propagate = std::function<std::vector<Real>(Real)>([](Real){ throw("Delay propagation is called but is uninitialized."); return std::vector<Real>(); });
    }
    
    template <typename Callable>
    Delay(const Callable& tau, const Callable& propagate) : tau(tau), propagate(propagate) {};
    
    Delay(const Real&& tau_) : tau([ tau_](Real){return tau_;}), propagate([ tau_](Real t){return std::vector<Real>{t + tau_};}) {}; // copy if rvalue
    Delay(const Real&  tau_) : tau([&tau_](Real){return tau_;}), propagate([&tau_](Real t){return std::vector<Real>{t + tau_};}) {}; // reference if lvalue
    
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


