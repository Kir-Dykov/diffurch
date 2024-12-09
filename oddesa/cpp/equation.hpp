#pragma once

#include <functional>

#include "utils/vec.hpp"
#include "utils/tuple.hpp"
#include "utils/RK_tables.hpp"
#include "utils/macros.hpp"
#include "utils/print.hpp"

using namespace std;


// struct Event {
    
//     bool detect();
//     bool condition();
//     void locate();
//     void save();
//     void change();
//     void change_above();
//     void change_below();
//     void print();
//     void action(); // step_on
// }


template <typename ReturnType>
struct Event {
    
    function<bool()> detect;
    function<bool()> locate;
    function<void()> change;
    function<ReturnType()> save;
};

template <typename... Types>
struct Event<tuple<Types...>> {
    function<bool()> detect;
    function<bool()> locate;
    function<void()> change;
    function<tuple<Types...>()> what_to_save;
    tuple<vector<Types>...>  container;
    void save() {
        auto X = what_to_save();
        push_back(container, X);
    }
};




// FIRSTLY
// ODE
// NO DISCONTINUITIES
// NO variable stepsize
// YES : events for return handlers
template <size_t n, typename Derived>
struct IVP {
    
    using rk = RK98_;
    
    // using var = variable_names<"x", "y", "z">;
    
    function<Vec<n>(double, const Vec<n>&)> lhs;
    function<Vec<n>(double)> initial_condition; // initial value depends on the initial time
    pair<double, double> integration_interval;
    
    double t; // current time
    Vec<n> X; // current state
    Vec<n> prev_X; // current state
    
    // Vec<m> Y; // descrete variables stae
    
    
    // Empty by default
    tuple<> step_events             = make_tuple();
    tuple<> stop_integration_events = make_tuple();
    tuple<> conditional_events      = make_tuple();
        
    auto solution() {
        
        t = integration_interval.first; // current time
        X = initial_condition(t); // current state
        
        auto& step_events             = static_cast<Derived*>(this)->step_events;
        auto& stop_integration_events = static_cast<Derived*>(this)->stop_integration_events;
        auto& conditional_events      = static_cast<Derived*>(this)->conditional_events;
        
        while (t < integration_interval.second) {
            
            double prev_t = t;
            t = prev_t + 0.001; // no stepsize controller yet
            
            double h = t - prev_t;
            
            /* Runge Kutta Step */ 
            {
                array<Vec<n>, rk::s> K;
                for (int i = 0; i < rk::s; i++) {
                    Vec<n> Ksum = dot(rk::A[i], K, i);
                    K[i] = lhs(t + rk::C[i]*h, X + h*Ksum);
                }
                Vec<n> Ksum  = dot(rk::B,     K, rk::s);
                // Vec<n> Ksum_ = dot(rk::B_hat, K, rk::s);
                // double error = h*norm(Ksum - Ksum_);
                X = X + h*Ksum;
            }
            // cout << t << endl;
           // cout << X << endl;                        
           for_each(step_events, [](auto& event){ event.save(); });
            
        }
        
        for_each(stop_integration_events, [](auto& event){ event.save(); });

        return tuple_cat(
            transform(step_events,             [](const auto& event){ return event.container; }),
            transform(stop_integration_events, [](const auto& event){ return event.container; }),
            transform(conditional_events,      [](const auto& event){ return event.container; })
        );
    }
};




struct Lorenz : IVP<3, Lorenz> {
    double sigma;
    double rho;
    double beta;
    
    Lorenz(const tuple<double, double, double>& params) {
        tie(sigma, rho, beta) = params;
        
        lhs = [this](double t, Vec<3> XX){
            const auto& [x, y, z] = XX;
            return Vec<3>{
                sigma*(y - x),
                rho*x - x*z,
                -beta*z + x*z
            };
        };
        
        initial_condition = [this](double t){return Vec<3>{1,2,3};};
        
        integration_interval = make_pair(0., 0.1);
    }

    tuple<Event<tuple<double, double, double, double>>> step_events = make_tuple(
        Event<tuple<double, double, double, double>> {
            .what_to_save = [this]() -> tuple<double, double, double, double> {
                const auto& [x, y, z] = X;
                return make_tuple(t, x, y, z);
            }
        }
    );
};

        
