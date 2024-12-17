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


template <size_t n, typename ReturnType>
struct Event {};

template <size_t n, typename... Types>
struct Event<n, tuple<Types...>> {
    function<bool(double, const Vec<n>&)> detect;
    function<bool(double, const Vec<n>&)> filter;
    function<void()> change;
    function<tuple<Types...>(double, const Vec<n>&)> what_to_save;
    
    using container_t = tuple<vector<Types>...>;
    container_t container;
    
    void save() {
        auto X = what_to_save();
        push_back(container, X);
    }
    
    Event(
        optional<function<bool(double, const Vec<n>&)>> detect_ = nullopt,
        optional<function<bool(double, const Vec<n>&)>> filter_ = nullopt,
        optional<function<void()>> change_ = nullopt,
        optional<function<tuple<Types...>(double, const Vec<n>&)>> what_to_save_ = nullopt    
    ) {
        if (detect_.has_value()) detect = detect_;
        if (filter_.has_value()) filter = filter_;
        if (change_.has_value()) change = change_;
        if (what_to_save_.has_value()) what_to_save = what_to_save_;
    }
};

template <size_t n>
struct Event<n, void> {
    function<bool(double, const Vec<n>&)> detect;
    function<bool(double, const Vec<n>&)> filter;
    function<void()> change;
    function<void(double, const Vec<n>&)> what_to_save;
    
    using container_t = void;
    
    void save() {}
    
    Event(
        optional<function<bool(double, const Vec<n>&)>> detect_ = nullopt,
        optional<function<bool(double, const Vec<n>&)>> filter_ = nullopt,
        optional<function<void()>> change_ = nullopt,
        optional<function<void(double, const Vec<n>&)>> what_to_save_ = nullopt    
    ) {
        if (detect_.has_value()) detect = detect_;
        if (filter_.has_value()) filter = filter_;
        if (change_.has_value()) change = change_;
        if (what_to_save_.has_value()) what_to_save = what_to_save_;
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
    
    double t; // current time
    double prev_t; // previous time
    Vec<n> X; // current state
    Vec<n> prev_X; // prevouous state
    
    // Vec<m> Y; // descrete variables stae
    
    
    // Empty by default
    tuple<> step_events             = make_tuple();
    tuple<> stop_integration_events = make_tuple();
    tuple<> conditional_events      = make_tuple();
        
    auto solution(double initial_time, double final_time) {
        
        t = initial_time; // current time
        X = initial_condition(t); // current state
        
        auto& step_events             = static_cast<Derived*>(this)->step_events;
        auto& stop_integration_events = static_cast<Derived*>(this)->stop_integration_events;
        auto& conditional_events      = static_cast<Derived*>(this)->conditional_events;
        
        while (t < final_time) {
            
            prev_t = t;
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


// vector<Event<n>>

// vector<vector<double>> Event<n>::container

// result = vector<vector<vector<double>>>

struct Lorenz : IVP<3, Lorenz>, Event<[](){}, [](){}> {
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
        
    }
    

    // tuple<Event<tuple<double, double, double, double>>> step_events = make_tuple(
    //     Event<tuple<double, double, double, double>> {
    //         .what_to_save = [this](){
    //             const auto& [x, y, z] = X;
    //             return make_tuple(t, x, y, z);
    //         }
    //     }
    // );
    
    Event(
        [this](double t, Vec<n> X){
                const auto& [x, y, z] = X;
                return make_tuple(t, x, y, z);
        }
    ).set_detect(
        [this](double t, Vec<n> X){
            const auto& [x, y, z] = X; 
            return x > 0;
        }
    )
        
    // Event<n>().set_detect(...).set_change(...)
        
        
            
    // NON_STATIC_AUTO_DECLARATION(step_events, make_tuple(
    //     Event(nullopt, 
    //           nullopt, 
    //           nullopt, 
    //           [this](){
    //             const auto& [x, y, z] = X;
    //             return make_tuple(t, x, y, z);
    //           })
    //     )
    // );
};

        
// tuple<Event<tuple<{", ".join(["double"]*(variable_n+1))}>>> step_events = make_tuple(
//         Event<tuple<{", ".join(["double"]*(variable_n+1))}>> {{
//             .what_to_save = [this]() {{
//                 const auto& [{variable_list}] = X;
//                 return make_tuple(t, {variable_list});
//             }}
//         }}
//     );