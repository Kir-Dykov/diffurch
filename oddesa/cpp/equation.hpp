Мне нравитсяМне нравится#pragma once

#include <functional>

#include "utils/vec.hpp"
#include "utils/tuple.hpp"
#include "utils/RK_tables.hpp"
#include "utils/macros.hpp"
#include "utils/print.hpp"

using namespace std;




template <typename Equation>
struct order_of_equation {
    static constexpr size_t value = std::tuple_size<decltype(std::declval<Equation>().initial_condition(0.))>::value;
};


template <typename Equation>
concept Solvable = requires (Equation eq, double t) { 
     { std::is_same_v<decltype(eq.initial_condition(t)), decltype(eq.lhs(State(t, &eq)))> };
}

template <Solvable Equation>
struct Solver {
    
    using rk = RK98_;
    
    
    // Vec<m> Y; // descrete variables stae
    
    
    // Empty by default
        
    auto solution(double initial_time, double final_time) {
    
        auto self = static_cast<Equation*>(this);
        
        constexpr size_t n = order_of_equation<Equation>::value;
        
        auto& step_events             = self->step_events;
        auto& stop_integration_events = self->stop_integration_events;
        auto& conditional_events      = self->conditional_events;
        auto& stepsize_controller     = self->stepsize_controller;
        
        auto state = State(initial_time);
        
        double h = stepsize_controller.initial_stepsize();
        
        while (t < final_time) {
            
            array<Vec<n>, rk::s> K;
            for (int i = 0; i < rk::s; i++) {
                state.t = prev_t + rk::C[i]*h;
                state.X = prev_X + h*dot(rk::A[i], K, i);
                K[i] -> self->lhs(state);
            }
            
            Vec<n> Ksum = h*dot(rk::B, K, rk::s);
            
            state.t = prev_t + h;
            state.X = prev_X + Ksum1;
            
            Vec<n> Ksum2 = h*dot(rk::B_hat, K, rk::s);
            double error = norm(Ksum-Ksum2);

            h = stepsize_controller.update_stepsize(h, error);
            
            if (stepsize_controller.reject_step(error)) {
                // rejected_step_events
                continue
            }
            
            state.prev_t = state.t;
            state.
            
            
            // /* Runge Kutta Step */ 
            // {
            //     array<Vec<n>, rk::s> K;
            //     for (int i = 0; i < rk::s; i++) {
            //         Vec<n> Ksum = dot(rk::A[i], K, i);
            //         K[i] = self->lhs(prev_t + rk::C[i]*h, X + h*Ksum);
            //     }
            //     Vec<n> Ksum  = dot(rk::B, K, rk::s);
            //     // Vec<n> Ksum_ = dot(rk::B_hat, K, rk::s);
            //     // double error = h*norm(Ksum - Ksum_);
            //     X = X + h*Ksum;
            // }
            // cout << t << endl;
           // cout << X << endl;                        
           for_each(step_events, [](auto& event){ event.save(); });
            
        }
        
        for_each(stop_integration_events, [](auto& event){ event.save(); });

        return tuple_cat(
            transform(step_events,             [](const auto& event){ return std::move(event.get()); }),
            transform(stop_integration_events, [](const auto& event){ return std::move(event.get()); }),
            transform(conditional_events,      [](const auto& event){ return std::move(event.get()); })
        );
    }
};


// vector<Event<n>>

// vector<vector<double>> Event<n>::container

// result = vector<vector<vector<double>>>

struct Lorenz : 
    Solver<Lorenz>, 
    RungeKutta::RK98_,
    Events<
        event::SaveSteps,
        Event<event::SaveWhenStep<[] STATE_FUNC((x,y,z), (return make_tuple(t, x, y, z);))>>,
        Event<
            event::When<[] STATE_FUNC((x,y,z), (return x > 0;)), 
                        [] STATE_FUNC((x,y,z), (return x_(t-1) > 1;))>,
            event::Change<[] STATE_FUNC_NON_CONST((x,y,z), (y = -y; z *= 0.9;))>,
            event::Save<[] STATE_FUNC((x,y,z), (return t;))
        >
    > 
{
                
    
    double sigma;
    double rho;
    double beta;
    
    // need to implement
    // state.but(state.t + c*h, state.X + c*x)
    // that makes shallow copy, with changed fields
    
    auto lhs(const auto& state) {
        const auto& [x, y, z] = state.X;
        return Vec<3>{
            sigma*(y - x),
            rho*x - x*z,
            -beta*z + x*z
        };
    }
    
    auto initial_condition(double t) {
        return Vec<3>{1,2,beta};
    }
    
    Lorenz(const tuple<double, double, double>& params) {
        tie(sigma, rho, beta) = params;
        initial_condition = [this](double t){return Vec<3>{1,2,3};};
    }
    
    StepsizeController stepsize_controller(sigma, rho, beta);
                
                
//                 auto lhs(const State_t& state) {
//         const auto& [x, y] = state.X;
//         return Vec<2>{-y, x};
//     }
    
//     template <size_t derivative = 0>
//     auto initial_condition(double t);
    
//     template <> 
//     auto initial_condition<0>(double t) {
//         return Vec<2>{cos(t), sin(t)};
//     }
//     template <> 
//     auto initial_condition<1>(double t) {
//         return Vec<2>{-sin(t), cos(t)};
//     }
    
    auto events() { 
        using namespace Events;
        return Events(
            Event(
                WhenStep{}, 
                Save(STATE_EVAL__(make_tuple(t, x, y, z))))
            ),
            Event(
                When([this] STATE_EVAL((x,y,z), (x_<1>(t) - a))),
                Save([this] STATE_EVAL((x,y,z), (t)))
            )
        )
    }

};


template <typename Equation>
struct ParametrizedEquation {
    
    // solve for parameter sets
    
    
}