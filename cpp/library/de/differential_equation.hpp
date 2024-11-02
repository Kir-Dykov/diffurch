#pragma once

// directly integration related

#include <tuple>
#include <limits>

#include "discoque.hpp"
#include "runge_kutta_tables.hpp"
#include "time_series.hpp"
#include "return_handler.hpp"
#include "return_handler_list.hpp"

#include "../utils/vec.hpp"
#include "../utils/math.hpp"
#include "../utils/find_root.hpp"

using namespace std;

template <size_t n, typename f_spec, typename b_spec, size_t b_levels_n = 0>
class DifferentialEquation {};

template <size_t n,
        bool   f_non_delayed, size_t f_delays_n,  size_t f_neutral_delays_n, size_t f_derivatives_n,
        bool   b_non_delayed,  size_t b_delays_n, size_t b_neutral_delays_n,  size_t b_derivatives_n, size_t b_levels_n>  
class DifferentialEquation<   n,
            ArgSpec<f_non_delayed, f_delays_n, f_neutral_delays_n, f_derivatives_n>,
            ArgSpec<b_non_delayed, b_delays_n, b_neutral_delays_n, b_derivatives_n>,
            b_levels_n > {
public:
                
    static constexpr size_t f_arg_n = n*(f_non_delayed + f_delays_n + f_neutral_delays_n);
    static constexpr size_t b_arg_n = n*(b_non_delayed + b_delays_n + b_neutral_delays_n);
    
    array<VecMapC<f_arg_n, n, f_derivatives_n>, 1+b_levels_n> f;
          VecMapC<b_arg_n, 1, b_derivatives_n>                b;

    array<Real, b_levels_n> b_levels;
    
    ArgSpec<f_non_delayed, f_delays_n, f_neutral_delays_n, f_derivatives_n> f_spec;
    ArgSpec<b_non_delayed, b_delays_n, b_neutral_delays_n, b_derivatives_n> b_spec;     
                
    vector<VecMapC<1,n>> analytic_solutions;
    // vector<VecMapC<1,n>> initial_conditions;
        
    template <typename ReturnHandler, size_t phi_derivatives_n = 0>
    auto solution(Real h, Real t_f, 
                  VecMapC<1, n, phi_derivatives_n> phi, ReturnHandler return_handler = ReturnNone{}) {
    
        Real prev_t, t = 0;
        size_t f_i = 0;
        
        
        DiscoQue<f_delays_n,          6> disco(                f_spec.delays); // Discontinuities Queue
        DiscoQue<f_neutral_delays_n, -1> disco_neutral(f_spec.neutral_delays); // Discontinuities Queue
        
        disco.push(t);         // initial discontinuity, always present
        disco_neutral.push(t); // initial discontinuity, always present
        
        RK_TimeSeries<n, phi_derivatives_n> TS(t, phi);
        
        function<Real(Real)> b_eval = [this, &TS](Real t) {
            return b(TS.eval_pack(t, b_spec));
        };
                
        Real b_val;
        Real b_val_prev;
        
        if constexpr (b_levels_n > 0) {
            b_val = b_eval(t);
            f_i = distance(b_levels.begin(), upper_bound(b_levels.begin(), b_levels.end(), b_val));
        }
        

        while (t < t_f) {
            prev_t = t;
            t = t + h;

            if (t >= t_f) { t = t_f; }
            
            if constexpr (f_delays_n > 0) {
                if (t >= disco.min)         {   t = disco.min; 	                disco.pop(); } // collect propagated discontinuities
            }
            if constexpr (f_neutral_delays_n > 0) {
                if (t >= disco_neutral.min) {   t = disco_neutral.min; 	disco_neutral.pop(); } // collect propagated discontinuities
            }
            
            TS.step_push(t, f[f_i], f_spec);
            
            if constexpr (b_levels_n > 0) { // Discontinuity detection
                b_val_prev = b_val;
                b_val = b_eval(t);          
                if (bool downward =  (f_i != 0          && b_val < b_levels[f_i - 1]); 
                         downward || (f_i != b_levels_n && b_val > b_levels[f_i    ])) 
                {
                    Real t_discontinuity = find_root(b_eval, prev_t, t, b_levels[f_i - downward]);

                            disco.push(t_discontinuity); // add propagated points to mesh in the future
                    disco_neutral.push(t_discontinuity); // add propagated points to mesh in the future

                    TS.pop_back();
                    if (t_discontinuity > std::nextafter(TS.t_ts.back(), std::numeric_limits<Real>::infinity())) {
                        TS.step_push(t_discontinuity, f[f_i], f_spec);
                        if constexpr(has_discontinuity_event_function<ReturnHandler, RK_TimeSeries<n, phi_derivatives_n>>) {
                            return_handler.event_function(TS);
                        }
                    }

                    if (downward) {f_i--;} else {f_i++;}

                    TS.step_push(t, f[f_i], f_spec);

                }
            }
             
            if constexpr(has_event_function<ReturnHandler, RK_TimeSeries<n, phi_derivatives_n>>) {
                return_handler.event_function(TS);
            }
            
        }
        
        if constexpr(has_Return<ReturnHandler, RK_TimeSeries<n, phi_derivatives_n>>) {
            return return_handler.Return(TS);
        }

    }
    
    
};