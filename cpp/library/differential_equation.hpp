#pragma once

#include <tuple>
#include <limits>

#include "discoque.hpp"
#include "types.hpp"
#include "vec.hpp"
#include "runge_kutta_tables.hpp"
#include "time_series.hpp"
#include "utils.hpp"
#include "find_root.hpp"



using namespace std;

/**
 * @brief Enumeration of possible return options for DifferentialEquation::solution function.
 */
enum class ReturnOption {
    None, All, AllInit, Last
};


     

struct ReturnSolution {
    template <size_t n, size_t phi_derivatives_n>
    static auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        return make_tuple(TS.t_ts, TS.x_ts);
    };
};

struct ReturnLastValue {
    template <size_t n, size_t phi_derivatives_n>
    static auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        return TS.x_ts.back(); 
    };
};

struct ReturnNone{
    template <size_t n, size_t phi_derivatives_n>
    static void Return(RK_TimeSeries<n, phi_derivatives_n>& TS) {};
};

/**
 * @brief The base class for DifferentialEquation derived objects.
 *
 * This class holds the right-hand-side function as well as other information about the differential equation.
 *
 * @tparam n The order of the system, namely the dimension of the solution function x(t).
 * @tparam f_spec ArgSpec templated type.
 */
template <size_t n, typename f_spec, typename b_spec, size_t b_levels_n = 0>
class DifferentialEquation {};


template <size_t n,          // order of the system
        bool   f_non_delayed, size_t f_delays_n,  size_t f_neutral_delays_n, size_t f_derivatives_n,
        bool   b_non_delayed,  size_t b_delays_n, size_t b_neutral_delays_n,  size_t b_derivatives_n, size_t b_levels_n>   // number of discontinuity surfaces
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

    // template <ReturnOption return_option = ReturnOption::None, size_t phi_derivatives_n = 0>
    // auto solution(Real h, Real t_f, 
    //               VecMapC<1, n, phi_derivatives_n> phi,
    //               function<void(RK_TimeSeries<n, phi_derivatives_n>&)> event_function = [](RK_TimeSeries<n, phi_derivatives_n>& TS){return;}) {    
        
    template <typename ReturnType, size_t phi_derivatives_n = 0>
    auto solution(Real h, Real t_f, 
                  VecMapC<1, n, phi_derivatives_n> phi,
                  function<void(RK_TimeSeries<n, phi_derivatives_n>&)> event_function = [](RK_TimeSeries<n, phi_derivatives_n>& TS){return;}) {    
        
    
        Real prev_t, t = 0;
        size_t f_i = 0;
        
        DiscoQue<f_delays_n,          6> disco(                f_spec.delays); // Discontinuities Queue
        DiscoQue<f_neutral_delays_n, -1> disco_neutral(f_spec.neutral_delays); // Discontinuities Queue
        disco.push(t);         // initial discontinuity, always present
        disco_neutral.push(t); // initial discontinuity, always present
        
        RK_TimeSeries<n, phi_derivatives_n> TS(t, phi);
        
        auto b_eval = [this, &TS](Real t) {
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
                if (bool upward =  (f_i != b_levels_n && b_val > b_levels[f_i    ]); 
                         upward || (f_i != 0          && b_val < b_levels[f_i - 1])) 
                {
                    Real t_discontinuity = find_root(b_eval, prev_t, t, b_levels[f_i - !upward]);

                            disco.push(t_discontinuity); // add propagated points to mesh in the future
                    disco_neutral.push(t_discontinuity); // add propagated points to mesh in the future

                    TS.pop_back();
                    if (t_discontinuity > std::nextafter(TS.t_ts.back(), std::numeric_limits<Real>::infinity())) {
                        TS.step_push(t_discontinuity, f[f_i], f_spec);
                        event_function(TS);
                    }

                    if (upward) {f_i++;} else {f_i--;}

                    TS.step_push(t, f[f_i], f_spec);

                }
            }
             
            event_function(TS);
        }
            
        return ReturnType::Return(TS);

//         if      constexpr (return_option == ReturnOption::Last) {
//             return TS.x_ts.back(); 
//         } else if constexpr (return_option == ReturnOption::All ) {
//             return make_tuple(TS.t_ts, TS.x_ts);
//         } else if constexpr (return_option == ReturnOption::AllInit ) {
//             size_t size = TS.t_ts.size();
//             vector<vector<Real>> result;
//             Real t = 0;
//             for (auto tau : f_spec.delays) {
//                 t = min(t, -tau);
//             }
//             for (auto tau : b_spec.delays) {
//                 t = min(t, -tau);
//             }
            
//             while(t <= t_f) {
//                 vector<Real> v(n+1);
//                 v[0] = t;
//                 auto x = TS.interpolation(t);
//                 VecCopy(x, v, 1);
//                 result.push_back(v);
//                 t += 0.0001;
//             }
//             return result;
//         } else {
//             return;
//         }
    }
    
    
};