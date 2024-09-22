#pragma once

#include <tuple>
#include <limits>

#include "discoque.hpp"
#include "types.hpp"
#include "vec.hpp"
#include "runge_kutta_tables.hpp"
#include "runge_kutta_time_series.hpp"
#include "utils.hpp"



using namespace std;

enum class ReturnOption {
    None,
    All,
    AllInit,
    Last,
    DiscontinuityCount
};



template <size_t n, typename f_spec, typename b_spec, size_t b_levels_n = 0>
class DE {};


template <size_t n,          // order of the system
        bool   f_non_delayed,
        size_t f_delays_n,   // number of delayed arguments in f
        size_t f_neutral_delays_n,   // number of delayed arguments in f
        size_t f_derivatives_n,
        bool   b_non_delayed, // is there a non-delayed argument in b
        size_t b_delays_n,   // number of delayed arguments in b
        size_t b_neutral_delays_n,   // number of delayed arguments in b
        size_t b_levels_n,
        size_t b_derivatives_n>   // number of discontinuity surfaces
class DE<   n,
            ArgSpec<f_non_delayed, f_delays_n, f_neutral_delays_n, f_derivatives_n>,
            ArgSpec<b_non_delayed, b_delays_n, b_neutral_delays_n, b_derivatives_n>,
            b_levels_n
        > {
public:
                
    // DE () = delete; // no constructor, is to be used as a base class only
                
    static constexpr size_t f_arg_n = n*(f_non_delayed + f_delays_n + f_neutral_delays_n);
    static constexpr size_t b_arg_n = n*(b_non_delayed + b_delays_n + b_neutral_delays_n);
                
    
    array<VecMapC<f_arg_n, n, f_derivatives_n>, 1+b_levels_n> f;
          VecMapC<b_arg_n, 1, b_derivatives_n>                b;

    array<double, b_levels_n> b_levels;
    
    ArgSpec<f_non_delayed, f_delays_n, f_neutral_delays_n, f_derivatives_n> f_spec;
    ArgSpec<b_non_delayed, b_delays_n, b_neutral_delays_n, b_derivatives_n> b_spec;

    
    template <ReturnOption return_option = ReturnOption::None, size_t phi_derivatives_n = 0>
    auto solution(double h, double t_f, VecFuncC<n, phi_derivatives_n> phi,
             function<void(RK_TimeSeries<n, phi_derivatives_n>&)> event_function = [](RK_TimeSeries<n, phi_derivatives_n>& TS){return;}) {    
        
            
        double prev_t, t;
        t = 0;
        
        DiscoQue<f_delays_n,         6> disco(                f_delays); // Discontinuities Queue
        DiscoQue<f_neutral_delays_n, 6> disco_neutral(f_neutral_delays); // Discontinuities Queue
        
        RK_TimeSeries<n, phi_derivatives_n> TS(t, phi);
        
        auto b_eval = [this, &TS](double t) {
            return b(TS.eval_pack(t, b_spec));
        };
                
        auto find_root = [this, &TS, &b_eval](double prev_t, double t, double b_level, bool is_upward){
            double t_discontinuity;
            double l = prev_t; double r = t;
            // cout << array{l, r} << " --> " << array{b_eval(l), b_eval(r)} << endl;
            for (int i = 0; i < 50; i++) {
                t_discontinuity = (r+l)*0.5;
                if ((b_eval(t_discontinuity) > b_level) == !is_upward) 
                { l = t_discontinuity; }  else  { r = t_discontinuity; }
            }
            // cout << "\t" << t_discontinuity << " --> " << b_eval(t_discontinuity) << endl;
            return t_discontinuity;
        };
            
        
        double b_val = b_eval(t);
        double b_val_prev = b_val;
        
        // the least index, such that b_val < b_levels[f_i]
        size_t f_i = distance(b_levels.begin(), upper_bound(b_levels.begin(), b_levels.end(), b_val));

        while (t < t_f) {
            prev_t = t;
            b_val_prev = b_val;
            
            t = t + h;

            if (t >= t_f) { t = t_f; }
            if (t >= disco.min)         {   t = disco.min; 	                disco.pop(); } // collect propagated discontinuities
            if (t >= disco_neutral.min) {   t = disco_neutral.min; 	disco_neutral.pop(); } // collect propagated discontinuities
            
            TS.step_push(t, f[f_i], f_spec);
            b_val = b_eval(t);
            
            if (bool upward =  (f_i != b_levels_n && b_val > b_levels[f_i    ]); 
                     upward || (f_i != 0          && b_val < b_levels[f_i - 1])) 
            {
                double t_discontinuity = find_root(prev_t, t, b_levels[f_i - !upward], upward);
                
                disco.push(t_discontinuity); // add propagated points to mesh in the future
                disco_neutral.push(t_discontinuity); // add propagated points to mesh in the future
                
                TS.pop_back();
                if (t_discontinuity > std::nextafter(TS.t_ts.back(), std::numeric_limits<double>::infinity())) {
                    TS.step_push(t_discontinuity, f[f_i], f_spec);
                    event_function(TS);
                }
                
                if (upward) {f_i++;} else {f_i--;}
                TS.step_push(t,               f[f_i], f_spec);
                
                // if constexpr (return_option == ReturnOption::DiscontinuityCount) {
                //     if (TS.t_ts.back() > 0.5*t_f)
                //         discontinuity_count++;
                // }
            }
            
            event_function(TS);
        }
            


        if      constexpr (return_option == ReturnOption::Last) {  return TS.x_ts.back(); } 
        else if constexpr (return_option == ReturnOption::All ) {
            size_t size = TS.t_ts.size();
            vector<vector<double>> result(size, vector<double>(1+n));
            for (int i = 0; i < size; i++) {
                result[i][0] = TS.t_ts[i];
                VecCopy(TS.x_ts[i], result[i], 1);
                // copy(TS.x_ts[i].begin(), TS.x_ts[i].end(), result[i].begin()+1);
            }
            return result;
        } else if constexpr (return_option == ReturnOption::AllInit ) {
            size_t size = TS.t_ts.size();
            vector<vector<double>> result;
            double t = 0;
            for (auto tau : f_spec.delays) {
                t = min(t, -tau);
            }
            for (auto tau : b_spec.delays) {
                t = min(t, -tau);
            }
            
            while(t <= t_f) {
                vector<double> v(n+1);
                v[0] = t;
                auto x = TS.interpolation(t);
                VecCopy(x, v, 1);
                result.push_back(v);
                t += 0.0001;
            }
            
            // for (int i = 0; i < size; i++) {
            //     vector<double> v(n+1);
            //     v[0] = TS.t_ts[i];
            //     VecCopy(TS.x_ts[i], v, 1);
            //     // copy(TS.x_ts[i].begin(), TS.x_ts[i].end(), v.begin()+1);
            //     result.push_back(v);
            // }
            return result;
        // } else if constexpr (return_option == ReturnOption::DiscontinuityCount) {
        //     return discontinuity_count;
        } else {
            return;
        }
    }
    
    
};