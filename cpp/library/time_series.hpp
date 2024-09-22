#pragma once

#include <vector>
#include "runge_kutta_tables.hpp"
#include "vec.hpp"

using namespace std;


template<bool non_delayed = true, size_t delays_n = 0, size_t neutral_delays_n = 0, size_t derivatives_n = 0>
struct ArgSpec {
    array<double, delays_n> delays;
    array<double, neutral_delays_n> neutral_delays;
    
    ArgSpec(array<double, delays_n> delays,
    array<double, neutral_delays_n> neutral_delays) : delays(delays), neutral_delays(neutral_delays) {};
};






template<size_t n, size_t phi_derivatives_n, 
RK_method rk_method = RK_method::RK5_4_7FM>
class RK_TimeSeries {
public:

    // inline static const RK< RK_method::RK4 > rk; // without `inline' rk is not initialized and linker cannot find it
    inline static const RK<rk_method> rk; // without `inline' rk is not initialized and linker cannot find it
   
    vector<double> t_ts;
    vector<Vec<n>> x_ts;
    vector<array<Vec<n>, rk.s>> K_ts;
    
    VecFuncC<n, phi_derivatives_n> phi;
    
    RK_TimeSeries(double t_0, VecFuncC<n, phi_derivatives_n> phi) : phi{phi}, t_ts({t_0}), x_ts({phi(t_0)}) {};
    
    template<size_t derivative_order = 0>
    Vec<n> eval(double t) {
        if (t <= t_ts[0]) {
            return phi.template eval<derivative_order>(t);
        } else {
            // `t_ts[t_i + 1]` is the first element of `t_ts` greater than `t`
            int t_i = distance(t_ts.begin(), upper_bound(t_ts.begin(), t_ts.end(), t)) - 1;
            
            if (t_i + 1 == t_ts.size()) { t_i--; }

            double h = t_ts[t_i+1] - t_ts[t_i]; // step size
            double theta = (t - t_ts[t_i])/h;
            Vec<n> y{}; //initialize to zero
            
            for (int j = 0; j < rk.s; j++) {
                y = y + rk.BB[j].template eval<derivative_order>(theta)*K_ts[t_i][j];
            }
            if constexpr (derivative_order == 0) {
                return x_ts[t_i] + h*y;
            } else {
                for (int i = 0; i < derivative_order - 1; i++) {
                    y = y * (1/h);
                }
                return y;
            }
        }
    }
      
    template<bool non_delayed = true, size_t delays_n = 0, size_t neutral_delays_n = 0, size_t derivatives_n = 0>
    auto eval_pack(double t, const ArgSpec<non_delayed, delays_n, neutral_delays_n, derivatives_n>& arg_spec) {
        Vec<n*(non_delayed + delays_n + neutral_delays_n)> result;
        size_t offset = 0;
        if constexpr (non_delayed) {
            VecCopy(eval(t), result, offset);
            offset += n;
        }
        for (size_t i = 0; i < delays_n; i++) {
            VecCopy(eval(t - arg_spec.delays[i]), result, offset);
            offset += n;
        }
        for (size_t i = 0; i < neutral_delays_n; i++) {
            VecCopy(eval<1>(t - arg_spec.neutral_delays[i]), result, offset);
            offset += n;
        }
    }
    
    template<bool non_delayed = true, size_t delays_n = 0, size_t neutral_delays_n = 0, size_t derivatives_n = 0>
    auto eval_pack(double t, const Vec<n>& x, const ArgSpec<non_delayed, delays_n, neutral_delays_n, derivatives_n>& arg_spec) {
        Vec<n*(non_delayed + delays_n + neutral_delays_n)> result;
        size_t offset = 0;
        if constexpr (non_delayed) {
            VecCopy(x, result, offset);
            offset += n;
        }
        for (size_t i = 0; i < delays_n; i++) {
            VecCopy(interpolation(t - arg_spec.delays[i]), result, offset);
            offset += n;
        }
        for (size_t i = 0; i < neutral_delays_n; i++) {
            VecCopy(interpolation<1>(t - arg_spec.neutral_delays[i]), result, offset);
            offset += n;
        }
    }
       
    
    template<bool non_delayed, size_t delays_n, size_t neutral_delays_n, size_t derivatives_n>
    void step_push(double t, const VecMapC<n*(non_delayed + delays_n + neutral_delays_n), n, derivatives_n>& f, const ArgSpec<non_delayed, delays_n, neutral_delays_n, derivatives_n>& f_spec) {
    
        const double h = t - t_ts.back();
        array<Vec<n>, rk.s> K;

        for (int i = 0; i < rk.s; i++) {
            Vec<n> Ksum = dot(rk.A[i], K, i);
            K[i] = f(delayed_arg(t + rk.C[i]*h, x_ts.back() + h*Ksum, f_spec));
        }
        
        Vec<n> Ksum = dot(rk.B, K, rk.s);

        Vec<n> x = x_ts.back() + h*Ksum;
        
        t_ts.push_back(t);
        x_ts.push_back(x);
        K_ts.push_back(K);
    }
    
    
    
    void pop_back() {
        t_ts.pop_back();
        x_ts.pop_back();
        K_ts.pop_back();
    }

        
};