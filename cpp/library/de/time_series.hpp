#pragma once

#include <vector>
#include "runge_kutta_tables.hpp"
#include "vec.hpp"


using namespace std;


#ifndef RK
#define RK RK98_
#endif


template<size_t n, size_t phi_derivatives_n, 
typename rk = RK>
class RK_TimeSeries {
public:

    // inline static const RK< RK_method::RK4 > rk; // without `inline' rk is not initialized and linker cannot find it
    // inline static const RK<rk_method> rk; // without `inline' rk is not initialized and linker cannot find it
   
    vector<Real> t_ts;
    vector<Vec<n>> x_ts;
    vector<array<Vec<n>, rk::s>> K_ts;
    
    VecMapC<1, n, phi_derivatives_n> phi;
    
    RK_TimeSeries(Real t_0, VecMapC<1, n, phi_derivatives_n> phi) : phi{phi}, t_ts({t_0}), x_ts({phi(t_0)}) {};
    
    template<size_t derivative_order = 0>
    Vec<n> eval(Real t) {
        if (t <= t_ts[0]) {
            return phi.template eval<derivative_order>(t);
        } else {
            // `t_ts[t_i + 1]` is the first element of `t_ts` greater than `t`
            int t_i = distance(t_ts.begin(), upper_bound(t_ts.begin(), t_ts.end(), t)) - 1;
            
            if (t_i + 1 == t_ts.size()) { t_i--; }
            
            // cout << "t : " << t << "    in " << Vec<2>{t_ts[t_i], t_ts[t_i+1]} << endl;

            Real h = t_ts[t_i+1] - t_ts[t_i]; // step size
            Real theta = (t - t_ts[t_i])/h;
            Vec<n> y{}; //initialize to zero
            
            for (int j = 0; j < rk::s; j++) {
                y = y + rk::BB[j].template eval<derivative_order>(theta)*K_ts[t_i][j];
                // y = y + rk::B[j]*theta*K_ts[t_i][j]; // linear interpolation default
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
    Vec<n*(non_delayed + delays_n + neutral_delays_n)> eval_pack(Real t, const ArgSpec<non_delayed, delays_n, neutral_delays_n, derivatives_n>& arg_spec) {
        Vec<n*(non_delayed + delays_n + neutral_delays_n)> result;
        size_t offset = 0;
        if constexpr (non_delayed) {
            VecCopy(eval(t), result, offset);
            offset += n;
        }
        for (size_t i = 0; i < delays_n; i++) {
            VecCopy(eval(t - arg_spec.delays[i](t)), result, offset);
            offset += n;
        }
        for (size_t i = 0; i < neutral_delays_n; i++) {
            VecCopy(eval<1>(t - arg_spec.neutral_delays[i](t)), result, offset);
            offset += n;
        }
        return result;
    }
    
    template<bool non_delayed = true, size_t delays_n = 0, size_t neutral_delays_n = 0, size_t derivatives_n = 0>
    auto eval_pack(Real t, const Vec<n>& x, const ArgSpec<non_delayed, delays_n, neutral_delays_n, derivatives_n>& arg_spec) {
        Vec<n*(non_delayed + delays_n + neutral_delays_n)> result;
        size_t offset = 0;
        if constexpr (non_delayed) {
            VecCopy(x, result, offset);
            offset += n;
        }
        for (size_t i = 0; i < delays_n; i++) {
            VecCopy(eval(t - arg_spec.delays[i](t)), result, offset);
            offset += n;
        }
        for (size_t i = 0; i < neutral_delays_n; i++) {
            VecCopy(eval<1>(t - arg_spec.neutral_delays[i](t)), result, offset);
            offset += n;
        }
        return result;
    }
       
    
    template<bool non_delayed, size_t delays_n, size_t neutral_delays_n, size_t derivatives_n>
    void step_push(Real t, const VecMapC<n*(non_delayed + delays_n + neutral_delays_n), n, derivatives_n>& f, const ArgSpec<non_delayed, delays_n, neutral_delays_n, derivatives_n>& f_spec) {
    
        const Real h = t - t_ts.back();
        array<Vec<n>, rk::s> K;

        for (int i = 0; i < rk::s; i++) {
            Vec<n> Ksum = dot(rk::A[i], K, i);
            K[i] = f(eval_pack(t_ts.back() + rk::C[i]*h, x_ts.back() + h*Ksum, f_spec));
            
        }
        
        Vec<n> Ksum = dot(rk::B, K, rk::s);

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