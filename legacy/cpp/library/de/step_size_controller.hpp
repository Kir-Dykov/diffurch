#pragma once

#include "time_series.hpp"

struct ConstantStepSize {
    
    Real h = 0.1;
    
    ConstantStepSize(Real h) : h(h) {};
    
    template <size_t n, size_t phi_derivatives_n>
    inline Real step_size(RK_TimeSeries<n, phi_derivatives_n>& TS) {
        return h;
    }
    
    template <size_t n, size_t phi_derivatives_n>
    inline bool is_to_be_rejected(RK_TimeSeries<n, phi_derivatives_n>& TS) {
        return false;
    }
    
};




struct VariableStepSize {
    
    Real max_step_size = 0.1;
    Real init_step_size = 0.1;
    
    Real atol = 0.001;
    Real rtol = 0.001;
    
    Real current_step_size;
    
    VariableStepSize(Real atol, Real rtol, Real max_step_size=0.1, Real init_step_size=0.1) : 
    atol(atol), rtol(rtol), max_step_size(max_step_size), init_step_size(init_step_size) {};
    
    template <size_t n, size_t phi_derivatives_n>
    inline Real step_size(RK_TimeSeries<n, phi_derivatives_n>& TS) {
        if (TS.t_ts.size() <= 1) {
            current_step_size = init_step_size;
        } if (TS.e_ts.back() == 0) {
            current_step_size *= 1.2;
        } else {
            Real tol = atol + rtol*norm(TS.x_ts.back());
            Real prev_step_size = (TS.t_ts.size() > 1) ? (TS.t_ts.end()[-1] - TS.t_ts.end()[-2]) : init_step_size;
            current_step_size = 0.9 * prev_step_size * pow(tol/TS.e_ts.back(), 1/(Real)(RK_TimeSeries<n, phi_derivatives_n>::rk::p));   
        }
        current_step_size = min(current_step_size, max_step_size);
        return current_step_size;
    }
    
    template <size_t n, size_t phi_derivatives_n>
    inline bool is_to_be_rejected(RK_TimeSeries<n, phi_derivatives_n>& TS) {
        return TS.e_ts.back() > atol + rtol*norm(TS.x_ts.back());
    }  
};

