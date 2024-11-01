#pragma once

#include "return_handler.hpp"

struct ReturnNone{
    
};

struct ReturnLast {
    template <size_t n, size_t phi_derivatives_n>
    auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        return TS.x_ts.back(); 
    };
};

struct ReturnSolution {
    template <size_t n, size_t phi_derivatives_n>
    static auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        return make_tuple(TS.t_ts, TS.x_ts);
    };
};

struct ReturnSolutionAndDerivative {
    template <size_t n, size_t phi_derivatives_n>
    static auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        vector<Vec<n>> dot_x_ts(TS.t_ts.size());
        std::transform(TS.t_ts.begin(), TS.t_ts.end(), dot_x_ts.begin(), [&TS](Real t){return TS.template eval<1>(t);});
        return make_tuple(TS.t_ts, TS.x_ts, dot_x_ts);
    };
};

struct ReturnDenseSolution{
    int factor = 1;
    ReturnDenseSolution(int factor_) : factor(factor_) {};

    template <size_t n, size_t phi_derivatives_n>
    auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        vector<Real> t_dense_ts = linspace(TS.t_ts[0], TS.t_ts.back(), TS.t_ts.size()*factor);
        vector<Vec<n>> x_dense_ts(t_dense_ts.size());
        std::transform(t_dense_ts.begin(), t_dense_ts.end(), x_dense_ts.begin(), [&TS](Real t){return TS.eval(t);});
        return make_tuple(t_dense_ts, x_dense_ts);
    };
};

struct ReturnSolutionAt {
    vector<Real> ts;
    ReturnSolutionAt(const vector<Real> ts_) ts(ts_) {};
    
    template <size_t n, size_t phi_derivatives_n>
    auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        vector<Vec<n>> result(ts.size());
        std::transform(ts.begin(), ts.end(), result.begin(), [&TS](Real t){return TS.eval(t);});
        return return;
    };
}















template <size_t n>
struct ReturnIntersections {
    function<Real(Vec<n>)> section_func;
    vector<Real> levels;
    
    ReturnIntersections(
        const function<Real(Vec<n>)>& sec_f, 
        const vector<Real> levels) : section_func(sec_f), levels(levels) {};
    
    template <size_t phi_derivatives_n>
    auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 

        function<Real(Real)> func = [&TS, this](Real t){
            return section_func(TS.eval(t));
        };
        vector<Real> ts = find_roots(func, TS.t_ts, b_levels);
        vector<Vec<n>> xs(ts.size());
        std::transform(ts.begin(), ts.end(), xs.begin(), [&TS](Real t){return TS.eval(t);});
        
        return make_pair(ts, xs);
    };
};



template <size_t n>
struct ReturnPeriod {
    
    Real window = 0.;
    function<Real(Vec<n>)> section_func;
    vector<Real> levels;
    Real eps;
    
    ReturnPeriod(const function<Real(Vec<n>)>& sec_f, const vector<Real> levels, Real w = 0., Real eps = 1.e-4) : section_func(sec_f), levels(levels), window(w), eps(eps) {};
    
    template <size_t phi_derivatives_n>
    Real Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        
        function<Real(Real)> func = [&TS, this](Real t){
            return section_func(TS.eval(t));
        };
        vector<Real> ts = find_roots(func, TS.t_ts, b_levels);
        vector<Vec<n>> xs(ts.size());
        std::transform(ts.begin(), ts.end(), xs.begin(), [&TS](Real t){return TS.eval(t);});
        
        int I = ts.size()-1;
        for (int i = I-1; i >= 0; i--) {
            Real distance = norm(xs[i] - xs[I])/(1 + norm(xs[i]));
            
            int j = 1;
            while (i - j >= 0 && (ts[i - j] >= ts[i] - window)) {
                distance += abs((ts[i] - ts[i - j]) - (ts[I] - ts[I - j]))/window;
                j++;
            }
            distance /= Real(j);
            
            if (distance < eps) {
                return ts[I] - ts[i];
            }
        }
        
        return 0;
        
    };
};




struct ReturnDebug {
    template <size_t n, size_t phi_derivatives_n>
    auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
        return TS;
    };
    
    int step_counter = 0;
    
    template <size_t n, size_t phi_derivatives_n>
    void event_function(RK_TimeSeries<n, phi_derivatives_n>& TS) {
        step_counter++;
    }
    
    int discontinuity_counter = 0;
      
    template <size_t n, size_t phi_derivatives_n>
    void discontinuity_event_function(RK_TimeSeries<n, phi_derivatives_n>& TS) {
        discontinuity_counter++;
    }
};