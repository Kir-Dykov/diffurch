#pragma once


template <typename Equation>
requires IsEquation<Equation>
struct State {
    
    // the order of Equation
    static constexpr size_t n = order_of_equation<Equation>::value;
    
    double t;
    double prev_t;
    vector<double> t_sequence;
   
    Vec<n> X;
    Vec<n> prev_X;  
    vector<Vec<n>> X_sequence;
    
    // vector<array<Vec<n>, rk::s>> K_sequence;
    
    Equation* eq_ptr;
    
    State(double initial_time, Equation* eq_ptr_) {
        eq_ptr = eq_ptr_;
        
        t = prev_t = initial_time;
        t_sequence.push_pack(t);
        
        X = prev_X = eq_ptr->initial_condition(initial_time);
        X_sequence.push_back(X);
        
        
    }
    
    template <size_t derivative_order = 0>
    auto eval(double t) {
        if (t <= t_sequence[0]) {
            if constexpr (derivative == 0):
                return eq_ptr->initial_condition(t);
            else
                return eq_ptr->initial_condition<derivative_order>(t);
        } else {
        
            // BAAAAD CODE
            
            // need to generalize for non-interpolating rk-methods
            
            
            
            
            // `t_ts[t_i + 1]` is the first element of `t_ts` greater than `t`
            int t_i = distance(t_sequence.begin(), 
                               upper_bound(t_sequence.begin(), 
                                           t_sequence.end(), t)
                              ) - 1;
            
            if (t_i + 1 == t_sequence.size()) { t_i--; }
            
            // cout << "t : " << t << "    in " << Vec<2>{t_ts[t_i], t_ts[t_i+1]} << endl;

            Real h = t_sequence[t_i+1] - t_sequence[t_i]; // step size
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
    
   
};