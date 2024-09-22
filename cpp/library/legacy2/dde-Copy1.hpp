#pragma once

#include <tuple>
#include "discoque.hpp"
#include "types.hpp"
#include "vec.hpp"
#include <math.h>
#include <limits>
#include "runge_kutta_tables.hpp"
#include "utils.hpp"

#define sNAN std::numeric_limits<double>::signaling_NaN()


#define DEBUG

#ifdef DEBUG
// if this ever happens, then binary search is to be implemented instead
#define DEBUG_DISCONTINUITY if (!(t <= prev_t + h && t >= prev_t)) {\
                    cout << "File:" << __FILE__ << ";    Warning: discontinuity time = " << t \
                    << "is outside of the range [prev_t, prev_t + h] = " << array<double,2>{prev_t, prev_t + h} << endl;  }
#endif

using namespace std;








/*
    How I want a solution function to look like.
    
    Depending on the problem, different things can be required to be computed, but
    the whole process is very similar. The following variations are needed:
    * Solution alone or with it's lyapunov exponents
    * At a fixed time interval [0, t_f] or until an event
    * Return only the end value or all intermediate steps
    * Custom Events??? Like if x[0] croses zero from above then dx[0] -> -dx[0]
    
    So we have the following information for the customization of the solving function:
    * How many Lyapunov exponents to compute: currently 0 or 1
    * When to stop the integration: 
        integrate while (! t > t_f )
        integrate while (! g(t) > 0 )
        integrate while (! g(t, X, X_tau, ...) > 0 ) ???
    * What to save and return: last value or all values
    
    solve<size_t LE_n, bool save_all>
    
*/



enum class ReturnOption {
    All,
    Last
};



template <size_t n,          // order of the system
        size_t f_delays_n,   // number of delayed arguments in f
        size_t b_delays_n,   // number of delayed arguments in b
        bool b_x_dependence, // is there a non-delayed argument in b
        size_t b_levels_n>   // number discontinuity surfaces
class DE;


template <size_t n,          // order of the system
        size_t f_delays_n,   // number of delayed arguments in f
        size_t b_delays_n,   // number of delayed arguments in b
        bool b_x_dependence, // is there a non-delayed argument in b
        size_t b_levels_n,   // number discontinuity surfaces
        ReturnOption return_option, 
        size_t le_n = 0,
        RK_method rk_method = RK_method::RK5_4_7FM
>   
class DE_IVP;



template <size_t n,          // order of the system
        size_t f_delays_n,   // number of delayed arguments in f
        size_t b_delays_n,   // number of delayed arguments in b
        bool b_x_dependence, // is there a non-delayed argument in b
        size_t b_levels_n>   // number of discontinuity surfaces
class DE {
public:

    static constexpr size_t f_arg_n = n*(1 + f_delays_n);
    static constexpr size_t b_arg_n = n*((size_t)b_x_dependence + b_delays_n);
    
    array<VecMap<f_arg_n, n>,           1+b_levels_n> f;
    array<VecMap2<f_arg_n, f_arg_n, n>, 1+b_levels_n> df;
    Vec<f_delays_n> f_tau;
    
    VecScalarFunc<b_arg_n>         b;
    VecMap<b_arg_n, b_arg_n> b_derivative;
    Vec<b_levels_n>          b_levels;
    Vec<b_delays_n>          b_tau;
       
    DE(
        array<VecMap<f_arg_n, n>, 1+b_levels_n>            f, 
        array<VecMap2<f_arg_n, f_arg_n, n>, 1+b_levels_n> df, 
        Vec<f_delays_n> f_tau, 
        VecScalarFunc<b_arg_n> b, 
        VecMap<b_arg_n, b_arg_n> b_derivative, 
        Vec<b_levels_n> b_levels,
        Vec<b_delays_n> b_tau)
    : f(f), df(df), f_tau(f_tau), b(b), b_derivative(b_derivative), b_levels(b_levels), b_tau(b_tau) {};
      
    template<ReturnOption r_o = ReturnOption::Last>
    auto solution(
             double h, double t_finish,
             const function<Vec<n>(double)>& phi, 
             const function<Vec<n>(double)>& phi_derivative) {
        
        DE_IVP<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n, r_o, 0> 
            ivp (*this, phi, phi_derivative);
        
        return (ivp.solution_function_maker())(h, t_finish);
        
        // if constexpr (r_o == ReturnOption::All) {
        //     return ivp.compute_solution_all(h, t_finish);
        // } else if constexpr (r_o == ReturnOption::Last) {
        //     return ivp.compute_solution_last(h, t_finish);
        // } else {
        //     return;
        // }
    }
};    

/***************************************
        C O N S T R U C T O R S
***************************************/
// Smooth ODE :     ODE( f,   df)
// Non-smooth ODE:  DDE_discontinuous([f], [df],          b, b', [b_levels])
// Smooth DDE :     DDE( f,   df,  [f_tau])
// Non-smooth DDE:  DDE_discontinuous([f], [df], [f_tau], b, b', [b_levels], [b_tau])
// Non-smooth DDE*: DDE_discontinuous_non_delayed_b([f], [df], [f_tau], b, b', [b_levels])
// Relay DDE:       DDE_relay([f], [df],          b, b', [b_levels], [b_tau])

template<size_t n>
auto ODE(VecMap<n, n> f, VecMap2<n, n, n> df) {
    return DE<n, 0, 0, false, 0>({f}, {df}, {}, {}, {}, {}, {});
}


template<size_t n, size_t b_levels_n>
auto ODE_discontinuous(array<VecMap<n, n>, 1+ b_levels_n> f, 
                       array<VecMap2<n, n, n>, 1+ b_levels_n> df, 
                       VecScalarFunc<n> b,
                       VecMap<n, n> b_derivative,
                       Vec<b_levels_n> b_levels) {
    return DE<n, 0, 0, true, b_levels_n>(f, df, {}, b, b_derivative, b_levels, {});
}

template<size_t n, size_t f_delays_n>
auto DDE(VecMap<n, n> f, VecMap2<n, n, n> df, Vec<f_delays_n> f_tau) {
    return DE<n, f_delays_n, 0, false, 0>({f}, {df}, f_tau, {}, {}, {}, {});
}


template<size_t n, size_t f_delays_n, size_t b_delays_n, bool b_x_dependence, size_t b_levels_n>
auto DDE_discontinuous(array<VecMap<n*(1 + f_delays_n), n*(1 + f_delays_n)>, 1+b_levels_n>            f, 
        array<VecMap2<n*(1 + f_delays_n), n*(1 + f_delays_n), n>, 1+b_levels_n> df, 
        Vec<f_delays_n> f_tau, 
        VecScalarFunc<n*((size_t)b_x_dependence + b_delays_n)> b, 
        VecMap<n*((size_t)b_x_dependence + b_delays_n), n*((size_t)b_x_dependence + b_delays_n)> b_derivative, 
        Vec<b_levels_n> b_levels,
        Vec<b_delays_n> b_tau) {
    return DE<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n>(f, df, f_tau, b, b_derivative, b_levels, b_tau);
}

template<size_t n, size_t f_delays_n, size_t b_levels_n>
auto DDE_discontinuous_non_delayed_b(array<VecMap<n*(1 + f_delays_n), n*(1 + f_delays_n)>, 1+b_levels_n>            f, 
        array<VecMap2<n*(1 + f_delays_n), n*(1 + f_delays_n), n>, 1+b_levels_n> df, 
        Vec<f_delays_n> f_tau, 
        VecScalarFunc<n*(1)> b, 
        VecMap<n*(1), n*(1)> b_derivative, 
        Vec<b_levels_n> b_levels) {
    return DE<n, f_delays_n, 0, true, b_levels_n>(f, df, f_tau, b, b_derivative, b_levels, {});
}

// when f does not include the dependence on retareded arguments
template<size_t n, size_t b_delays_n, size_t b_arg_n, size_t b_levels_n>
auto DDE_relay(array<VecMap<n, n>, 1+b_levels_n>            f, 
        array<VecMap2<n, n, n>, 1+b_levels_n> df, 
        VecScalarFunc<(b_arg_n)> b, 
        VecMap<(b_arg_n), (b_arg_n)> b_derivative, 
        Vec<b_levels_n> b_levels,
        Vec<b_delays_n> b_tau) {
    static_assert(b_arg_n % n == 0 && b_arg_n/n >= b_delays_n && b_arg_n/n <= b_delays_n +1, "b-function doesn't have a correct number of arguments");
    return DE<n, 0, b_delays_n, (bool)(b_arg_n/n == b_delays_n + 1), b_levels_n>(f, df, {}, b, b_derivative, b_levels, b_tau);
}














template<typename DE_IVP_Type, bool is_derivative>
struct interpolation_function_maker {};

template<size_t n,          // order of the system
        size_t f_delays_n,   // number of delayed arguments in f
        size_t b_delays_n,   // number of delayed arguments in b
        bool b_x_dependence, // is there a non-delayed argument in b
        size_t b_levels_n,   // number discontinuity surfaces
        ReturnOption r_o, 
        size_t le_n,
        RK_method rk_method,
        bool is_derivative>
struct interpolation_function_maker<
    DE_IVP<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n, r_o, le_n, rk_method>, 
    is_derivative> 
{ 

    static const RK<rk_method> rk;

    using A =  ArgTuple< 
        double, // t
        vector<double>&, // t_ts
        conditional_t<!is_derivative, vector<Vec<n>>&, void>, // x_ts
        vector<array<Vec<n>, rk.s>>&, // K_ts
        VecFunc<n>& //phi or phi_derivative
    >;
    using Arg_t = A::filtered_type;    

    using Ret_t = Vec<n>;

    const function<Ret_t(Arg_t)> func = FuncTupleToArgsTransform([](Arg_t args) {

        auto [t, t_ts, x_ts, K_ts, phi] = A::unpack(args);
        auto& phi_derivative = phi;


        if (t <= 0) {
            if constexpr (is_derivative) {
                return phi_derivative(t);
            } else {
                return phi(t);
            }
        } else {
            // `t_ts[t_i + 1]` is the first element of `t_ts` greater than `t`
            int t_i = distance(t_ts.begin(), upper_bound(t_ts.begin(), t_ts.end(), t)) - 1;
            double h = t_ts[t_i+1] - t_ts[t_i]; // step size
            double theta = (t - t_ts[t_i])/h;
            Vec<n> y{}; //initialize to zero

            if constexpr (is_derivative) {
                for (int i = 0; i < rk.s; i++) {
                    double bi_derivative = rk.BB[i].eval_derivative(theta);
                    y = y + bi_derivative*K_ts[t_i][i];
                }
                return y; // h * SUM_j b_j'(theta)/h * K_j
            } else {
                for (int j = 0; j < rk.s; j++) {
                    double bj = rk.BB[j].eval(theta);
                    y = y + bj*K_ts[t_i][j];
                }
                return x_ts[t_i] + h*y; // x_i + h * SUM_j b_j(theta) * K_j 
            }
        }
    });
};








template<typename DE_IVP_Type, size_t tau_n, bool has_non_retarded_argument, bool is_non_retarded_argument_supplied, bool is_derivative>
struct retarded_argument_function_maker;


template<size_t n,          // order of the system
        size_t f_delays_n,   // number of delayed arguments in f
        size_t b_delays_n,   // number of delayed arguments in b
        bool b_x_dependence, // is there a non-delayed argument in b
        size_t b_levels_n,   // number discontinuity surfaces
        ReturnOption r_o, 
        size_t le_n,
        RK_method rk_method,
        size_t tau_n, bool has_non_retarded_argument, bool is_non_retarded_argument_supplied, bool is_derivative>
struct retarded_argument_function_maker<DE_IVP<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n, r_o, le_n, rk_method>, tau_n, has_non_retarded_argument, is_non_retarded_argument_supplied, is_derivative> {
    static_assert ((has_non_retarded_argument + tau_n) > 0, "retarded argument with zero dimension");
    static_assert (!is_non_retarded_argument_supplied || has_non_retarded_argument, 
                   "Can't supply absent non_retarded_argument");

    static const RK<rk_method> rk;

    using de_ivp_type = DE_IVP<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n, r_o, le_n, rk_method>;

    const TRUE_AUTO_EQ(interpolation,            (interpolation_function_maker<de_ivp_type, false>::func));
    const TRUE_AUTO_EQ(interpolation_derivative, (interpolation_function_maker<de_ivp_type, true >::func));

    using A =  ArgTuple< 
        conditional_t<is_non_retarded_argument_supplied, Vec<n>, void>, // x
        double, // t
        vector<double>&, // t_ts
        conditional_t<!is_derivative, vector<Vec<n>>&, void>, // x_ts
        vector<array<Vec<n>, rk.s>>&, // K_ts
        VecFunc<n>&, //phi or phi_derivative
        conditional_t<tau_n != 0, Vec<tau_n>, void> //tau>,
    >;
    using Arg_t = A::filtered_type;    

    using Ret_t = Vec<n*(has_non_retarded_argument + tau_n)>;

    const function<Ret_t(Arg_t)> func = FuncTupleToArgsTransform([this](Arg_t args) {

        auto [x, t, t_ts, x_ts, K_ts, phi, tau] = A::unpack(args);
        auto& phi_derivative = phi;
        // auto& dx = x;



        Ret_t result;
        auto result_iter = result.begin();  

        if constexpr (has_non_retarded_argument) {
            if constexpr (is_non_retarded_argument_supplied) {
                copy(x.begin(), x.end(), result_iter);
            } else {
                Vec<n> x_tau_j;
                if constexpr (is_derivative) {
                    x_tau_j =            interpolation(t, t_ts, x_ts, K_ts, phi);
                } else {
                    x_tau_j = interpolation_derivative(t, t_ts,       K_ts, phi_derivative);
                }
                copy(x_tau_j.begin(), x_tau_j.end(), result_iter);
            }
            result_iter += n;
        }

        for (int j = 0; j < tau_n; j++) {
            Vec<n> x_tau_j;
            if constexpr (is_derivative) {
                x_tau_j = interpolation(t - tau[j], t_ts, x_ts, K_ts, phi);
            } else {
                x_tau_j = interpolation_derivative(t - tau[j], t_ts,       K_ts, phi_derivative);
            }
            copy(x_tau_j.begin(), x_tau_j.end(), result_iter);
            result_iter += n;
        }
        return result;
    });
};



















template <size_t n,          // order of the system
        size_t f_delays_n,   // number of delayed arguments in f
        size_t b_delays_n,   // number of delayed arguments in b
        bool b_x_dependence, // is there a non-delayed argument in b
        size_t b_levels_n,   // number discontinuity surfaces
        ReturnOption r_o, 
        size_t le_n,
        RK_method rk_method
>   
class DE_IVP :  public DE<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n> {
private:
    /////////////////////   F I E L D S   //////////////////////////////
    
    //#######################
    //#   Equation data     #
    //#######################
    using de_ivp_type = DE_IVP<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n, r_o, le_n, rk_method>;
    using DE_ = DE<n, f_delays_n, b_delays_n, b_x_dependence, b_levels_n>;
    using DE_::f_arg_n;
    using DE_::b_arg_n;
    using DE_::f;
    using DE_::df;
    using DE_::f_tau;
    using DE_::b;
    using DE_::b_derivative;
    using DE_::b_levels;
    using DE_::b_tau;
    //#   RK method data   #
    static const RK<rk_method> rk;
    
    size_t f_i = 0; // index of currently active branch of f, that is f[f_i]
    
    //#######################
    //#   Solution data     #
    //#######################
    VecFunc<n> phi;
    VecFunc<n> phi_derivative;
    vector<double> t_ts; // t-(time-series)
    vector<Vec<n>> x_ts; // x time-series
    vector<array<Vec<n>, rk.s>> K_ts; // k time-series
    
    array<double, f_delays_n + b_delays_n> disco_tau;
    DiscoQue<f_delays_n + b_delays_n, 6> disco; // Discontinuities Queue
    
    // //#######################
    // //#     LE's data       #
    // //#######################    
    // array<VecFunc<n>, le_n> psi;
    // vector<Vec<n>> w_ts; // w time-series
    // vector<array<Vec<n>, rk.s>> k_w_ts; // k time-series for w
    // array<vector<double>, le_n> le_ts;
    
public:
    ////////////////////     C O N S T R U C T O R S     ///////////////////////
    DE_IVP(const DE_& dde, 
           const function<Vec<n>(double)>& phi, 
             const function<Vec<n>(double)>& phi_derivative
           ) : DE_(dde), phi(phi), phi_derivative(phi_derivative), disco(concatenate(f_tau, b_tau)) {
        static_assert(le_n==0, "Error: Initial functions for Lyapunov exponent calculation were not specified.");
    };
    
    
private:
    
    DE_IVP() {}; // Private constructor for use in function-makers
    
    ////////////////////   A U X I L I A R Y   M E T H O D S   ///////////////////////

    //###########################################
    //#   Interpolation at delayed arguments    #
    //###########################################    
    
    
    
    // const decltype(interpolation_function_maker<false>())  interpolation              = interpolation_function_maker<false>();
    // const decltype(interpolation_function_maker<true >())  interpolation_derivative   = interpolation_function_maker<true>();

    // const function<Vec<n>(double, vector<double>&, vector<Vec<n>>&, vector<array<Vec<n>, rk.s>>&, VecFunc<n>&)> interpolation = interpolation_function_maker<false>();
    // const function<Vec<n>(double, vector<double>&, vector<array<Vec<n>, rk.s>>&, VecFunc<n>&)> interpolation_derivative = interpolation_function_maker<true>();
    
    
    
    
    const TRUE_AUTO_EQ(x_in_b,      (retarded_argument_function_maker<de_ivp_type, b_delays_n, b_x_dependence, false, false>::func));
    const TRUE_AUTO_EQ(dot_x_in_b,  (retarded_argument_function_maker<de_ivp_type, b_delays_n, b_x_dependence, false,  true>::func));
    const TRUE_AUTO_EQ(x_in_f,      (retarded_argument_function_maker<de_ivp_type, f_delays_n, true,            true, false>::func));
    
    const function<double(double)> evaluate_beta =        [this](double t) {
        return b(x_in_b(t, t_ts, x_ts, K_ts, phi, b_tau)); 
    };
    const function<double(double)> beta_newton_fraction = [this](double t){ 
        auto     x_in_b_value =     x_in_b(t, t_ts, x_ts, K_ts, phi,            b_tau); 
        auto dot_x_in_b_value = dot_x_in_b(t, t_ts,       K_ts, phi_derivative, b_tau);
        return b(x_in_b_value) / (b_derivative(x_in_b_value) * dot_x_in_b_value);
    };

    /*
    #######################
    #   RUNGE KUTTA STEP  #
    #######################
    */ 
    void step_x(double t, 
                double t_next,
         const Vec<n>& x, 
               Vec<n>& xNext,
     array<Vec<n>, rk.s>& K) const {  
           
        const double h = t_next - t;

        for (int i = 0; i < rk.s; i++) {
            // sum(rk.A[i][j]*K[j], j, 0, i)
            
            Vec<n> Ksum = dot(rk.A[i], K, i);
            // for (int j = 0; j < i; j++) {
                // Ksum = Ksum + rk.A[i][j]*K[j];
            // }
            
            if constexpr (f_delays_n == 0) 
                 {  K[i] = f[f_i](x + h*Ksum);  } 
            else {  K[i] = f[f_i](x_in_f(x + h*Ksum, t + rk.C[i]*h, t_ts, x_ts, K_ts, phi, f_tau) ); }
        }
        
        Vec<n> Ksum = dot(rk.B, K, rk.s);
        
        // Vec<n> Ksum =  rk.B[0]*K[0];
        // for (int i = 1; i < rk.s; i++) {
            // Ksum = Ksum + rk.B[i]*K[i];
        // }
        xNext = x + h*Ksum;
    }
public:             
                                    
    // template <ReturnOption r_o, bool has_event_function = false>
    auto solution_function_maker() {
        
        // static_assert (!has_event_function, "Event function is not implemented.");
                
        using A =  ArgTuple<
            double, // h
            double, // t_finish
            conditional_t<false, VecScalarFunc<n>&, void> // event_function
        >;
        using Arg_t = A::filtered_type;    
        
        using Ret_t = switch_case_t<r_o, void,
            case_pair<ReturnOption::All, tuple<vector<double>, vector<Vec<n>>>>,
            case_pair<ReturnOption::Last, Vec<n>>
        >;
        
        
        function<Ret_t(Arg_t)> func = [this](Arg_t args) {
            auto [h, t_finish, event_function] = A::unpack(args);
            
            double prev_t = 0, t = 0;
            Vec<n> X = phi(0);   
            array<Vec<n>, rk.s> K;
            
            DiscoQue<f_delays_n, 6> disco(f_tau); // Discontinuities Queue
            
            
            double b_val_prev;
            double b_val = evaluate_beta(0);
        
            // the least index, such that b_levels[f_i] > b_val
            size_t f_i = distance(b_levels.begin(), upper_bound(b_levels.begin(), b_levels.end(), b_val));
        
                
            x_ts.push_back(X);
            t_ts.push_back(t);
            auto push_all = [&](){ t_ts.push_back(t); x_ts.push_back(X); K_ts.push_back(K); };
            auto pop_all  = [&](){ t_ts.pop_back(); x_ts.pop_back(); K_ts.pop_back(); };

            while (t < t_finish) {
                prev_t = t;
                b_val_prev = b_val;

                t = t + h;
                
                if (t >= t_finish) { t = t_finish; }
                
                if (t >= disco.min) {   t = disco.min; 	disco.pop();   }
                
                
                step_x(prev_t, t, X, X, K); 
                push_all();

                b_val = evaluate_beta(t);
                
                for (size_t b_levels_i = f_i; b_levels_i >= f_i - 1 && b_levels_i >= 0; b_levels_i--) {
                    
                    double b_level = b_levels[b_levels_i];
                    
                    if ((b_val - b_level)*(b_val_prev - b_level) <= 0. && b_val != b_level) {     
                        /* REFINE DISCONTINUITY */
                        double t_discontinuity = 0.5*(t+prev_t);
                        for (int i = 0; i < 20; i++) {  t_discontinuity -= beta_newton_fraction(t_discontinuity); }
                        if (t_discontinuity < prev_t-0.000000001 || t_discontinuity > t + 0.000000001) { 
                            // then binary search
                            double r = prev_t; double l = t;
                            for (int i = 0; i < 50; i++) {
                                t_discontinuity = (r+l)*0.5;
                                if ((b_val - b_level)*(evaluate_beta(t_discontinuity) - b_level) <= 0) 
                                { l = t;}  else {r = t;}
                            }
                        }
                        DEBUG_DISCONTINUITY
                        disco.push(t_discontinuity); // add propagated points to mesh in the future
                        
                        //reverse last step
                        pop_all();
                        X = x_ts.back();
                        
                        //step from prev_t to t_discontinuity
                        step_x(prev_t, t_discontinuity, X, X, K); 
                        push_all();
                        
                        // change of active function branch
                        if (b_val < b_level) { f_i = b_levels_i; }
                        else                 { f_i = b_levels_i + 1; }
                        
                        // if (b_val < b_level) { f_i = f_i - 1; }
                        // else                 { f_i = f_i + 1; }
                        
                        step_x(t_discontinuity, t, X, X, K); 
                        push_all();
                        
                        break; // only one discontinuity surface must be crossed at each step
                    } 
                }
            }      
            
            if        constexpr (r_o == ReturnOption::Last) {
                return X;
            } else if constexpr (r_o == ReturnOption::All) {
                return make_tuple(t_ts, x_ts);
            } else {
                return;
            }
        };
        
        return FuncTupleToArgsTransform(func);
    }

    
    
    // const TRUE_AUTO_EQ(compute_solution_all , solution_function_maker<ReturnOption::All >());
    // const TRUE_AUTO_EQ(compute_solution_last, solution_function_maker<ReturnOption::Last>());
    
    // const function<tuple<vector<double>, vector<Vec<n>>>(double, double)> compute_solution_all  = solution_function_maker<ReturnOption::All>();
    // const function<Vec<n>(double, double)> compute_solution_last = solution_function_maker<ReturnOption::Last>();

};
   
    
    
// NOTE in interpolation function the use of upper_bound can be optimized: instead of binary search, where at each iteration the array splits in half, we use binary search, where array splits in ratios (t - minimum(ts)) : (maximum(ts) - t), abusing that values in ts are almost an arithmetical progression