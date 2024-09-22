#pragma once

#include "discoque.hpp"
#include "vec.hpp"
#include <math.h>
#include <limits>
#include "runge_kutta_tables.hpp"

#define sNAN std::numeric_limits<double>::signaling_NaN()


#define DEBUG

#ifdef DEBUG
// if this ever happens, then binary search is to be implemented instead
#define DEBUG_DISCONTINUITY if (!(t <= prev_t + h && t >= prev_t)) {\
                    cout << "File:" << __FILE__ << ";    Warning: discontinuity time = " << t \
                    << "is outside of the range [prev_t, prev_t + h] = " << array<double,2>{prev_t, prev_t + h} << endl;  }
#endif

using namespace std;


double eval_polynomial(double theta, const vector<double>& coefs) {
    // coefs[i] is the coefficient for theta^(i)
    int deg = coefs.size() - 1;
    double res = coefs[deg];
    for (int j = deg - 1; j >= 0; j--) {
        res *= theta;
        res += coefs[j];
    }
    return res;
}
double eval_polynomial_derivative(double theta, const vector<double>& coefs) {
    // coefs[i] is the coefficient for theta^(i)
    int deg = coefs.size() - 1;
    double res = deg * coefs[deg];
    for (int j = deg - 1; j > 0; j--) {
        res *= theta;
        res += j*coefs[j];
    }
    return res;
}





/*

What I want it to look like:

we define
tau : either vector or a single double, or, perhaps, some bad value to sign an absence of delays

[f_0, f_1, f_2, ..., f_l]
[f_0, f_1, f_2, ..., f_l]

b : R^n_tau -> R
b_derivative
b_levels = []

ODE:
    DDE(f, df)
Smooth DDE:
    DDE(tau or [tau_i], f, df)
Non-smooth DDE with one discontinuity function:
    DDE(tau or [tau_i], [f_i] of size b_n, [df_i] of size b_n, b, b_derivative, b_n)
Non-smooth DDE with multiple discontinuity functions:
    here either we have a sum of Heaviside functions or IDK there are too many f's then


b : if n_b_levels > 0
    from definition of b m_tau can be determined
n_b_levels
f_i : R^(n(1+m_tau)) -> R^ n
    from that n and m_tau can be determined

DDE<n>(f, df, [tau]) --- no b
DDE<n>(f, df)
*/



/*
structure:
    DDE 
        contains
            f,df,b,b',tau
        uses none
    IVP
        inherits
            DDE
            RK_Table
        contains
            phi,phi',psi_i
            ts,Ks,Xs,disco,Ws,WKs,LEs -> for each solution it is unique
        functions:
            interpolation(t, ts, Xs, Ks), interpolation_derivative(t, ts, Xs, Ks)
                used both for x and w
            Z_tau

for each solution Ks, Xs, etc need to be global across functions, separate IVP class exists
*/




// template <int n, int n_tau, int m_tau>
// class IVP : public DDE<n, n_tau, m_tau>;

// const double Nan = numeric_limits<double>::quiet_NaN();

// n - order of the system
// delays - number of delays in the system
// true_delays - number of delays such that delayed argument appears in f0 or f1
template <int n, int n_tau, int m_tau = n_tau, int b_levels = 0>
class DDE {
public:
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    // discontinuity surface function
    Func_b b; 
    Func_db b_derivative; 
    Func_f f0;
    Func_f f1;
    
    Func_df df0;
    Func_df df1;
    
    // delay values them-selfs
    array<double, n_tau> tau;
     
    DDE( Func_b b, Func_db b_derivative, Func_f f0, Func_f f1, Func_df df0, Func_df df1, array<double, n_tau> tau) : b(b), b_derivative(b_derivative), f0(f0), f1(f1), df0(df0), df1(df1), tau(tau) {};
        
    DDE( Func_b b, Func_db b_derivative, Func_f f0, Func_f f1, Func_df df0, Func_df df1, double tau) : b(b), b_derivative(b_derivative), f0(f0), f1(f1), df0(df0), df1(df1), tau({tau}) {};
        
    DDE(const DDE<n, n_tau, m_tau>& other) : b(other.b), b_derivative(other.b_derivative), f0(other.f0), f1(other.f1), df0(other.df0), df1(other.df1), tau(other.tau) {
        std::cout << "Copy constructor called" << std::endl;
    }
    
};


template <int n, int n_tau, int m_tau>
class InitialValueProblem : public DDE<n, n_tau, m_tau>, 
                            public RK<RK5_4_7FM> {
public:
    using DDE<n, n_tau, m_tau>::tau;
    using DDE<n, n_tau, m_tau>::b;
    using DDE<n, n_tau, m_tau>::b_derivative;
    using DDE<n, n_tau, m_tau>::f0;
    using DDE<n, n_tau, m_tau>::f1;
    using DDE<n, n_tau, m_tau>::df0;
    using DDE<n, n_tau, m_tau>::df1;
                                
    static const int n_arg_f = n*(1 + m_tau);
    using Func_b = Func<n*n_tau, 1>;
    using Func_db = Func<n*n_tau, n*n_tau>;
    using Func_f = Func<n_arg_f, n>;
    using Func_df= Func2<n_arg_f,n_arg_f, n>;
    
    Func_f   f;
    Func_df df;   
    
    double max_tau;
    /*
    #######################
    #   Initial data   #
    #######################
    */
    Func<1, n> phi; // initial function defined for t <= 0
    Func<1, n> phi_derivative; // initial function defined for t <= 0

    Func<1, n> psi; // initial function defined for t <= 0

    /*
    #######################
    #   Solution data     #
    #######################
    */
    vector<double> ts;
    vector<Vec<n>> Xs;
    vector<array<Vec<n>, s>> Ks;
    DiscoQue<n_tau, 6> disco;
                         
    /*
    ###############################
    #   Lyapunov exponents data   #
    ##############################
    */
    vector<double> LEs; // finite time lyapunov exponent (only one for now)
    vector<Vec<n>> Ws;
    vector<array<Vec<n>, s>> WKs;
                                
                                
    
    /* Constructor */
    InitialValueProblem(DDE<n, n_tau, m_tau> dde, Func<1, n> phi, Func<1, n> phi_derivative) : DDE<n, n_tau, m_tau>(dde), phi(phi), phi_derivative(phi_derivative), disco(dde.tau) {
        if constexpr (n_tau == 1) { max_tau = tau[0]; }
        else {max_tau = *max_element(tau.begin(), tau.end());}
    };

                                
     /*
    ###########################################
    #   Interpolation at delayed arguments    #
    ###########################################
    */                            
    // interpolation is given by the formula
    // x( t_i(1- \theta) + t_{i+1} \theta ) =
    //    = x_i + theta h \sum_j b_j(theta) K_j
    Vec<n> interpolation(double t, const vector<double>& ts, const vector<Vec<n>>& Xs, const vector<array<Vec<n>, s>>& Ks, const Func<1, n>& phi) const {    
        if (t <= 0) {
            return phi(t);
        } 
        else {
            // upper_bound can be optimized here:
            //      instead of binary search, where at each iteration the array splits in half
            //      we use binary search, where array splits in ratios (t - minimum(ts)) : (maximum(ts) - t)
            //      abusing that values in ts are almost an arithmetical progression
            int t_i = distance(ts.begin(), upper_bound(ts.begin(), ts.end(), t)) - 1;
            if (t_i < 0) {
                cout << "t_i < 0!!!" << endl; //debug
            }
            double h = ts[t_i+1] - ts[t_i];
            double theta = (t - ts[t_i])/h;
            Vec<n> y{}; //initialize to zero
            for (int i = 0; i < s; i++) {
                double bi = eval_polynomial(theta, B_polynomial_coefs[i]);
                y = y + bi*Ks[t_i][i];
            }
            return Xs[t_i] + theta*h*y;
        }
    };
    
    Vec<n> interpolation_x(double t) const {
        return interpolation(t, ts, Xs, Ks, phi);
    }
    
    Vec<n> interpolation_w(double t) const {
        return interpolation(t, ts, Ws, WKs, psi);
    }

    Vec<n> interpolation_derivative(double t, const vector<double>& ts, const vector<Vec<n>>& Xs, const vector<array<Vec<n>, s>>& Ks, const Func<1, n>& phi_derivative) const {      // dx/dt = dx/d\theta * d\theta / dt
        if (t <= 0) {
            return phi_derivative(t); // t_tau < 0 // INITIAL INTERVAL ACCESSED
        } 
        else {
            int t_i = distance(ts.begin(), upper_bound(ts.begin(), ts.end(), t)) - 1;
            if (t_i < 0) {
                cout << "t_i < 0!!!" << endl; //debug
            }
            double h = ts[t_i+1] - ts[t_i];
            double theta = (t - ts[t_i])/h;
            Vec<n> y{};
            for (int i = 0; i < s; i++) {
                double bi = eval_polynomial(theta, B_polynomial_coefs[i]);
                double bi_derivative = eval_polynomial_derivative(theta, B_polynomial_coefs[i]);
                // CAN BE OPTIMIZED
                y = y + (bi + theta*bi_derivative)*Ks[t_i][i];
            }
            return y;
        }
    };
                                
    Vec<n> interpolation_derivative_x(double t) const {
        return interpolation_derivative(t, ts, Xs, Ks, phi_derivative);
    }
                                
                                
    /*
    #########################
    #   Delayed argumetns   #
    #########################
    */                                  
    template <int k_tau>
    Vec<n*k_tau> Z_tau(double t, const vector<double>& ts, const vector<Vec<n>>& Xs, const vector<array<Vec<n>, s>>& Ks, const Func<1, n>& phi) const {
        Vec<n*n_tau> argument;
        
        if constexpr (n == 1 && k_tau == 1) { 
            argument = interpolation(t - tau[0], ts, Xs, Ks, phi); 
        } else if constexpr (n == 1) { // n_tau > 1
            for (int j = 0; j < k_tau; j++) {
                argument[j] = interpolation(t - tau[j], ts, Xs, Ks, phi); 
            }
        } else {
            for (int j = 0; j < k_tau; j++) {
                Vec<n> x_tau_j = interpolation(t - tau[j], ts, Xs, Ks, phi);
                copy(x_tau_j.begin(), x_tau_j.end(), argument.begin() + j*n);
            }
        }

        return argument;
    }
                                
    Vec<n*n_tau> x_tau(double t) const {
        return Z_tau<n_tau>(t, ts, Xs, Ks, phi);
    }
                                
    Vec<n*n_tau> w_tau(double t) const {
        return Z_tau<n_tau>(t, ts, Ws, WKs, psi);
    }                                
                                
    Vec<n*n_tau> dot_x_tau(double t) const {
        Vec<n*n_tau> derivative;
        
        if constexpr (n == 1 && n_tau == 1) { 
            derivative = interpolation_derivative_x(t - tau[0]); 
        } else if constexpr (n == 1) { // n_tau > 1
            for (int j = 0; j < n_tau; j++) {
                derivative[j] = interpolation_derivative_x(t - tau[j]); 
            }
        } else {
            for (int j = 0; j < n_tau; j++) {
                Vec<n> dot_x_tau_j = interpolation_derivative_x(t - tau[j]);
                copy(dot_x_tau_j.begin(), dot_x_tau_j.end(), derivative.begin() + j*n);
            }
        }
        return derivative;
    }
    
    template <int N>
    Vec<n*(1+m_tau)> f_arg(const Vec<n>& X, const Vec<N>& X_tau) const {
        if constexpr (m_tau == 0) {
            return X;
        } else {
            static_assert(N >= n*m_tau, "Error: X_tau argument in f_arg is too small!");
            if constexpr (n == 1 && m_tau == 1) { // scalar case
                return {X, X_tau};
            } else {
                Vec<n*(1+m_tau)> X_f;
                if constexpr (n == 1) { // scalar case
                    X_f[0] = X;
                    for (int j = 0; j < m_tau; j++) {
                        X_f[j+1] = X_tau[j];
                    }
                } else { // vector case
                    copy(X.begin(), X.end(), X_f.begin());
                    copy(X_tau.begin(), X_tau.begin()+n*m_tau, X_f.begin() + n);
                }

                return X_f;
            }
            
        }
    }     
                                                   
                                
    // Vec<n*(1+m_tau)> f_arg(const Vec<n>& X, const Vec<1>& X_tau) const {
    //     if constexpr (m_tau == 0) {
    //         return X;
    //     } else {
    //         static_assert(n == 1 && m_tau==1, "Error: X_tau argument in f_arg is too small!");
    //         return {X, X_tau};
    //     }
    // }  
    // template Vector<int> Vector<int>::operator+(const Vector<short> & other) const;
                                
    Vec<n*(1+m_tau)> x_f(double t, const Vec<n>& X) const {
        if constexpr (m_tau == 0) {
            return X;
        } else {
            Vec<n*m_tau> X_tau = Z_tau<m_tau>(t, ts, Xs, Ks, phi);
            return f_arg<n*m_tau>(X, X_tau);
        }
    }     
                                
    Vec<n*(1+m_tau)> w_f(double t, const Vec<n>& W) const {
        if constexpr (m_tau == 0) {
            return W;
        } else {
            Vec<n*m_tau> W_tau = Z_tau<m_tau>(t, ts, Ws, WKs, psi);
            return f_arg<n*m_tau>(W, W_tau);
        }
    }  
                                
                                    
    
                                
    /*
    ####################
    #   b evaluations  #
    ####################
    */ 
                                
    double evaluate_beta(double t) const
    {
        return b(x_tau(t));
    }

    double evaluate_beta_derivative(double t) const
    {
        return b_derivative(x_tau(t)) * dot_x_tau(t); //dot product
    }

    double beta_newton_fraction(double t) const
    {
        Vec<n*n_tau> argument = x_tau(t);
        return b(argument)/(b_derivative(argument)*dot_x_tau(t));
    }

    /*
    #######################
    #   RUNGE KUTTA STEP  #
    #######################
    */ 
    void step_x(double t, 
                double t_next,
         const Vec<n>& x, 
               Vec<n>& xNext,
     array<Vec<n>, s>& K) const {    
        double h = t_next - t;

        for (int i = 0; i < s; i++) {
            Vec<n> Ksum{};
            for (int j = 0; j < i; j++) {
                Ksum = Ksum + A[i][j]*K[j];
            }
            K[i] = f(x_f(t + C[i]*h, x + h*Ksum));
        }
        
        Vec<n> Ksum =  B[0]*K[0];
        for (int i = 1; i < s; i++) {
            Ksum = Ksum + B[i]*K[i];
        }
        xNext = x + h*Ksum;
    }
                                
    void step_xw(double t, 
                double t_next,
         const Vec<n>& x, 
               Vec<n>& xNext,
     array<Vec<n>, s>& K,
         const Vec<n>& w, 
               Vec<n>& wNext,
     array<Vec<n>, s>& WK) const {  
        
        
        double h = t_next - t;

        for (int i = 0; i < s; i++) {
            Vec<n> Ksum{};
            Vec<n> WKsum{};
            for (int j = 0; j < i; j++) {
                 Ksum = Ksum  + A[i][j]*K[j];
                WKsum = WKsum + A[i][j]*WK[j];
            }
            Vec<n*(1 + m_tau)> X_f = x_f(t + C[i]*h, x + h*Ksum);
            Vec<n*(1 + m_tau)> W_f = w_f(t + C[i]*h, w + h*WKsum);
            K[i] = f(X_f);
            WK[i] = df(X_f, W_f);
        }
        
        Vec<n> Ksum =  B[0]*K[0];
        Vec<n> WKsum =  B[0]*WK[0];
        for (int i = 1; i < s; i++) {
            Ksum = Ksum + B[i]*K[i];
            WKsum = WKsum + B[i]*WK[i];
        }
        xNext = x + h*Ksum;
        wNext = w + h*WKsum;
    }

    
                                
    void jump_w(double t, const Vec<n>& X, const Vec<n>& W, Vec<n>& Wnext) const {        
        Vec<n*n_tau> X_tau = x_tau(t);
        Vec<n*n_tau> W_tau = w_tau(t);
        Vec<n*n_tau> dot_X_tau = dot_x_tau(t);
        Vec<n*(1 + m_tau)> X_f = f_arg<n*n_tau>(X, X_tau);
        auto b_der = b_derivative(X_tau);
        Wnext = W + (f1(X_f) - f0(X_f)) * (b_der * W_tau)/(abs(b_der * dot_X_tau));
    }

    

    void compute_solution(double h, double t_finish) {
        double prev_t, t = 0;
        Vec<n> X = phi(0);   
        array<Vec<n>, s> K;
        double b_val_prev, b_val = evaluate_beta(0);
        if (b_val < 0) { f = f0; }
        else           { f = f1; }
        Xs.push_back(X);
        ts.push_back(t);
        auto push_all = [&](){ ts.push_back(t); Xs.push_back(X); Ks.push_back(K); };
        
        while (t < t_finish) {
            prev_t = t;
            b_val_prev = b_val;
            
            t = t + h;
			if (t >= disco.min) {	t = disco.min; 	disco.pop();	}

            b_val = evaluate_beta(t);
            if (b_val*b_val_prev <= 0. && b_val != 0) {     
                /* REFINE DISCONTINUITY */
                t = 0.5*(t+prev_t);
                for (int i = 0; i < 20; i++) {  t -= beta_newton_fraction(t); }
                if (t < prev_t-0.000000001 || t > prev_t + h + 0.000000001) { // then binary search
                    double r = prev_t;          double l = prev_t + h;
                    for (int i = 0; i < 50; i++) {
                        t = (r+l)*0.5;
                        if (evaluate_beta(t)*b_val <= 0) { l = t;}  else {r = t;}
                    }
                }
                DEBUG_DISCONTINUITY
                disco.push(t); // add propagated points to mesh in the future

                step_x(prev_t, t, X, X, K); push_all();
                
                prev_t = t;
                t = prev_t + h;
                
                if (b_val < 0) { f = f0; }
                else           { f = f1; }
                step_x(prev_t, t, X, X, K); push_all();
                
            } 
            else {
                step_x(prev_t, t, X, X, K); push_all();
            }

        }        
    };
                                
    
    void compute_solution_and_lyapunov_exponent(double h, double t_finish, Func<1, n> psi_) {
        double prev_t, t = 0;
        
        Vec<n> X = phi(0);
        array<Vec<n>, s> K;
        
        psi = psi_;
        Vec<n> W = psi(0);
        array<Vec<n>, s> WK;
        
        double LEsum = 0;
        double LE = 0;

        double b_val_prev, b_val = evaluate_beta(0);

        if (b_val < 0) { f = f0; df = df0; }
        else           { f = f1; df = df1; }

        Xs.push_back(X); ts.push_back(t); Ws.push_back(W); LEs.push_back(LE);
        
        auto push_all = [&](){ ts.push_back(t); Xs.push_back(X); Ks.push_back(K); Ws.push_back(W); WKs.push_back(WK); LEs.push_back((LEsum + 0.5*log(W*W))/t);  };
                              
        while (t < t_finish) {
            
            prev_t = t;
            b_val_prev = b_val;
            
            t = t + h;
			if (t >= disco.min) { t = disco.min; disco.pop(); }
            b_val = evaluate_beta(t);
            if (b_val*b_val_prev <= 0. && b_val != 0) {                  
                t = 0.5*(t+prev_t);
                for (int i = 0; i < 20; i++) {  t -= beta_newton_fraction(t); }  
                  if (t < prev_t-0.000000001 || t > prev_t + h + 0.000000001) { // then binary search
                    double r = prev_t;          double l = prev_t + h;
                    for (int i = 0; i < 50; i++) {
                        t = (r+l)*0.5;
                        if (evaluate_beta(t)*b_val <= 0) { l = t;}  else {r = t;}
                    }
                }
                DEBUG_DISCONTINUITY; 
                disco.push(t);

                step_xw(prev_t, t, X, X, K, W, W, WK); 
                
                push_all();
                
                jump_w(t, X, W, W); 
                push_all(); 

                prev_t = t;
                t = prev_t + h;
                
                if (b_val < 0) { f = f0; df = df0; }
                else           { f = f1; df = df1; }

                step_xw(prev_t, t, X, X, K, W, W, WK); 
                push_all();
            } 
            else {                
                step_xw(prev_t, t, X, X, K, W, W, WK);      
                push_all();
            }
            
            // cout << W << endl;
            // Renormalization whooo
            if (0.5*log(W*W) > log(100) || 0.5*log(W*W + WK[0]*WK[0]) < -log(100)) { // if either ||w(t)|| > 10^{3} or ||w(t)|| < 10^{-3}
                // cout << "Renormalization" << endl;
                // cout << W << endl;
                
                double norm = sqrt(W*W);
                double scale = 1/norm;
                
                LEsum += log(norm);
                                
                W = W*scale;
                
                int ti = ts.size()-1;    
                
                // cout << WK << endl;
                
                while (ti >= 0 && ts[ti] >= t - max_tau) {
                    for (int j = 0; j < s; j++) {
                        WKs[ti][j] = WKs[ti][j]*scale;
                    }
                    ti--;
                }
            }

        }        
    };
                                
                                
    vector<vector<double>> get_solution() const {
        vector<vector<double>> output(ts.size(), vector<double>(1 + n));
        for (int i = 0; i < ts.size(); i++) {
            output[i][0] = ts[i];
            if constexpr (n == 1) {
                output[i][1] = Xs[i];
            } else {
                copy(Xs[i].begin(), Xs[i].end(), output[i].begin()+1);
            }
        }
        
        return output;
    }
                                
    vector<vector<double>> get_solution_and_lyapunov_exponent() const {
        vector<vector<double>> output(ts.size(), vector<double>(1 + 1 + n));
        for (int i = 0; i < ts.size(); i++) {
            output[i][0] = ts[i];
            output[i][1] = LEs[i];
            if constexpr (n == 1) {
                output[i][2] = Xs[i];
            } else {
                copy(Xs[i].begin(), Xs[i].end(), output[i].begin()+2);
            }
        }
        
        return output;
    }

};

// template <int n, int n_tau, int m_tau>
// InitialValueProblem<n, n_tau, m_tau> IVP (DDE<n, n_tau, m_tau>& dde, 
//      function<array<double, n>(double)> phi,
//     function<array<double, n>(double)> phi_derivative) 
// {
//     return IVP<n, n_tau, m_tau> (dde, phi, phi_derivative);
// }


