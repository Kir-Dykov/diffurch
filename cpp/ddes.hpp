#pragma once

#include "discoque.hpp"
#include "array_operations.hpp"
#include <math.h>
#include <limits>
#include "runge_kutta_tables.hpp"

#define sNAN std::numeric_limits<double>::signaling_NaN()


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



template <int n, int delays_number, int true_delays_number>
class RelayDD;

// template <int n, int delays_number, int true_delays_number>
// class IVP : public RelayDDE<n, delays_number, true_delays_number>;

// const double Nan = numeric_limits<double>::quiet_NaN();

// n - order of the system
// delays - number of delays in the system
// true_delays - number of delays such that delayed argument appears in f0 or f1
template <int n, int delays_number, int true_delays_number>
class RelayDDE {
public:
    const int disco_que_order = 5;
    
    // discontinuity surface function
    function<double(const array<double, n * delays_number>&)> b; 
    function<array<double, n * delays_number>(const array<double, n * delays_number>&)> b_derivative; 
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f0;
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f1;
    
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&, const array<double, n*(true_delays_number + 1)>&)> df0;
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&, const array<double, n*(true_delays_number + 1)>&)> df1;
    
    // delay values them-selfs
    array<double, delays_number> tau;
     
    RelayDDE(
        function<double(const array<double,  n * delays_number>&)> b, 
        function<array<double, n * delays_number>(const array<double, n * delays_number>&)> b_derivative,
        function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f0,
        function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f1,
        array<double, delays_number> tau) : b(b), b_derivative(b_derivative), f0(f0), f1(f1), tau(tau) {};
        
    RelayDDE(const RelayDDE<n, delays_number, true_delays_number>& other) : b(other.b), b_derivative(other.b_derivative), f0(other.f0), f1(other.f1), tau(other.tau) {
        std::cout << "Copy constructor called" << std::endl;
    }
    
        
    void step_RK4(double h, 
        function<array<double, n>(const array<double, n>&)> f,
        const array<double, n>& x, 
        array<double, n>& xNext) {

        auto K1 = f(x);
        auto K2 = f(x + h*0.5*K1);
        auto K3 = f(x + h*0.5*K2);
        auto K4 = f(x +     h*K3);
        
        xNext = x + h * (K1 * (1./6.) + K2 * (1./3.) + K3 * (1./3.) + K4 * (1./6.));
    }
    
        

    
};


template <int n, int delays_number, int true_delays_number>
class InitialValueProblem : public RelayDDE<n, delays_number, true_delays_number>, 
                            public RK<RK4> {
public:
    using RelayDDE<n, delays_number, true_delays_number>::tau;
    using RelayDDE<n, delays_number, true_delays_number>::b;
    using RelayDDE<n, delays_number, true_delays_number>::b_derivative;
    using RelayDDE<n, delays_number, true_delays_number>::f0;
    using RelayDDE<n, delays_number, true_delays_number>::f1;
    // using RelayDDE<n, delays_number, true_delays_number>::RelayDDE;
    
    
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f;
    
    
    /*
    #######################
    #   Initial data   #
    #######################
    */
    function<array<double, n>(double)> phi; // initial function defined for t <= 0
    function<array<double, n>(double)> phi_derivative; // initial function defined for t <= 0


    /*
    #######################
    #   Solution data     #
    #######################
    */
    vector<double> ts;
    vector<array<double, n>> Xs;
    vector<array<array<double, n>, s>> Ks;
    DiscoQue<delays_number, 6> disco;
    
    // Constructor
    InitialValueProblem(RelayDDE<n, delays_number, true_delays_number> dde,
        function<array<double, n>(double)> phi,
        function<array<double, n>(double)> phi_derivative
        ) : RelayDDE<n, delays_number, true_delays_number>(dde), phi(phi), phi_derivative(phi_derivative), disco(dde.tau) {};


    

    // interpolation is given by the formula
    // x( t_i(1- \theta) + t_{i+1} \theta ) =
    //    = x_i + theta h \sum_j b_j(theta) K_j
    array<double, n> interpolation(double t) const {    
        if (t <= 0) {
            return phi(t);
        } 
        else {
            int t_i = distance(ts.begin(), upper_bound(ts.begin(), ts.end(), t)) - 1;
            if (t_i < 0) {
                cout << "t_i < 0!!!" << endl; //debug
            }
            double h = ts[t_i+1] - ts[t_i];
            double theta = (t - ts[t_i])/h;
            array<double, n> y;
            y.fill(0);
            for (int i = 0; i < s; i++) {
                double bi = eval_polynomial(theta, B_polynomial_coefs[i]);
                y = y + bi*Ks[t_i][i];
            }
            return Xs[t_i] + theta*h*y;
        }
    };

    array<double, n> interpolation_derivative(double t) const {    // dx/dt = dx/d\theta * d\theta / dt
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
            array<double, n> y;
            y.fill(0);
            for (int i = 0; i < s; i++) {
                double bi = eval_polynomial(theta, B_polynomial_coefs[i]);
                double bi_derivative = eval_polynomial(theta, B_polynomial_coefs[i]);

                y = y + (bi + theta*bi_derivative)*Ks[t_i][i];
            }
            // cout << "derivative is " << y << endl;
            return y;
        }
    };

    void step(double t, double t_next,
            const array<double, n>& x, array<double, n>& xNext,
            array<array<double, n>, s>& K) const
    {    
        double h = t_next - t;
        // Define each K[i]
        for (int i = 0; i < s; i++) {
            K[i].fill(0.);
            for (int j = 0; j < i; j++) {
                K[i] = K[i] + A[i][j]*K[j];
            }

            if constexpr (true_delays_number == 0) {
                K[i] = f(x + h*K[i]);
                
            }
            else {
                array<double, n*(true_delays_number)> argument;
                auto y = x + h*K[i];
                copy(y.begin(), y.end(), argument.begin());
                for (int j = 0; j < true_delays_number; j++) {
                    y = interpolation(t + C[i]*h - tau[j]); // DELAYED THINGY
                    copy(y.begin(), y.end(), argument.begin() + (j+1)*n);
                }
                K[i] = f(argument);
            } 

        }
        auto y =  B[0]*K[0];
        for (int i = 1; i < s; i++) {
            y = y + B[i]*K[i];
        }
        xNext = x + h*y;
    }

    double evaluate_beta(double t) const
    {
        array<double, n*(delays_number)> argument;

        for (int j = 0; j < delays_number; j++) {
            auto y = interpolation(t - tau[j]); // DELAYED THINGY
            copy(y.begin(), y.end(), argument.begin() + j*n);
        }
        return b(argument);
    }

    double evaluate_beta_derivative(double t) const
    {
        array<double, n*(delays_number)> argument;
        array<double, n*(delays_number)> derivative;

        for (int j = 0; j < delays_number; j++) {
            auto y = interpolation(t - tau[j]); // DELAYED THINGY
            copy(y.begin(), y.end(), argument.begin() + j*n);
        }

        for (int j = 0; j < delays_number; j++) {
            auto y = interpolation_derivative(t - tau[j]); // DELAYED THINGY
            copy(y.begin(), y.end(), derivative.begin() + j*n);
        }

        return b_derivative(argument) * derivative; //dot product
    }

    double beta_newton_fraction(double t) const
    {
        array<double, n*(delays_number)> argument;
        array<double, n*(delays_number)> derivative;

        for (int j = 0; j < delays_number; j++) {
            auto y = interpolation(t - tau[j]); // DELAYED THINGY
            copy(y.begin(), y.end(), argument.begin() + j*n);
        }

        for (int j = 0; j < delays_number; j++) {
            auto y = interpolation_derivative(t - tau[j]); // DELAYED THINGY
            copy(y.begin(), y.end(), derivative.begin() + j*n);
        }

        return b(argument)/(b_derivative(argument)*derivative);
    }


    void solution(double h, double t_finish) {
        // cout << "start frm" << endl;
        // zero_count = 0;
 
        double t = 0;
        array<double, n> X = phi(0);
        array<array<double, n>, s> K;


        double b_val, b_val_prev;
        b_val = b_val_prev = evaluate_beta(0);

        if (b_val < 0) { f = f0; }
        else           { f = f1; }


         Xs.push_back(X);
         ts.push_back(t);

        double prev_t;
        while (t < t_finish) {
        
            prev_t = t;
            
            t = t + h;

			if (t >= disco.min) {
				t = disco.min;
				disco.pop();
			}

            

            b_val = evaluate_beta(t);
            if (b_val*b_val_prev <= 0. && b_val != 0) {                  
                double t_discontinuity = 0.5*(t+prev_t);
                for (int i = 0; i < 20; i++) {
                    double beta = evaluate_beta(t_discontinuity);
                    double d_beta = evaluate_beta_derivative(t_discontinuity);
//                     cout << "Newton iteration:" << endl;
//                     cout << "\tt_discontinuity = " << t_discontinuity << endl;
//                     cout << "\tbeta = " << beta << endl;
//                     cout << "\tbeta' = " << d_beta << endl;
//                     cout << "\tdelta = " << -beta/d_beta << endl;
                    
                    t_discontinuity -= beta_newton_fraction(t_discontinuity);
                }
                
            
                if (!(t_discontinuity <= t && t_discontinuity >= prev_t)) {
                    cout << "Warning: t_discontinuity = " << t_discontinuity 
                    << "is outside of the range [prev_t, t] = " << array<double,2>{prev_t, t} << endl;
                    // if this happens, then binary search is to be implemented instead
                }
                
                // t_discontinuity = min(t, t_discontinuity);
                // t_discontinuity = max(prev_t, t_discontinuity);
                
                disco.push(t_discontinuity);

                step(prev_t, t_discontinuity, X, X, K);
                ts.push_back(t);
                Xs.push_back(X);
                Ks.push_back(K);

                if (b_val < 0) { f = f0; }
                else           { f = f1; }

                step(t_discontinuity, t, X, X, K);
                ts.push_back(t);
                Xs.push_back(X);
                Ks.push_back(K);
            } 
            else {
            
                step(prev_t, t, X, X, K);

                ts.push_back(t);
                Xs.push_back(X);
                Ks.push_back(K);
            }

            b_val_prev = b_val;
        }        
    };

};

// template <int n, int delays_number, int true_delays_number>
// InitialValueProblem<n, delays_number, true_delays_number> IVP (RelayDDE<n, delays_number, true_delays_number>& dde, 
//      function<array<double, n>(double)> phi,
//     function<array<double, n>(double)> phi_derivative) 
// {
//     return IVP<n, delays_number, true_delays_number> (dde, phi, phi_derivative);
// }