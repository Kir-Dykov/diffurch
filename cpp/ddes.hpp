#pragma once

#include "discoque.hpp"
#include "array_operations.hpp"
#include <math.h>
#include <limits>

#define sNAN std::numeric_limits<double>::signaling_NaN()


using namespace std;



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
    // RHS for B < 0
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f0;
    // RHS for B < 0
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f1;
    // delay values them-selfs
    
    // array<double, n> (f*)(const array<double, n*(true_delays_number + 1)>&);
    
    array<double, delays_number> tau;
     
    RelayDDE(function<double(const array<double,  n * delays_number>&)> b, 
    function<array<double, n * delays_number>(const array<double, n * delays_number>&)> b_derivative,
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f0,
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f1,
    array<double, delays_number> tau) : b(b), b_derivative(b_derivative), f0(f0), f1(f1), tau(tau) {};
    
        
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
    
    const int s = 4;
    const vector<vector<double>> A = 
        {{}, 
         {1./2.}, 
         {0.,     1./2.}, 
         {0.,     0.,      1.}};
    const vector<double> C =
        {0,
         1./2.,
         1./2.,
         1.};
    const vector<double> B =
        {1./6., 1./3., 1./3., 1./6.};
    const vector<vector<double>> B_polynomial_coefs =
        {{1./6.}, {1./3.}, {1./3.}, {1./6.}}; // linear interpolation coefficients from high to low power
        
//     void interpolant(array<array<double, n>, s> K, double h, double theta) {
//         array<double, n> result;
//         result.fill(0);
//         for (int i = 0; i < s; i++) {
//             result = result + b[i](theta)*K[i];
//         }
        
//     }
            
    class IVP {
        function<array<double, n>(double))> phi; // initial function defined for t <= 0
        vector<double>& ts;
        vector<array<double, n>>& Xs;
        vector<array<array<double, n>, s>>& Ks;
        
        IVP(function<array<double, n>(double))> phi) : phi(phi) {};
        
        array<double, n> interpolation(double t, function<array<double, n>(double))> phi,  const vector<double>& ts, const vector<array<double, n>>& Xs, const vector<array<array<double, n>, s>>& Ks) const {    
        if (t < 0) {
            y = X0; // t_tau < 0 // INITIAL INTERVAL ACCESSED
        } 
        else {
            int t_i = distance(ts.begin(), upper_bound(ts.begin(), ts.end(), t)) - 1;
            if (t_i < 0) {
                cout << "t_i < 0!!!" << endl; //debug
            }
            double h = ts[t_i+1] - ts[t_i];
            double theta = (t_tau - ts[t_i])/h;
            array<double, n> y;
            y.fill(0);
            for (int i = 0; i < s; i++) {
                bi = B_polynomial_coefs[i][0];
                for (int j = 1; j < B_polynomial_coefs[i].size(); j++) {
                    bi *= theta;
                    bi += B_polynomial_coefs[i][j];
                }
                y = y + bi*Ks[t_i][i];
            }
            return Xs[t_i] + theta*h*y
        }
    }
        
    };

    
    
    void step(double t, double t_next,
            const array<double, n>& x, 
            array<double, n>& xNext, array<array<double, n>, s>& K,
            , const vector<double>& ts, const vector<array<double, n>>& Xs, const vector<array<array<double, n>, s>>& Ks) const
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
                    y = interpolation(t + C[i]*h - tau[j], phi, ts, Xs, Ks); // DELAYED THINGY
                    copy(y.begin(), y.end(), argument.begin() + (j+1)*n);
                }
                K[i] = f(argument);
            } 

        }
        xNext.fill(0.);
        for (int i = 0; i < s; i++) {
            xNext = xNext + b[i]*K[i];
        }
        xNext = x + h*xNext;
    }
    
    double evaluate_b(double t,
            function<array<double, n>(double))> phi, 
            const vector<double>& ts, const vector<array<double, n>>& Xs, const vector<array<array<double, n>, s>>& Ks) const
    {
        array<double, n*(delays_number)> argument;
        
        for (int j = 0; j < delays_number; j++) {
            y = interpolation(t - tau[j], phi, ts, Xs, Ks); // DELAYED THINGY
            copy(y.begin(), y.end(), argument.begin() + j*n);
        }
        return b(argument);
    }
    
    void solution(function<array<double, n>(double))> phi, double h, double t_finish, vector<double>& ts, vector<array<double, n>>& Xs) {
        // cout << "start frm" << endl;
        // zero_count = 0;
        array<double, n> X0 = phi(0);
        
        double t = 0;
        array<double, n> X = X0;
        array<array<double, n>, s>& K;
        
        function<array<double, n>(const array<double, n*(1+true_delays_number)>&)> f;
        
        array<double, n*(delays_number)> phi; 
        for (int i = 0; i < delays_number; i++) {
            copy(X0.begin(), X0.end(), phi.begin()+i*n);
        }
        
        if (b(X0) < 0) { f = f0; }
        else           { f = f1; }
        
        
         Xs.push_back(X);
		 ts.push_back(t);
         
         vector<array<array<double, n>, s>>& Ks;
        
		DiscoQue<delays_number, 6, int> discoque(tau);
                
		double prev_t0 = 0;
		bool sign_switch_flag = false;
        // double dw_step = 0;
        
        
        
        

        
        while (t < t_finish) {
			auto prev_X = X;
			double prev_t = t;

			t += h;

// 			if (sign_switch_flag) {
// 				x_tau = -x_tau;
// 				sign_switch_flag = false;
// 			}

// 			if (t >= disco.top[0]) {
// 				t = disco.top[0];
// 				if (disco.top_order == 0) { 
//                     sign_switch_flag = true; 
//                 }
// 				disco.pop();
// 			}
            
                       
            step(prev_t, t, X, X, K, phi, ts, Xs, Ks);
            
            ts.push_back(t);
            Xs.push_back(X);
            Ks.push_back(K);
            
// 			if (x*prev_x <= 0. && x != 0.) {
				
// 				double x0 = x, dx0 = dx;
// 				double h0 = (t - prev_t);
                
//                 if (10* abs(x0) < abs(dx0)) { // if x0 is small and derivative is not small
//                     h0 -= x0/dx0;
//                     for (int newton_i = 0; newton_i < 20; newton_i++) {
//                         step(prev_x, prev_dx, x_tau, h0, x0, dx0);
//                         h0 -= x0/dx0;
// 				    }
//                 } else {
//                     double h_l = 0;
//                     double h_r = t - prev_t;
//                     double h_m, x_m, dx_m;
//                     for (int iii=0; iii < 50; iii++) {
//                         h_m = (h_l+h_r)*0.5;
//                         step(prev_x, prev_dx, x_tau, h_m, x_m, dx_m);
//                         if (x_m*prev_x <= 0) h_r = h_m;
//                         else                 h_l = h_m;
//                     }
//                     h0 = h_m;
//                     x0 = x_m;
//                     dx0 = dx_m;
//                 }
                
//                 double t0 = prev_t + h0;
//                 // cout << t0 << endl;
// 				disco.push({t0});
// 				prev_t0 = t0;
// 			}
//             x_.push_back(x);
//             dx_.push_back(dx);
//             t_.push_back(t);
		}        
    };
    
};

