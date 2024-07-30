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
    function<double(const array<double, n * delays_number>&)> B; 
    // RHS for B < 0
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f0;
    // RHS for B < 0
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f1;
    // delay values them-selfs
    
    // array<double, n> (f*)(const array<double, n*(true_delays_number + 1)>&);
    
    array<double, delays_number> tau;
     
    RelayDDE(function<double(const array<double,  n * delays_number>&)> B, 
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f0,
    function<array<double, n>(const array<double, n*(true_delays_number + 1)>&)> f1,
    array<double, delays_number> tau) : B(B), f0(f0), f1(f1), tau(tau) {};
    
        
    void step_RK4_relay(double h, 
        function<array<double, n>(const array<double, n>&)> f,
        const array<double, n>& x, 
        array<double, n>& xNext) {

        auto K1 = f(x);
        auto K2 = f(x + h*0.5*K1);
        auto K3 = f(x + h*0.5*K2);
        auto K4 = f(x +     h*K3);
        
        xNext = x + h * (K1 * (1./6.) + K2 * (1./3.) + K3 * (1./3.) + K4 * (1./6.));
    }
    
    void step_DP(double h, 
        function<array<double, n>(const array<double, n>&)> f,
        const array<double, n>& x, 
        array<double, n>& xNext) {

        auto K1 = f(x);
        auto K2 = f(x + h*
        (1./5.)*K1);
        auto K3 = f(x + h*(
        (3./40.)*K1 + (9./40.)*K2));
        auto K4 = f(x + h*(
        (44./45.)*K1 + (-56./15.)*K2 + (32./9.)*K3));
        auto K5 = f(x + h*(
        (19372./6561.)*K1 + (-25360./2187.)*K2 
        + (64448./6561.)*K3 + (-212./729.)*K4));
        auto K6  = 0/0;
        // ... 
        // xNext = x + h * (K1 * (1./6.) + K2 * (1./3.) + K3 * (1./3.) + K4 * (1./6.));
    }
    
    
    void solution(array<double, n> X0, vector<double>& ts, vector<array<double, n>>& Xs) {
        // cout << "start frm" << endl;
        // zero_count = 0;

        double t = 0;
        auto X = X0;
        
         Xs.push_back(X);
		 ts.push_back(t);
        
        double x_tau = -dx;

		DiscoQue disco(tau, disco_que_order);
        

		double prev_t0 = 0;
		bool sign_switch_flag = false;
        // double dw_step = 0;
        
        while (t < T) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;
            
            // double prev_w = w;
            // double prev_dw = dw;

			t += h;

			if (sign_switch_flag) {
				x_tau = -x_tau;
				sign_switch_flag = false;
			}

			if (t >= disco.top[0]) {
				t = disco.top[0];
				if (disco.top_order == 0) { 
                    sign_switch_flag = true; 
                }
				disco.pop();
			}
            

			step(x, dx, x_tau, t - prev_t, x, dx);

			if (x*prev_x <= 0. && x != 0.) {
				
				double x0 = x, dx0 = dx;
				double h0 = (t - prev_t);
                
                if (10* abs(x0) < abs(dx0)) { // if x0 is small and derivative is not small
                    h0 -= x0/dx0;
                    for (int newton_i = 0; newton_i < 20; newton_i++) {
                        step(prev_x, prev_dx, x_tau, h0, x0, dx0);
                        h0 -= x0/dx0;
				    }
                } else {
                    double h_l = 0;
                    double h_r = t - prev_t;
                    double h_m, x_m, dx_m;
                    for (int iii=0; iii < 50; iii++) {
                        h_m = (h_l+h_r)*0.5;
                        step(prev_x, prev_dx, x_tau, h_m, x_m, dx_m);
                        if (x_m*prev_x <= 0) h_r = h_m;
                        else                 h_l = h_m;
                    }
                    h0 = h_m;
                    x0 = x_m;
                    dx0 = dx_m;
                }
                
                double t0 = prev_t + h0;
                // cout << t0 << endl;
				disco.push({t0});
				prev_t0 = t0;
			}
            x_.push_back(x);
            dx_.push_back(dx);
            t_.push_back(t);
		}        
    };
    
};

