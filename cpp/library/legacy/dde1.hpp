#pragma once

#include "discoque.hpp"
#include <math.h>
#include <limits>

#define sNAN std::numeric_limits<double>::signaling_NaN()


using namespace std;
// const double Nan = numeric_limits<double>::quiet_NaN();

inline double dot(double x1, double dx1, double x2, double dx2) {
    return x1*x2 + dx1*dx2;
}
inline double dot(double x1, double dx1) {
    return x1*x1 + dx1*dx1;
}

class DDE1 {
public:
	double b=sNAN;
    double c=sNAN;
    double d=sNAN;
    double tau=sNAN;
	double h=sNAN;
    double T=sNAN;
    
    int disco_que_order = 2;

	inline void eval_f(const double& x, const double& dx, const double& s, double& to_x, double& to_dx ){
		to_dx = -b * dx - c * x + d * s;
		to_x  = dx;
	}

	inline void RKstep(double x, double dx, double x_tau, double h, double & xNext, double & dxNext) {
		double s = sign(x_tau); 

		double dK1 = (-b * dx - c * x + d*s); 
		double  K1 = dx; 

		double dK2 = (-b * (dx + dK1*h*0.5) - c * (x + K1*h*0.5) + d*s);
		double  K2 = (dx +dK1*h*0.5);

		double dK3 = (-b * (dx + dK2*h*0.5) - c * (x + K2*h*0.5) + d*s);
		double  K3 = (dx +dK2*h*0.5);

		double dK4 = (-b * (dx + dK3*h) - c * (x + K3*h) + d*s);
		double  K4 = (dx + dK3*h); 

		dxNext = dx + h * (1./6.) * (dK1 + 2. * dK2 + 2. * dK3 + dK4);
		xNext  =  x + h * (1./6.) * ( K1 + 2. *  K2 + 2. *  K3 +  K4);
	}

	inline void DOPRI5step(double x, double dx, double x_tau, double h, double & xNext, double & dxNext) {
		double s = sign(x_tau); 

		double  k1, k2, k3, k4, k5, k6, k7;
		double dk1,dk2,dk3,dk4,dk5,dk6,dk7;

		eval_f(x, dx, s, k1, dk1);
		eval_f( x + h*(1./5.)* k1,
				dx + h*(1./5.)*dk1,
				s, k2, dk2);
		eval_f( x + h*((3./40.)* k1 + (9./40.)* k2),
				dx + h*((3./40.)*dk1 + (9./40.)*dk2),
				s, k3, dk3);
		eval_f( x + h*((44./45.)* k1 + (-56./15.)* k2 + (32./9.)* k3),
				dx + h*((44./45.)*dk1 + (-56./15.)*dk2 + (32./9.)*dk3),
				s, k4, dk4);

		eval_f( 
				x + h*((19372./6561.)* k1 + (-25360./2187.)* k2 + (64448./6561.)* k3 + (-212./729.)* k4),
				dx + h*((19372./6561.)*dk1 + (-25360./2187.)*dk2 + (64448./6561.)*dk3 + (-212./729.)*dk4),
				s, k5, dk5);

		eval_f( 
				x + h*((9017./3168.)* k1 + (-355./33.)* k2 + (46732./5247.)* k3 + (49./176.)* k4 + (-5103./18656)* k5),
				dx + h*((9017./3168.)*dk1 + (-355./33.)*dk2 + (46732./5247.)*dk3 + (49./176.)*dk4 + (-5103./18656)*dk5),
				s, k6, dk6);


		xNext =  x + h*((35./384.)* k1 + (500./1113.)* k3 + (125./192.)* k4 + (-2187./6784)* k5 + (11./84.)* k6);
		dxNext = dx + h*((35./384.)*dk1 + (500./1113.)*dk3 + (125./192.)*dk4 + (-2187./6784)*dk5 + (11./84.)*dk6);

		eval_f( 
				xNext,
				dxNext,
				s, k7, dk7);

	}

	void step(double x, double dx, double x_tau, double h, double & xNext, double & dxNext) {
		return DOPRI5step(x,dx,x_tau,h,xNext,dxNext);
	}
    
    void first_return_map(double v, double& p, double& dp, double&t_return, double& zero_count) {
        // cout << "start frm" << endl;
        zero_count = 0;

        double t = 0;
        
        double x = 0;
        double dx = v;
        
        double w  = 0;
        double dw = 1;

        
        
        double x_tau = -dx;

		DiscoQue disco(tau, disco_que_order);
        

		double prev_t0 = 0;
		bool sign_switch_flag = false;
        double dw_step = 0;
        
        while (t < T) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;
            
            double prev_w = w;
            double prev_dw = dw;

			t += h;

			if (sign_switch_flag) {
				x_tau = -x_tau;
                dw += dw_step;
				sign_switch_flag = false;
			}

			if (t >= disco.top[0]) {
				t = disco.top[0];
				if (disco.top_order == 0) { 
                    sign_switch_flag = true; 
                    dw_step = 2*d*disco.top[2]/abs(disco.top[1]);
                }
				disco.pop();
			}
            

			step(x, dx, x_tau, t - prev_t, x, dx);
			step(w, dw, 0,     t - prev_t, w, dw);
            
            if (abs(x) > 1000) {
                t = T+1;
                break;
            }

			if (x*prev_x <= 0. && x != 0.) {
				
				zero_count++;
				double x0 = x, dx0 = dx;
				double h0 = (t - prev_t);
                double w0, dw0;
                
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
                
                step(prev_w, prev_dw, 0, h0, w0, dw0);
                
				double t0 = prev_t + h0;
				
				if (t0 - prev_t0 > tau) {
					p = abs(dx0);
                    double ddx0 = (-b * dx0 - c * x0 + d*sign(x_tau));
                    dp = sign(dx0)*(-ddx0*w0/dx0 + dw0);
					t_return = t0;
					break;
				}

				disco.push({t0, dx0, w0});
				prev_t0 = t0;
			}
		}        

		if (t >= T) {
            zero_count = NAN;
			p = NAN;
            dp = NAN;
			t_return = NAN;
		} 
        
        // cout << "finsih frm" << endl << endl;
        
    };
    
    void solution(double v, vector<double>& t_, vector<double>& x_, vector<double>& dx_) {
        // cout << "start frm" << endl;
        // zero_count = 0;

        double t = 0;
        double x = 0;
        double dx = v;
        
         x_.push_back(x);
		dx_.push_back(dx);
		 t_.push_back(t);
        
        // double w  = 0;
        // double dw = 1;

        
        
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
    
    
    void solution_value(double v, double& x, double& dx) {
        double t = 0;
        x = 0;
        dx = v;        
        
        double x_tau = -dx;

		DiscoQue disco(tau, disco_que_order);
        
		double prev_t0 = 0;
		bool sign_switch_flag = false;
        // double dw_step = 0;
        
        while (t < T) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;
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
            
            if (t >= T) {
                t = T;
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
				disco.push({t0});
				prev_t0 = t0;
			}
		}        
        
    };
    
    
    
    void benettin(double v, vector<double> &le_t, vector<double> & le1, vector<double> & le2) {
        // cout << "start frm" << endl;

        double t = 0;
        
        double x = 0;
        double dx = v;
        
        double w1  = 0;
        double dw1 = 1;
        
        double w2  = 1;
        double dw2 = 0;
        
        double le1_numerator = 0;
        double le2_numerator = 0;
        
        
        double x_tau = -dx;

		DiscoQue disco(tau, disco_que_order);
        

		double prev_t0 = 0;
		bool sign_switch_flag = false;
        double dw1_step = 0;
        double dw2_step = 0;
        
        double t_ortho_step = tau;
        double next_t_ortho = t + t_ortho_step;
        
        while (t < T) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;
            
            double prev_w1 = w1;
            double prev_dw1 = dw1;
            double prev_w2 = w2;
            double prev_dw2 = dw2;

			t += h;

			if (sign_switch_flag) {
				x_tau = -x_tau;
                dw1 += dw1_step;
                dw2 += dw2_step;
				sign_switch_flag = false;
			}

			if (t >= disco.top[0]) {
				t = disco.top[0];
				if (disco.top_order == 0) { 
                    sign_switch_flag = true; 
                    dw1_step = 2*d*disco.top[2]/abs(disco.top[1]);
                    dw2_step = 2*d*disco.top[3]/abs(disco.top[1]);
                }
				disco.pop();
			}
            

			step(x, dx, x_tau, t - prev_t, x, dx);
			step(w1, dw1, 0,     t - prev_t, w1, dw1);
			step(w2, dw2, 0,     t - prev_t, w2, dw2);
            
            // t_.push_back(t);
            // w_.push_back(w1);
            // dw_.push_back(dw1);
            
            if (t >= next_t_ortho) {
                next_t_ortho = t + t_ortho_step;
                
                // orthogonalization process
                
                double dot12 = dot(w1, dw1, w2, dw2);
                double dot11 = dot(w1, dw1);
                w2  -= dot12/dot11 *  w1;
                dw2 -= dot12/dot11 * dw1;
                
                
                for (int i = disco.disco_i[0]; i < disco.disco[0].size(); i++) {
                    disco.disco[0][i][3] -= dot12/dot11 * disco.disco[0][i][2];
                    if (disco.top_order == 0) { 
                       disco.top[3] -= dot12/dot11 * disco.top[2];
                    }
                }
                
                double dot22 = dot(w2, dw2);
                
                double norm1 = sqrt(dot11);
                double norm2 = sqrt(dot22);
                            
                w1 /= norm1;
                dw1 /= norm1;
                w2 /= norm2;
                dw2 /= norm2;
                
                for (int i = disco.disco_i[0]; i < disco.disco[0].size(); i++) {
                    disco.disco[0][i][2] /= norm1;
                    disco.disco[0][i][3] /= norm2;
                    if (disco.top_order == 0) { 
                       disco.top[2] /= norm1;
                       disco.top[3] /= norm2;
                    }
                }
                
                le1_numerator += log(norm1);
                le2_numerator += log(norm2);
                
                le_t.push_back(t);
                le1.push_back(le1_numerator / t);
                le2.push_back(le2_numerator / t);
                
                
            }
            

			if (x*prev_x <= 0. && x != 0.) {
				
				// zero_count++;
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
                
                step(prev_w1, prev_dw1, 0, h0, prev_w1, prev_dw1);
                step(prev_w2, prev_dw2, 0, h0, prev_w2, prev_dw2);
                
				double t0 = prev_t + h0;
                
                

				disco.push({t0, dx0, prev_w1, prev_w2});
				prev_t0 = t0;
                
                
                // cout << "??????????????"<< endl;
                // // for (int order=0; order < disco.disco.size(); order++) {
                // for (int order=0; order < 1; order++) {
                //     cout << "\t order=" << order << endl;
                //     for (int number = 0; number < disco.disco[order].size(); number++) {
                //     cout << "\t\t number=" << number << endl;
                //         cout << "\t\t\t t=" << disco.disco[order][number][0] << 
                //              "\t x=" << disco.disco[order][number][1] << 
                //              "\t w1=" << disco.disco[order][number][2] << 
                //              "\t w2=" << disco.disco[order][number][3] << endl;
                //     }
                // }
			}
		}        


        
        // cout << "finsih frm" << endl << endl;
        
    };
    
    
    
    
    
    
    
    
    
    
    void benettin_value(double v, double & le1, double & le2) {
        // cout << "start frm" << endl;

        double t = 0;
        
        double x = 0;
        double dx = v;
        
        double w1  = 0;
        double dw1 = 1;
        
        double w2  = 1;
        double dw2 = 0;
        
        double le1_numerator = 0;
        double le2_numerator = 0;
        
        
        double x_tau = -dx;

		DiscoQue disco(tau, disco_que_order);
        

		double prev_t0 = 0;
		bool sign_switch_flag = false;
        double dw1_step = 0;
        double dw2_step = 0;
        
        double t_ortho_step = tau;
        double next_t_ortho = t + t_ortho_step;
        
        while (t < T) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;
            
            double prev_w1 = w1;
            double prev_dw1 = dw1;
            double prev_w2 = w2;
            double prev_dw2 = dw2;

			t += h;

			if (sign_switch_flag) {
				x_tau = -x_tau;
                dw1 += dw1_step;
                dw2 += dw2_step;
				sign_switch_flag = false;
			}

			if (t >= disco.top[0]) {
				t = disco.top[0];
				if (disco.top_order == 0) { 
                    sign_switch_flag = true; 
                    dw1_step = 2*d*disco.top[2]/abs(disco.top[1]);
                    dw2_step = 2*d*disco.top[3]/abs(disco.top[1]);
                }
				disco.pop();
			}
            

			step(x, dx, x_tau, t - prev_t, x, dx);
			step(w1, dw1, 0,     t - prev_t, w1, dw1);
			step(w2, dw2, 0,     t - prev_t, w2, dw2);
            
            // t_.push_back(t);
            // w_.push_back(w1);
            // dw_.push_back(dw1);
            
            if (t >= next_t_ortho) {
                next_t_ortho = t + t_ortho_step;
                
                // orthogonalization process
                
                double dot12 = dot(w1, dw1, w2, dw2);
                double dot11 = dot(w1, dw1);
                w2  -= dot12/dot11 *  w1;
                dw2 -= dot12/dot11 * dw1;
                
                
                for (int i = disco.disco_i[0]; i < disco.disco[0].size(); i++) {
                    disco.disco[0][i][3] -= dot12/dot11 * disco.disco[0][i][2];
                    if (disco.top_order == 0) { 
                       disco.top[3] -= dot12/dot11 * disco.top[2];
                    }
                }
                
                double dot22 = dot(w2, dw2);
                
                double norm1 = sqrt(dot11);
                double norm2 = sqrt(dot22);
                            
                w1 /= norm1;
                dw1 /= norm1;
                w2 /= norm2;
                dw2 /= norm2;
                
                for (int i = disco.disco_i[0]; i < disco.disco[0].size(); i++) {
                    disco.disco[0][i][2] /= norm1;
                    disco.disco[0][i][3] /= norm2;
                    if (disco.top_order == 0) { 
                       disco.top[2] /= norm1;
                       disco.top[3] /= norm2;
                    }
                }
                
                le1_numerator += log(norm1);
                le2_numerator += log(norm2);
                
                // le_t.push_back(t);
                // le1.push_back(le1_numerator / t);
                // le2.push_back(le2_numerator / t);
                
                
            }
            

			if (x*prev_x <= 0. && x != 0.) {
				
				// zero_count++;
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
                
                step(prev_w1, prev_dw1, 0, h0, prev_w1, prev_dw1);
                step(prev_w2, prev_dw2, 0, h0, prev_w2, prev_dw2);
                
				double t0 = prev_t + h0;
                
                

				disco.push({t0, dx0, prev_w1, prev_w2});
				prev_t0 = t0;
                
                
			}
		}        
        le1= le1_numerator / t;
        le2= le2_numerator / t;
        
    };
    
    
    
    
};
