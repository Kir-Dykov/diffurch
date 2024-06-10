#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

#include "json.hpp"
#include <boost/numeric/interval.hpp>
using namespace boost::numeric;

#include "utils.hpp"
#include "save.hpp"

using namespace std;

using json = nlohmann::json;

double b,c,d,tau,v,t_finish,h;


interval<double> sign(interval<double> val) {
	return interval<double>((double)((val.lower() > 0.) - (val.upper() < 0.)));
}

template <typename T>
inline void eval_f(const T& x, const T& dx, const T& s, T& to_x, T& to_dx){
	to_dx = -b * dx - c * x + d * s;
	to_x  = dx;
}

template <typename T>
void RK4step(T x, T dx, T x_tau, T h, T & xNext, T & dxNext) {
        T s = sign(x_tau); 

		T  k1, k2, k3, k4;
		T dk1,dk2,dk3,dk4;

		eval_f(x, dx, s, k1, dk1);
		eval_f( x + h*(0.5)* k1,
			   dx + h*(0.5)*dk1,
			   s, k2, dk2);
		eval_f( x + h*(0.5)* k2,
			   dx + h*(0.5)*dk2,
			   s, k3, dk3);
		eval_f( x + h* k3,
			   dx + h*dk3,
			   s, k4, dk4);

		// T dK1 = (-b * dx - c * x + d*s); 
  //       T  K1 = dx; 
  //       
  //       T dK2 = (-b * (dx + dK1*h*0.5) - c * (x + K1*h*0.5) + d*s);
  //       T  K2 = (dx +dK1*h*0.5);
  //       
  //       T dK3 = (-b * (dx + dK2*h*0.5) - c * (x + K2*h*0.5) + d*s);
  //       T  K3 = (dx +dK2*h*0.5);
  //       
  //       T dK4 = (-b * (dx + dK3*h) - c * (x + K3*h) + d*s);
  //       T  K4 = (dx + dK3*h); 
        
		xNext  =  x + h * (1./6.) * ( k1 + 2. *  k2 + 2. *  k3 +  k4);
        dxNext = dx + h * (1./6.) * (dk1 + 2. * dk2 + 2. * dk3 + dk4);
    }

template <typename T>
void DOPRI5step(T x, T dx, T x_tau, T h, T & xNext, T & dxNext) {
        T s = sign(x_tau); 

		T  k1, k2, k3, k4, k5, k6, k7;
		T dk1,dk2,dk3,dk4,dk5,dk6,dk7;

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

int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;

	string params = argv[1];
	string output_prefix = argv[2];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;

	b = (double)json_params["b"];
	c = (double)json_params["c"];
	d = (double)json_params["d"];
	tau=(double)json_params["tau"];
	v = (double)json_params["v"];
	t_finish = (double)json_params["t_finish"];
	h = (double)json_params["h"];

	// interval<double> x = 0;
	// interval<double>dx = v;
	// interval<double> t = 0;

	double x = 0.000000000001 * v;
	double dx = v;
	double t = 0;

	vector<double> x_;
	vector<double> dx_;
	vector<double> t_;

	// vector<interval<double>> x_;
	// vector<interval<double>> dx_;
	// vector<interval<double>> t_;

		 x_.push_back(-dx*h);
		dx_.push_back(dx);
		 t_.push_back(-h);

		 x_.push_back(x);
		dx_.push_back(dx);
		 t_.push_back(t);


	const int disco_order = 0;
	vector<int> i_tau(disco_order+1, 0);
	vector<int> prev_i_tau(disco_order+1, 0);


	while (t < t_finish) {

		double prev_x = x;
		double prev_dx = dx;
		double prev_t = t;
		double x_tau;

		t += h;

		for (int ii = 0; ii <= disco_order; ii++) {
			prev_i_tau[ii] = i_tau[ii];
			while (t_[i_tau[ii] + 1] < (t -(ii+1)*tau)) {
				i_tau[ii]++;
			}
		}
		x_tau = x_[i_tau[0]];	

		for (int ii = 0; ii < disco_order; ii++) {
			// break;
			if (x_[i_tau[ii]] * x_[prev_i_tau[ii]] <= 0) {
				// cout << ii << " opa " << endl;
				double x0 = x_[prev_i_tau[ii]];
				double dx0=dx_[prev_i_tau[ii]];
				double h0 = 0 - x0/dx0;
				double x1,dx1;

				for (int newton_i = 0; newton_i < 20; newton_i++) {
					DOPRI5step(x0, dx0, x_[prev_i_tau[ii+1]], h0, x1, dx1);
					h0 -= x1/dx1;
				}
				t = t_[prev_i_tau[ii]] + h0 + (ii+1) * tau;
				x_tau = x_[prev_i_tau[ii]];	
				cout << t_[prev_i_tau[ii]] << "	" << h0 << "   " << x1 << "  " << dx1 << " " << t << endl;
				break;
			}
		}
		DOPRI5step(x, dx, x_tau, t - prev_t, x, dx);

		x_.push_back(x);
		dx_.push_back(dx);
		t_.push_back(t);
	}

	vector<vector<double>> output {x_, dx_, t_};

	save(output, "output_bin/" + output_prefix + " " + params + ".bin");
	// save(output, "solution.bin");
}
