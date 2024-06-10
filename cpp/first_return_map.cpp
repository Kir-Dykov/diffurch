#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>
#include <limits>

#include "json.hpp"
#include <boost/numeric/interval.hpp>
using namespace boost::numeric;

#include "utils.hpp"
#include "save.hpp"
#include "discosque.hpp"
#include "dde1.hpp"

using namespace std;

using json = nlohmann::json;

double b,c,d,tau,v_start,v_finish,t_finish,h;

int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	string params = argv[1];
	string output_prefix = argv[2];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;

	b = (double)json_params["b"];
	c = (double)json_params["c"];
	d = (double)json_params["d"];
	tau=(double)json_params["tau"];
	double v_start = (double)json_params["v_start"];
	double v_finish = (double)json_params["v_finish"];
	int v_n = (int)json_params["v_n"];
	t_finish = (double)json_params["t_finish"];
	h = (double)json_params["h"];

	vector<double> v_ = linspace(v_start, v_finish, v_n);
	vector<double> p_(v_n);
	vector<double> t_return(v_n);
	vector<double> zero_count(v_n, 0);
	for (int i = 0; i < v_n; i++) {

		DDE1 dde;
		dde.b = b; dde.c = c; 
		dde.d = d; dde.tau = tau;
		dde.h = h; dde.t_finish = t_finish;


		double x = 0;
		double dx = v_[i];
		double t = 0;

		double x_tau = -dx;

		// vector<double> x_;
		// vector<double> dx_;
		// vector<double> t_;

		// // initial history
		//  x_.push_back(-dx*h);
		// dx_.push_back(dx);
		//  t_.push_back(-h);
		//
		//  // initial value
		//  x_.push_back(x);
		// dx_.push_back(dx);
		//  t_.push_back(t);

		DiscosQue discos(tau, 6);

		// int t_retarded_i = 0;
		double prev_t0 = 0;
		bool sign_switch_flag = false;

		while (t < t_finish) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;
			// double x_tau;

			t += h;

			if (sign_switch_flag) {
				x_tau = -x_tau;
				sign_switch_flag = false;
			}

			// while (t_[t_retarded_i] < t - tau) t_retarded_i++;

			if (t >= discos.top) {
					// cout << discos.top << endl;
				t = discos.top;
				if (discos.top_order == 0) { sign_switch_flag = true; }
				discos.pop();
				// while (t_[t_retarded_i] >= t-tau) t_retarded_i--;
			}

			// x_tau = x_[t_retarded_i];

			dde.step(x, dx, x_tau, t - prev_t, x, dx);

			if (x*prev_x <= 0. && x != 0.) {
				
				zero_count[i]++;
				double x0 = x, dx0 = dx;
				double h0 = (t - prev_t) - x0/dx0;
				for (int newton_i = 0; newton_i < 20; newton_i++) {
					dde.step(prev_x, prev_dx, x_tau, h0, x0, dx0);
					h0 -= x0/dx0;
				}
				double t0 = prev_t + h0;
				
				// cout << "zero at " <<	t0 << endl;

				if (t0 - prev_t0 > tau) {
					// cout << t0 << "    "<< h0 << "   " << x0 << "  " << dx0 << endl;
					p_[i] = abs(dx0);
					t_return[i] = t0;
					break;
				}

				discos.push(t0);
				// cout << discos.top << " " << discos.top_order << endl;
				prev_t0 = t0;
			}

			// x_.push_back(x);
			// dx_.push_back(dx);
			// t_.push_back(t);
		} 

		if (t >= t_finish) {
			p_[i] = 0;
			t_return[i] = 0;
		} 
	}

	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

	vector<vector<double>> output {v_, p_, t_return, zero_count};
	
	save(output, "output_bin/" + output_prefix + " " + params + ".bin");
	// save(output, "solution.bin");
}
