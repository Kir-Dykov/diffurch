#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>
#include <thread>

#include "json.hpp"
#include <boost/numeric/interval.hpp>
using namespace boost::numeric;

#include "utils.hpp"
#include "save.hpp"
#include "discosque.hpp"
#include "dde1.hpp"

using namespace std;

using json = nlohmann::json;

double b_global,c_global,d_global,tau_global;
double t_finish,h;

double v_start,v_finish,v_n;
vector<double> v_;

void set_as_l1l2(double x, double y, double& b, double& c, double& d, double& tau) {
	b = -(x+y);
	c = x*y;
	d = d_global;
	tau = tau_global;
}


void set_as_reim(double x, double y, double& b, double& c, double& d, double& tau) {
	b = -2*x;
	c = x*x+y*y;
	d = d_global;
	tau = tau_global;
}

void set_as_taub(double x, double y, double& b, double& c, double& d, double& tau) {
	c = c_global;
	d = d_global;
	tau = x;
	b = y;
}

void (*set_params)(double x, double y, double& b, double& c, double& d, double& tau);

void count_fixed_points(double x, double y, int& count_stable, int& count_unstable) {
	DDE1 dde;

	dde.h = h; 
	dde.t_finish = t_finish;

	set_params(x, y, dde.b, dde.c, dde.d, dde.tau);


	vector<double> p_(v_n);
	vector<double> t_return(v_n);
	vector<double> zero_count(v_n, 0);

	for (int i = 0; i < v_n; i++) {
		double x = 0;
		double dx = v_[i];
		double t = 0;

		double x_tau = -dx;

		DiscosQue discos(dde.tau, 6);

		double prev_t0 = 0;
		bool sign_switch_flag = false;

		while (t < t_finish) {
			double prev_x = x;
			double prev_dx = dx;
			double prev_t = t;

			t += h;


			if (sign_switch_flag) {
				x_tau = -x_tau;
				sign_switch_flag = false;
			}


			if (t >= discos.top) {
				t = discos.top;
				if (discos.top_order == 0) { sign_switch_flag = true; }
				discos.pop();
			}

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
				
				if (t0 - prev_t0 > dde.tau) {
					p_[i] = abs(dx0);
					t_return[i] = t0;
					break;
				}

				discos.push(t0);
				prev_t0 = t0;
			}
		} 

		if (t >= t_finish) {
			p_[i] = 0;
			t_return[i] = 0;
		} 
	}

	for (int i = 1; i < v_n; i++) {
		if ((p_[i-1]-v_[i-1])*(p_[i]-v_[i]) <= 0) {
				count_stable++;
			// if (abs((p_[i]-p_[i-1])/(v_[i]-v_[i-1])) < 1) {
			// 	count_stable++;
			// } else {
			// 	count_unstable++;
			// }
		}

	}
}

int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	string params = argv[1];
	string output_prefix = argv[2];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;

	// d = (double)json_params["d"];
	// tau=(double)json_params["tau"];
	v_start = (double)json_params["v_start"];
	v_finish = (double)json_params["v_finish"];
	v_n = (int)json_params["v_n"];
	v_ = linspace(v_start, v_finish, v_n);

	if (json_params.contains("b"))   b_global   = (double)json_params["b"];
	if (json_params.contains("c"))   c_global   = (double)json_params["c"];
	if (json_params.contains("d"))   d_global   = (double)json_params["d"];
	if (json_params.contains("tau")) tau_global = (double)json_params["tau"];

	h = (double)json_params["h"];
	t_finish = (double)json_params["t_finish"];


	string mode = (string)json_params["mode"];
	if (mode == "l1l2") {
		set_params = set_as_l1l2;
	} else if (mode == "reim") {
		set_params = set_as_reim;
	} else if (mode == "taub") {
		set_params = set_as_taub;
	}


	int size = (int)json_params["size"];

	vector<vector<int>> count_stable (size, vector<int>(size, 0));
	vector<vector<int>> count_unstable (size, vector<int>(size, 0));

	double x_start = (double)json_params["x_start"];
	double x_finish = (double)json_params["x_finish"];
	vector<double> x = linspace(x_start, x_finish, size);
	double y_start = (double)json_params["y_start"];
	double y_finish = (double)json_params["y_finish"];
	vector<double> y = linspace(y_start, y_finish, size);

	vector<thread> threads;

	for (int i = 0; i < size; i++) {
		threads.push_back(thread( [&, i] {
			for (int j = 0; j < size; j++) {
				count_fixed_points(x[i],y[j], count_stable[i][j], count_unstable[i][j]);
			}
			cout << i << "th / " << size << " is done" << endl; 
		}));
	}

	for (int i=0; i < threads.size(); i++) {
		threads[i].join();
	}

	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

	// vector<vector<double>> output {v_, p_, t_return, zero_count};
	
	save(count_stable, "output_bin/" + output_prefix + " " + params + ".bin");
	// save(output, "solution.bin");
}
