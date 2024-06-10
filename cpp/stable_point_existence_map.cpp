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

using namespace std;

using json = nlohmann::json;

double v_start,v_finish,v_n,t_finish,h;
vector<double> v_;

int number_of_stable_points(double b, double c, double d, double tau, double h) {
	DDE dde;
	dde.b = b;
	dde.c = c;
	dde.d = d;
	dde.tau = tau;
	dde.h = h;

	vector<double> p_(v_n);

	for (int i = 0; i < v_n; i++) {
		double x = 0.0000000001;
		double dx = v_[i];
		double t = 0;

		vector<double> x_;
		vector<double> dx_;
		vector<double> t_;

		 x_.push_back(-dx*h);
		dx_.push_back(dx);
		 t_.push_back(-h);

		 x_.push_back(x);
		dx_.push_back(dx);
		 t_.push_back(t);

		const int disco_order = 0;
		vector<int> i_tau(disco_order+1, 0);
		vector<int> prev_i_tau(disco_order+1, 0);
		double t_last_switch = t;

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

			dde.DOPRI5step(x, dx, x_tau, t - prev_t, x, dx);


			if (sign(x) != sign(prev_x)) {
				double dx0 = dx;
				double t_switch = t;

				if (t_switch - t_last_switch > tau) {
					p_[i] = abs(dx0);
					break;
				}
				t_last_switch = t_switch;
			}

			x_.push_back(x);
			dx_.push_back(dx);
			t_.push_back(t);
		} 

		if (t >= t_finish) {
			p_[i] = 0;
		} 

	}
	
	int count_stable_points = 0;
	for (int i = 0; i < v_n - 1; i++) {
		if (
				p_[i] != 0 && p_[i+1] != 0 && // existence
				(p_[i] - v_[i])*(p_[i+1] - v_[i+1]) <= 0 && // fixed point
				abs(p_[i+1] - p_[i]) < abs(v_[i+1] - v_[i]) // stability
		   )
			count_stable_points++;
	}

	return count_stable_points;
}

int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	string params = argv[1];
	string output_prefix = argv[2];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;

	double d = (double)json_params["d"];
	double tau=(double)json_params["tau"];
	v_start = (double)json_params["v_start"];
	v_finish = (double)json_params["v_finish"];
	v_n = (int)json_params["v_n"];
	t_finish = (double)json_params["t_finish"];
	v_ = linspace(v_start, v_finish, v_n);

	double h = (double)json_params["h"];

	string mode = (string)json_params["mode"];


	int size = (int)json_params["size"];

	vector<vector<int>> N (size, vector<int>(size, 0));

	vector<thread> threads;

	if (mode == "real") {
		double l1_start  = (double)json_params["l1_start"];
		double l1_finish = (double)json_params["l1_finish"];
		double l2_start  = (double)json_params["l2_start"];
		double l2_finish = (double)json_params["l2_finish"];
		vector<double> l1 = linspace(l1_start, l1_finish, size);
		vector<double> l2 = linspace(l2_start, l2_finish, size);
		for (int i = 0; i < size; i++) {
			threads.push_back(thread( [&, i] {
				cout << i << "/" << size << endl;
				for (int j = 0; j < size; j++) {
					double b = -(l1[i]+l2[j]);
					double c = l1[i]*l2[j];
					N[i][j] = number_of_stable_points(b,c,d,tau,h);
				}
			}));
		}
	} else if (mode == "complex") {
		double re_start  = (double)json_params["re_start"];
		double re_finish = (double)json_params["re_finish"];
		double im_start  = (double)json_params["im_start"];
		double im_finish = (double)json_params["im_finish"];
		vector<double> re = linspace(re_start, re_finish, size);
		vector<double> im = linspace(im_start, im_finish, size);
		for (int i = 0; i < size; i++) {
			// cout << i << "/" << size << endl;
			for (int j = 0; j < size; j++) {
				double b = -(2*re[i]);
				double c = re[i]*re[i] + im[j]*im[j];
					N[i][j] = number_of_stable_points(b,c,d,tau,h);
					cout << i << "/" << size << endl;
				}
		}
	}
	for (int i=0; i < threads.size(); i++) {
		threads[i].join();
	}

	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

	// vector<vector<double>> output {v_, p_, t_return, zero_count};
	
	save(N, "output_bin/" + output_prefix + " " + params + ".bin");
	// save(output, "solution.bin");
}
