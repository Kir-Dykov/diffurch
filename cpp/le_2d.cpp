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
#include "discoque.hpp"
#include "dde1.hpp"
#include "progressbar.hpp"

using namespace std;

using json = nlohmann::json;

double b_global,c_global,d_global,tau_global;
double T,h;

double v_l,v_r,v_n;
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

void le_s(double x, double y, double& le1, double& le2, double& le_p) {
	DDE1 dde;

	dde.h = h; 
	dde.T = T;

	set_params(x, y, dde.b, dde.c, dde.d, dde.tau);
     
    double _;
    
    DDE1 dde;
		dde.b = b; dde.c = c; 
		dde.d = d; dde.tau = tau;
		dde.h = h; dde.T = T;
    
    
    double t = 0;
    double k = 0;
    double v_k = v;
    double dp, T_k;
    double _;
    double sum_ln_dp = 0;
    
    while (t < T) {
        dde.first_return_map(v_k, v_k, dp, T_k, _);
        t += T_k;
        sum_ln_dp += log(abs(dp));
    }
    le_p = sum_ln_dp/t;
    
    // vector<double> t_w, w_, dw_;
    
    dde.benettin(v, lambda_B_t, lambda_B_1, lambda_B_2);
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
	v_l = (double)json_params["v_l"];
	v_r = (double)json_params["v_r"];
	v_n = (int)json_params["v_n"];
    
    if (json_params.contains("expspace") && (bool)json_params["expspace"]) {
        v_ = expspace(v_l, v_r, v_n);
    }
    else
        v_ = linspace(v_l, v_r, v_n);

	if (json_params.contains("b"))   b_global   = (double)json_params["b"];
	if (json_params.contains("c"))   c_global   = (double)json_params["c"];
	if (json_params.contains("d"))   d_global   = (double)json_params["d"];
	if (json_params.contains("tau")) tau_global = (double)json_params["tau"];

	h = (double)json_params["h"];
	T = (double)json_params["T"];


	string mode = (string)json_params["mode"];
	if (mode == "l1l2") {
		set_params = set_as_l1l2;
	} else if (mode == "reim") {
		set_params = set_as_reim;
	} else if (mode == "taub") {
		set_params = set_as_taub;
	}


	int size = (int)json_params["size"];
    
    

	vector<vector<double>> le_1 (size, vector<int>(size, 0));
	vector<vector<double>> le_2 (size, vector<int>(size, 0));
	vector<vector<double>> le_p (size, vector<int>(size, 0));

	double x_l = (double)json_params["x_l"];
	double x_r = (double)json_params["x_r"];
	vector<double> x = linspace(x_l, x_r, size);
	double y_l = (double)json_params["y_l"];
	double y_r = (double)json_params["y_r"];
	vector<double> y = linspace(y_l, y_r, size);

	vector<thread> threads;

    
    ProgressBar progress_bar (size*size);
    
	for (int i = 0; i < size; i++) {
		threads.push_back(thread( [&, i] {
			for (int j = 0; j < size; j++) {                
				le_s(x[j],y[i], le_1[i][j], le_2[i][j], le_p[i][j]);
                progress_bar.progress += 1;
                // progress_bar.update();
                // progress_bar.print();
			}
            progress_bar.update();
            progress_bar.print();
			// cout << i << "th / " << size << " is done" << endl; 
		}));
	}

	for (int i=0; i < threads.size(); i++) {
		threads[i].join();
	}

	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

	vector<vector<vector<int>>> output {count_stable, count_unstable};
	
    string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(output, "output_bin/" + filename + ".bin");
	// save(output, "solution.bin");
}
