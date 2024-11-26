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

void count_fixed_points(double x, double y, int& count_stable, int& count_unstable) {
	DDE1 dde;

	dde.h = h; 
	dde.T = T;

	set_params(x, y, dde.b, dde.c, dde.d, dde.tau);
     
    double _;
    
    double p, dp, prev_p, prev_dp;
	dde.first_return_map(v_[0], p, dp, _, _);
    
	for (int i = 1; i < v_n; i++) {
        prev_p = p;
        prev_dp = dp;
        dde.first_return_map(v_[i], p, dp, _, _);
        
        if (isnan(p) != isnan(prev_p)) {
            
            double l_v = v_[i-1], r_v = v_[i];
            double l_p = prev_p, r_p = p;
            
            bool left_is_nan =  l_p;
            
            double _;
            for (int bs_i = 0; bs_i < 50; bs_i++) {
                double m_v = (l_v+r_v)*0.5;
                double m_p;
                dde.first_return_map(m_v, m_p, _, _, _);
                if (isnan(m_p) == left_is_nan)
                     {l_v = m_v; l_p = m_p;}
                else {r_v = m_v; r_p = m_p;}
            }
            // v_[i-1] = l_v;
            // v_[i]   = r_v;
            dde.first_return_map(l_v, prev_p, prev_dp, _, _);     
            dde.first_return_map(r_v, p, dp, _, _);
        } else if (isnan(p)) continue;
        
		if (!isnan(p) && !isnan(prev_p) && (prev_p-v_[i-1])*(p-v_[i]) <= 0) {
            double v = v_[i];
            for (int nwtn_i = 0; nwtn_i < 20; nwtn_i++) {
                v -= (p - v)/(dp - 1);
                dde.first_return_map(v, p, dp, _, _);
            }
            
            if (abs(p - v) < 0.001 && abs(v - v_[i]) <= abs(v_[i] - v_[i-1]) ) {
                if (abs(dp) < 1) {
                    count_stable++;
                } else {
                    count_unstable++;
                }
            }
		}
        
        // find extremum
        if (dp * prev_dp < 0) {
            
            double l_v = v_[i-1], r_v = v_[i];
            double l_p = prev_p, r_p = p;
            
            
            double _;
            for (int bs_i = 0; bs_i < 50; bs_i++) {
                double m_v = (l_v+r_v)*0.5;
                double m_p;
                double m_dp;
                dde.first_return_map(m_v, m_p, m_dp, _, _);
                if (m_dp * prev_dp > 0)
                     {l_v = m_v; l_p = m_p;}
                else {r_v = m_v; r_p = m_p;}
            }
            // v_[i-1] = l_v;
            // v_[i]   = r_v;
            dde.first_return_map(l_v, prev_p, prev_dp, _, _);     
            dde.first_return_map(r_v, p, dp, _, _);
        }

	}
}


            
// 			double left_sign = (prev_p-v_[i-1]);
            
//             double l_v = v_[i-1];
//             double r_v = v_[i];
            
//             double l_p = prev_p;
//             double r_p = p;
            
            
//             for (int bs_i = 0; bs_i < 50; bs_i++) {
//                 double m_v = (l_v+r_v)*0.5;
//                 double m_p;
//                 dde.first_return_map(m_v, m_p, dp, _, _);
//                 if (left_sign*(m_p - m_v) > 0) {
//                     l_v = m_v;
//                     l_p = m_p;
//                 } else {
//                     r_v = m_v;
//                     r_p = m_p;
//                 }
//             }
            
//             if  (abs(l_p - l_v) > 0.001) {
//                 // we regard this as discontinuity
//             } else {
//                 // we regard this as fixed point
//                 if (abs(dp) < 1) {
//                     count_stable++;
//                 } else {
//                     count_unstable++;
//                 }
                
//             }

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
    
    

	vector<vector<int>> count_stable (size, vector<int>(size, 0));
	vector<vector<int>> count_unstable (size, vector<int>(size, 0));

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
				count_fixed_points(x[j],y[i], count_stable[i][j], count_unstable[i][j]);
                progress_bar.progress += 1;
                progress_bar.update();
                progress_bar.print();
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
