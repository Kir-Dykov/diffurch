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
#include "discoque.hpp"
#include "dde1.hpp"

using namespace std;

using json = nlohmann::json;

double b,c,d,tau,v,T,h;


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
	T = (double)json_params["T"];
	h = (double)json_params["h"];


	vector<double> x_;
	vector<double> dx_;
	vector<double> t_;
    
    vector<double> v_;
    vector<double> lambda_p;
    vector<double> lambda_T;
    vector<double> lambda_T_t;
    vector<double> lambda_B_1;
    vector<double> lambda_B_2;
    vector<double> lambda_B_t;
    
    


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
        k += 1.;
        sum_ln_dp += log(abs(dp));
        
        v_.push_back(v_k);
        lambda_p.push_back(sum_ln_dp/k);
        lambda_T.push_back(sum_ln_dp/t);
        lambda_T_t.push_back(t);
    }
    
    // vector<double> t_w, w_, dw_;
    
    dde.benettin(v, lambda_B_t, lambda_B_1, lambda_B_2);
    

    // dde.solution(v, t_, x_, dx_);


	vector<vector<double>> output {v_, lambda_p, lambda_T, lambda_T_t};
    
    vector<vector<double>> output2 {lambda_B_t, lambda_B_1, lambda_B_2};
    
    // vector<vector<double>> output3 {t_w, w_, dw_};
    
    
    string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(output, "output_bin/" + filename + ".p.bin");
	save(output2, "output_bin/" + filename + ".benettin.bin");
	// save(output3, "output_bin/" + filename + ".w.bin");
    
	// save(output, "solution.bin");
}
