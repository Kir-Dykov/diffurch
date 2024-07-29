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

double b,c,d,tau,v,T,h_l,h_r,h_n;


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
	h_l = (double)json_params["h_l"];
	h_r = (double)json_params["h_r"];
	h_n = (double)json_params["h_n"];
    int order = (int)json_params["order"];



	DDE1 dde;
		dde.b = b; dde.c = c; 
		dde.d = d; dde.tau = tau;
		// dde.h = h; 
        dde.T = T;
    
    vector<double> h_ = expspace(h_l, h_r, h_n);
    
    vector<vector<double>> error(order, vector<double>(h_n));
    
    for (int o = 0; o < order; o++) {
        dde.disco_que_order = o;
        for (int i = 0; i < h_n; i++) {
            double x, dx;
            double x_2, dx_2;

            dde.h = h_[i];
            dde.solution_value(v, x, dx);

            dde.h = 2*h_[i];
            dde.solution_value(v, x_2, dx_2);

            error[o][i] = sqrt((x-x_2)*(x-x_2) + (dx-dx_2)*(dx-dx_2));
        }
    }
    
    
    


	// vector<vector<double>> output {h_, error};
    
    
    string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(h_, "output_bin/" + filename + ".h.bin");
	save(error, "output_bin/" + filename + ".error.bin");
}
