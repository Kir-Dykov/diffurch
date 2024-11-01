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


	DDE1 dde;
		dde.b = b; dde.c = c; 
		dde.d = d; dde.tau = tau;
		dde.h = h; dde.T = T;

    dde.solution(v, t_, x_, dx_);


	vector<vector<double>> output {x_, dx_, t_};
    
    string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(output, "output_bin/" + filename + ".bin");
    
	// save(output, "solution.bin");
}
