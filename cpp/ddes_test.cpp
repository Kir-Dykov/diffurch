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
// #include "dde1.hpp"
#include "ddes.hpp"

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
    
    const int n = 2;
    
    auto B = [b,c,d](const array<double, n>& X) {return X[0];};
    auto B_derivative = [b,c,d](const array<double, n>& X) {return array<double, n>{1, 0};};
    auto f0 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] - d};};
    auto f1 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] + d};};
    RelayDDE<n, 1, 0> dde(B, B_derivative, f0, f1, {tau});
    
    auto phi            = [b,c,d](double t) {return array<double, n>{0.8*d/c, 0};};
    auto phi_derivative = [b,c,d](double t) {return array<double, n>{0, 0};};
    

    
    
    auto ivp = InitialValueProblem<n, 1, 0>(dde, phi, phi_derivative);
    
    ivp.solution(h, T);
    
  
    vector<vector<double>> output = vector<vector<double>>(ivp.ts.size(), vector<double>(1 + n));
    for (int i = 0; i < ivp.ts.size(); i++) {
        output[i][0] = ivp.ts[i];
        copy(ivp.Xs[i].begin(), ivp.Xs[i].end(), output[i].begin()+1);
    }
    
    string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(output, "output_bin/" + filename + ".bin");
    
	// save(output, "solution.bin");
    
    cout << "~~~ " << __FILE__ << " is finished ~~~" << endl;
}
