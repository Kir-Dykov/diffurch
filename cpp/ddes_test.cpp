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
    
    
    
    
    
    const int n = 2;

	b = (double)json_params["b"];
	c = (double)json_params["c"];
	d = (double)json_params["d"];
	tau=(double)json_params["tau"];
	v = (double)json_params["v"];
	T = (double)json_params["T"];
	h = (double)json_params["h"];
    
    auto B = [b,c,d](const array<double, n>& X) {return X[0];};
    auto f0 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] - d};};
    auto f1 =[b,c,d](const array<double, n>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] + d};};
    
    RelayDDE<n, 1, 0> dde(B, f0, f1, {tau});
    
    array<double, 2> X0 = {3,4};
    
    cout << X0 * 4. << endl;
    
    // dde.print();
    
//     {
//         auto B = [b,c,d](const array<double, n>& X) {return X[0];};
//         auto f0 =[b,c,d](const array<double, n*2>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] - d};};
//         auto f1 =[b,c,d](const array<double, n*2>& X) {return array<double, n>{X[1], -b*X[1] - c*X[0] + d};};
//         RelayDDE<n, 1, 1> dde(B, f0, f1, {tau});
//     // dde.print();
        
//     }
    

    
	// save(output, "solution.bin");
}
