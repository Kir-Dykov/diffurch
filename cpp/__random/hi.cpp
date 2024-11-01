#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

#include <boost/numeric/interval.hpp>

#include "json.hpp"

using namespace boost::numeric;

using namespace std;

using json = nlohmann::json;

int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;

	string params = argv[1];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~";

	double b = json_params["b"];
	double c = json_params["c"];
	double d = json_params["d"];
	double tau=json_params["tau"];


	interval<double> x(-1,1);
	interval<double> y = 4;
}
