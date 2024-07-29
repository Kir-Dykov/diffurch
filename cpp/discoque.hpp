#pragma once

#include <array>

using namespace std;


class DiscoQue {
public:
	double tau = 1;
	int order_limit = 5;
	vector<vector<vector<double>>> disco;
	vector<int> disco_i;

	DiscoQue (double tau_, int order_) {
		tau = tau_; 
		order_limit = order_;
		disco = vector<vector<vector<double>>>(order_limit);
		disco_i = vector<int>(order_limit, 0);
	};


	vector<double> top {numeric_limits<double>::max()};
	int    top_order = -1;

	double next(double t) { // solve t* from t = t* - tau(t*) 
		return t + tau;
	}

	void update_top() {
        // cout << "HMM1" << endl;
        
		top[0] = numeric_limits<double>::max();
		top_order = -1;
		for (int order_i = 0; order_i < order_limit; order_i++) {
			if (disco_i[order_i] < disco[order_i].size() &&
					disco[order_i][disco_i[order_i]][0] < top[0]) {
				top = disco[order_i][disco_i[order_i]];
				top_order = order_i;
			}
		}
        // cout << "HMM1-" << endl;
        
	}

	void push(vector<double> values) {
        // cout << "HMM2" << endl;
        
		for (int order_i = 0; order_i < order_limit; order_i++) {
			values[0] = next(values[0]);
			disco[order_i].push_back(values);
		}
		if (top_order == -1)
			update_top();
        
        // cout << "HMM2-" << endl;
        
	}

	void pop() {
        // cout << "HMM3" << endl;
        
		disco_i[top_order]++;
		update_top();
        
        // cout << "HMM3-" << endl;
        
	}
};
