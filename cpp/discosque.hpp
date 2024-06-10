#pragma once

class DiscosQue {
public:
	double tau = 1;
	int order_limit = 5;
	vector<vector<double>> disco;
	vector<int> disco_i;

	DiscosQue (double tau_, int order_) {
		tau = tau_; 
		order_limit = order_;
		disco = vector<vector<double>>(order_limit);
		disco_i = vector<int>(order_limit, 0);
	};


	double top   = numeric_limits<double>::max();
	int    top_order = -1;

	double next(double t) { // solve t* from t = t* - tau(t*) 
		return t + tau;
	}

	void update_top() {
		top = numeric_limits<double>::max();
		top_order = -1;
		for (int order_i = 0; order_i < order_limit; order_i++) {
			if (disco_i[order_i] < disco[order_i].size() &&
					disco[order_i][disco_i[order_i]] < top) {
				top = disco[order_i][disco_i[order_i]];
				top_order = order_i;
			}
		}
	}

	void push(double t) {
		for (int order_i = 0; order_i < order_limit; order_i++) {
			t = next(t);
			disco[order_i].push_back(t);
		}
		if (top_order == -1)
			update_top();
	}

	void pop() {
		disco_i[top_order]++;
		update_top();
	}
};
