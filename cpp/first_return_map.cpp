#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>
#include <limits>

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

double b,c,d,tau,v_l,v_r,T,h;

int main(int argc, char* argv[]) {
    cout << "~~~ " << __FILE__ << " is executed ~~~" << endl;
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	string params = argv[1];
	string output_prefix = argv[2];
	json json_params = json::parse(params);

	cout << "~~~  parameters: " << params << " ~~~" << endl;

	b = (double)json_params["b"];
	c = (double)json_params["c"];
	d = (double)json_params["d"];
	tau=(double)json_params["tau"];
	double v_l = (double)json_params["v_l"];
	double v_r = (double)json_params["v_r"];
	int v_n = (int)json_params["v_n"];
	T = (double)json_params["T"];
	h = (double)json_params["h"];


	vector<double> v_;
    
    if (json_params.contains("expspace") && (bool)json_params["expspace"])
        v_ = expspace(v_l, v_r, v_n);
    else
        v_ = linspace(v_l, v_r, v_n);
    
    DDE1 dde;
    dde.b = b; dde.c = c; 
    dde.d = d; dde.tau = tau;
    dde.h = h; dde.T = T;
    
    ProgressBar progress_bar (v_n);

    
	vector<double> p_(v_n);
    vector<double> dp_(v_n);
	vector<double> t_return(v_n);
	vector<double> zero_count(v_n, 0);
    // base graph
	for (int i = 0; i < v_n; i++) {
        dde.first_return_map(v_[i], p_[i], dp_[i], t_return[i], zero_count[i]);      
        
        progress_bar.progress += 1;
        if ((progress_bar.progress & 63) == 0) {
            progress_bar.update();
            progress_bar.print();
        }
        
	}
    progress_bar.update();
    progress_bar.print();    
    
    vector<double> stable_fixed_points;
    vector<double> unstable_fixed_points;
    vector<double> all_fixed_points;
    vector<int> disco_index;
    // vector<double> discontinuities;
    
    // discontinuities graph correction and fixed point detection
    for (int i = 1; i < v_n; i++) {
        // cout << "!" << endl;
        if (isnan(p_[i]) != isnan(p_[i-1])) {
            bool left_is_nan = isnan(p_[i-1]);
            
            double l_v = v_[i-1], r_v = v_[i];
            double l_p = p_[i-1], r_p = p_[i];
            double _;
            for (int bs_i = 0; bs_i < 50; bs_i++) {
                double m_v = (l_v+r_v)*0.5;
                double m_p;
                dde.first_return_map(m_v, m_p, _, _, _);
                if (isnan(m_p) == left_is_nan)
                     {l_v = m_v; l_p = m_p;}
                else {r_v = m_v; r_p = m_p;}
            }
            v_[i-1] = l_v;
            v_[i]   = r_v;
            dde.first_return_map(v_[i-1], p_[i-1], dp_[i-1], t_return[i-1], zero_count[i-1]);     
            dde.first_return_map(v_[i], p_[i], dp_[i], t_return[i], zero_count[i]);
        } else if (isnan(p_[i])) continue;
        
        
        double dp_finite_difference = (p_[i]-p_[i-1])/(v_[i] - v_[i-1]);
        double dp_mean            = 0.5*(dp_[i]+dp_[i-1]);
        
        if (zero_count[i]!=zero_count[i-1] || (!isnan(dp_mean) && abs(dp_finite_difference - dp_mean)/abs(dp_mean) > 1)) { // disco
            // cout << "disco at " << v_[i] << "!" << endl;
            // cout << "\tdp_fd = " << dp_finite_difference << endl;
            // cout << "\tdp_mean = " << dp_mean << endl;
            
            double l_v = v_[i-1], r_v = v_[i];
            double l_p = p_[i-1], r_p = p_[i];
            
            // cout << "\tl_p = " << l_p << endl;
            // cout << "\tr_p = " << r_p << endl;
            double _;
            for (int bs_i = 0; bs_i < 50; bs_i++) {
                double m_v = (l_v+r_v)*0.5;
                double m_p;
                dde.first_return_map(m_v, m_p, _, _, _);
                if (abs(m_p - p_[i-1]) < abs(m_p - p_[i]))
                     {l_v = m_v; l_p = m_p;}
                else {r_v = m_v; r_p = m_p;}
            }
            v_[i-1] = l_v;
            v_[i]   = r_v;
            dde.first_return_map(v_[i-1], p_[i-1], dp_[i-1], t_return[i-1], zero_count[i-1]);     
            dde.first_return_map(v_[i], p_[i], dp_[i], t_return[i], zero_count[i]);
            
            if (p_[i-1] < 0.001 && dp_[i-1] < -10) {
                p_[i-1] = 0;
            } else
            if (p_[i] < 0.001 && dp_[i] > 10) {
                p_[i] = 0;
            }
            
            disco_index.push_back(i);
            
//             cout << "\tafter iterations:" << endl;
//             cout << "\tl_p = " << l_p << endl;
//             cout << "\tr_p = " << r_p << endl;
//             cout << "\twidth = " << r_v - l_v << endl;
        } else if (!isnan(p_[i]) && !isnan(p_[i-1]) && (p_[i-1]-v_[i-1])*(p_[i]-v_[i]) <= 0) { // disco or fixed point
			double left_sign = (p_[i-1]-v_[i-1]);
            double l_v = v_[i-1], r_v = v_[i];
            double l_p = p_[i-1], r_p = p_[i];
            double m_dp;
            double _;
            
            for (int bs_i = 0; bs_i < 50; bs_i++) {
                double m_v = (l_v+r_v)*0.5, m_p;
                dde.first_return_map(m_v, m_p, m_dp, _, _);
                if (left_sign*(m_p - m_v) > 0) 
                     {l_v = m_v; l_p = m_p;}
                else {r_v = m_v; r_p = m_p;}
            }
            
            if  (abs(l_p - l_v) <= 0.001) {
                if (abs(m_dp) < 1)
                       stable_fixed_points.push_back(l_v);
                else unstable_fixed_points.push_back(l_v);  
                
                all_fixed_points.push_back(l_v);
            }
		}
	}
    
    
    vector<double> invariant_interval_l;
    vector<double> invariant_interval_r;
    
    for (int i = 0; i < all_fixed_points.size(); i++) {
        int j = 1;
        while (v_[j] < all_fixed_points[i]) j++;        
        int j_l = j-1;
        int j_r = j;
        
        int prev_j_l = j_l;
        int prev_j_r = j_r;
        
        double max_p = max(p_[j_l], p_[j_r]);
        double min_p = min(p_[j_l], p_[j_r]);
        
        while (j_r < v_n - 1 && 
               (max_p > v_[j_r] || (j_l > 0 && min_p < v_[j_l]))) {
            // what if min_p < v_[0]???
          while (j_r < v_n - 1 && v_[j_r] < max_p) j_r++;
          while (j_l > 0       && v_[j_l] > min_p) j_l--;
            
          for (int k = prev_j_r + 1; k <= j_r; k++) {
              max_p = max(max_p, p_[k]);
              min_p = min(min_p, p_[k]);
          }
          for (int k = prev_j_l - 1; k >= j_l; k--) {
              max_p = max(max_p, p_[k]);
              min_p = min(min_p, p_[k]);
          }
            
          prev_j_l = j_l;
          prev_j_r = j_r;
        }
        
        
        if (j_r != v_n - 1) {
            double l, r;
            
            if (j_l == 0) 
                 l = 0.;
            else l = v_[j_l];
            r = v_[j_r];
            
            invariant_interval_l.push_back(l);
            invariant_interval_r.push_back(r);
            
            // for (int ii = 0; ii < invariant_interval_l.size(); ii++) {
            //     if (l <= invariant_interval_l[ii] && r >= invariant_interval_r[ii]) {
            //         invariant_interval_l = l;
            //         invariant_interval_l = r;
            //     }
            // }
        }
        
        
    }
    

	chrono::steady_clock::time_point end = chrono::steady_clock::now();
	int seconds = chrono::duration_cast<chrono::seconds>(end - begin).count();
	cout << "~~~ Computation took " << (seconds / 3600) << ":" << (seconds / 60) % 60 << ":" << seconds % 60 << " (hh:mm:ss) ~~~" << endl;

	vector<vector<double>> output {v_, p_, dp_, t_return, zero_count};
    vector<vector<double>> output2 {invariant_interval_l, invariant_interval_r};
	
	string filename = output_prefix + " " + params;
    if (filename.size() > 200)
        filename.erase(200, string::npos);

	save(output, "output_bin/" + filename + ".bin");
	save(stable_fixed_points, "output_bin/" + filename + ".stable_fixed_points.bin");
	save(unstable_fixed_points, "output_bin/" + filename + ".unstable_fixed_points.bin");
	save(disco_index, "output_bin/" + filename + ".disco_index.bin");
    save(output2, "output_bin/" + filename + ".invariant_interval.bin");
    	// save(output, "solution.bin");
}
