#pragma once

#include <array>

using namespace std;

// also implement additional container

// there is a n_tau*order_limit of queues, 
// in which elements are added strictly in ascending order
// adding is implemented
// storing all-min value after each update
// removing all-min value 

// Q[i_delay][i_order]


template <int n_tau, int order_limit, typename T = void>
class DiscoQue {
private:
	
    
    // std::enable_if_t<!std::is_void<T>::value, 
    // array<array<vector<T>,      order_limit>, n_tau>> Q_data;
    
    std::conditional_t<std::is_void<T>::value, bool, 
    array<array<vector<double>,order_limit>, n_tau>> Q_data; 
    
    
    array<array<vector<double>, order_limit>, n_tau> Q_t; 
	array<array<     int      , order_limit>, n_tau> Q_fisrt_i; 
    
    array<double, n_tau> delays;
    
    int all_min_delay = -1;
    int all_min_order = -1;
    
    double all_min = numeric_limits<double>::max();
    
    // std::enable_if_t<!std::is_void<T>::value, T> all_min_data;
    
    std::conditional_t<std::is_void<T>::value, bool, T> all_min_data;
    
public:
    const double& min = all_min;
    const std::conditional_t<std::is_void<T>::value,  bool, T>& min_data = all_min_data; 
    
    
	DiscoQue (array<double, n_tau> delays) : delays(delays) {
        for (auto& row : Q_fisrt_i) { row.fill(0); }
    };

	void update_all_min() {
		all_min = numeric_limits<double>::max();
		all_min_order = -1;
        all_min_delay = -1;
        for (int delay_i = 0; delay_i < n_tau; delay_i++) {
            for (int order_i = 0; order_i < order_limit; order_i++) {
                int I = Q_fisrt_i[delay_i][order_i];
                if ( I < Q_t[delay_i][order_i].size() &&
                        Q_t[delay_i][order_i][I] < all_min) {
                    all_min = Q_t[delay_i][order_i][I];
                    all_min_order = order_i;
                    all_min_delay = delay_i;
                    if constexpr (!std::is_void<T>::value) {
                        all_min_data = Q_data[delay_i][order_i][I];                      
                    }
                }
            }
        }      
	}
    
    template<typename U = T>
    typename std::enable_if_t<std::is_void<T>::value && std::is_same<U, T>::value>
    push(double t) {        // !!! 
        for (int delay_i = 0; delay_i < n_tau; delay_i++) {
            for (int order_i = 0; order_i < order_limit; order_i++) {
                   Q_t[delay_i][order_i].push_back(t + delays[delay_i]*(order_i+1));
            }
        }
        if (all_min_delay == -1) {
            update_all_min();
        }
	}

    template<typename U = T>
	typename std::enable_if_t<!std::is_void<T>::value && std::is_same<U, T>::value> push(double t, U data) {  
        for (int delay_i = 0; delay_i < n_tau; delay_i++) {
            for (int order_i = 0; order_i < order_limit; order_i++) {
                   Q_t[delay_i][order_i].push_back(t + delays[delay_i]*(order_i+1));
                Q_data[delay_i][order_i].push_back(data);  // !!!
            }
        }
        if (all_min_delay == -1) {
            update_all_min();
        }
	}

	void pop() {
		Q_fisrt_i[all_min_delay][all_min_order]++;
		update_all_min();
	}
};








// template <int n_tau, int order_limit>
// class DiscoQue {
// // class DiscoQue<n_tau, order_limit> {
// private:
// 	array<array<vector<double>, order_limit>, n_tau> Q_t; 
// 	array<array<         int  , order_limit>, n_tau> Q_fisrt_i; 
//     array<double, n_tau> delays;
    
//     int all_min_delay = -1;
//     int all_min_order = -1;
//     double all_min = numeric_limits<double>::max();
// public:
//     const double& min = all_min;
    
// 	DiscoQue (array<double, n_tau> delays) : delays(delays) {
//         for (auto& row : Q_fisrt_i) { row.fill(0); }
//     };

// 	void update_all_min() {
// 		all_min = numeric_limits<double>::max();
// 		all_min_order = -1;
//         all_min_delay = -1;
//         for (int delay_i = 0; delay_i < n_tau; delay_i++) {
//             for (int order_i = 0; order_i < order_limit; order_i++) {
//                 int I = Q_fisrt_i[delay_i][order_i];
//                 if ( I < Q_t[delay_i][order_i].size() &&
//                         Q_t[delay_i][order_i][I] < all_min) {
//                     all_min = Q_t[delay_i][order_i][I];
//                     all_min_order = order_i;
//                     all_min_delay = delay_i;
//                 }
//             }
//         }      
// 	}

// 	void push(double t) {       
//         for (int delay_i = 0; delay_i < n_tau; delay_i++) {
//             for (int order_i = 0; order_i < order_limit; order_i++) {
//                    Q_t[delay_i][order_i].push_back(t + delays[delay_i]*(order_i+1));
//             }
//         }
//         if (all_min_delay == -1) {
//             update_all_min();
//         }
// 	}

// 	void pop() {
//         // NO CHECK IF EMPTY
// 		Q_fisrt_i[all_min_delay][all_min_order]++;
// 		update_all_min();
// 	};
    
// };

