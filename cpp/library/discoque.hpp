#pragma once

#include <array>
#include "delay.hpp"

using namespace std;

/**
 * @brief Class that implements the queue of discontinuities.
 @tparam delays_n Number of delays
 @tparam order_limit How many propagated points to add to the queue. Setting order_limit to -1 makes the number of propagated points virtually infinite, because unsigned type is used.
 */
template <size_t delays_n, size_t order_limit>
class DiscoQue {
private:
    array<Delay, delays_n> delays;
    
    vector<pair<Real, size_t>> queue;
    
    size_t top_i = 0;
    
    Real min_ = numeric_limits<Real>::max();
    
    
public:

    DiscoQue(array<Delay, delays_n> delays) : delays(delays) {};

    const Real& min = min_;

    void push(Real t, size_t order = order_limit) {
        if constexpr (delays_n == 0 || order_limit == 0) {
            return;
        }
        if (order == 0) {
            if (top_i == queue.size()) {
                min_ = numeric_limits<Real>::max();
            } 
        } else  {
            
            for (int delay_i = 0; delay_i < delays_n; delay_i++) {
                auto propagated_values = delays[delay_i].propagate(t);
                for (auto propagated_value : propagated_values) {
                    queue.push_back(make_pair(propagated_value, order));
                }
            }
            sort(queue.begin() + top_i, queue.end());
            min_ = queue[top_i].first;
        }    
    }

    void pop() {
        while ((top_i+1)<queue.size() && queue[top_i].first == queue[top_i+1].first) {
            top_i++;
        }
        
        Real new_t = queue[top_i].first;
        size_t new_order = queue[top_i].second - 1;
        top_i++;
        push(new_t, new_order); 
	}
    
    
    /**
     * @brief Method is used to test the functionality of this class.
     * It calls push(0.) and calls pop() until queue is empty or values exceed `t_limit`. Then, it prints the contents of queue history.
     */
    void __test(Real t_limit = 100) {
        push(0.);        
        while(top_i < queue.size() && min < t_limit){
            pop();
        }
        
        cout << "Discoque Test : " << endl;
        for (auto a : queue) {
            cout << "\t(" << a.first << ",\t" << a.second << ")" << endl;
        }
    }
};








// template <int delays_n, int order_limit>
// class DiscoQue {
// // class DiscoQue<delays_n, order_limit> {
// private:
// 	array<array<vector<Real>, order_limit>, delays_n> Q_t; 
// 	array<array<         int  , order_limit>, delays_n> Q_fisrt_i; 
//     array<Real, delays_n> delays;
    
//     int all_min_delay = -1;
//     int all_min_order = -1;
//     Real all_min = numeric_limits<Real>::max();
// public:
//     const Real& min = all_min;
    
// 	DiscoQue (array<Real, delays_n> delays) : delays(delays) {
//         for (auto& row : Q_fisrt_i) { row.fill(0); }
//     };

// 	void update_all_min() {
// 		all_min = numeric_limits<Real>::max();
// 		all_min_order = -1;
//         all_min_delay = -1;
//         for (int delay_i = 0; delay_i < delays_n; delay_i++) {
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

// 	void push(Real t) {       
//         for (int delay_i = 0; delay_i < delays_n; delay_i++) {
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

