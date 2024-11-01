#pragma once

// directly integration related


#include <array>
#include "delay.hpp"

using namespace std;


template <size_t delays_n, size_t order_limit>
class DiscoQue {
public:
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