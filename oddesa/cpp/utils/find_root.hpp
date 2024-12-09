#pragma once

// tool for discontinuity detection

using namespace std;

template <typename T>
inline T find_root(const function<T(T)>& f, T l, T r, T level = 0.) {
    bool as_l = f(l) > level;
    T m;
    for (int i = 0; i < 50; i++) {
        m = (r+l)*0.5;
        if ((f(m) > level) == as_l) 
        { l = m; }  else  { r = m; }
    }
    return m;
};



template <typename T>
inline vector<T> find_root_ts(const function<T(T)>& f, const vector<T>& ts, T level = 0.) {
    vector<T> result;
    
    T prev_t;
    T t = ts[0];
    
    T prev_val;
    T val = f(t);
    
    for (int i = 1; i < ts.size(); i++) {
        prev_t = t;
        prev_val = val;
        
        t = ts[i];
        val = f(t);
        
        if ((val - level)*(prev_val - level) <= 0 && val != level) {
            result.push_back(find_root(f, prev_t, t, level));
        }
    }
    
    return result;
};

template <typename T>
inline vector<T> find_root_ts(const function<T(T)>& f, const vector<T>& ts, const vector<T>& levels) {
    vector<T> result;
    
    T prev_t;
    T t = ts[0];
    
    T prev_val;
    T val = f(t);
    
    size_t level_n = levels.size();
    
    // level_i : at each iteration val is in (levels[level_i - 1], levels[level_i])
    size_t level_i = distance(levels.begin(), upper_bound(levels.begin(), levels.end(), val)); //first level_i such that val < levels[level_i] 
    
    size_t ts_size = ts.size();
    
    for (int i = 1; i < ts_size; i++) {
        prev_t = t;
        prev_val = val;
        
        t = ts[i];
        val = f(t);
        
        if (level_i != 0 && val < levels[level_i - 1]) {
            result.push_back(find_root(f, prev_t, t, levels[level_i-1]));
            level_i--;
        } else if (level_i != level_n && val > levels[level_i]) {
            result.push_back(find_root(f, prev_t, t, levels[level_i]));
            level_i++;
        }
    }
    
    return result;
};


