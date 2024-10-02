#pragma once

using namespace std;

template <typename T>
inline T find_root(const function<T(T)>& f, T l, T r, T level = 0.){
    bool as_l = f(l) > level;
    T m;
    for (int i = 0; i < 50; i++) {
        m = (r+l)*0.5;
        if ((f(m) > level) == as_l) 
        { l = m; }  else  { r = m; }
    }
    return m;
};



