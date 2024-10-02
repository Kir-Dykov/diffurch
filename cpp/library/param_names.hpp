#pragma once

#include "string_consteval.hpp"
#include "json_unpack.hpp" // for from_json function

// using namespace std;

template <string_consteval... parameter_names>
struct ParamNames {};

template <string_consteval s, typename P>
struct param_index {
    static constexpr int value = -1;
};

template <string_consteval target_name, string_consteval name, string_consteval... rest_names>
struct param_index<target_name, ParamNames<name, rest_names...>> {
    static constexpr int value = [](){
            if constexpr (name == target_name) {return 0;} 
            else if constexpr (param_index<target_name, ParamNames<rest_names...>>::value == -1) {return -1;}
            else {  return 1 + param_index<target_name, ParamNames<rest_names...>>::value; }
        }();
};


    
template <string_consteval parameter_name, typename T>
auto& get_param(T& obj) { 
    using param_names   = T::param_names;
    return std::get<param_index<parameter_name, param_names>::value>(obj.param_tuple); 
}


template <typename T>
T from_json(json json_params) {
    using param_names   = T::param_names;
    using param_tuple_t = tuple_decay<decltype(T::param_tuple)>;
    return T(json_unpack<param_tuple_t, param_names>(json_params));
}


