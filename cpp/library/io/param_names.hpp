#pragma once

// tool for dde definition

#include "string_consteval.hpp"
#include "json_unpack.hpp" // for from_json function
#include "macros.hpp"



// using namespace std;

#define PARAMS(ClassName, ...) \
    Real __VA_ARGS__; \
    using param_names = ParamNames<TO_STRINGS(__VA_ARGS__)>; \
    AUTO(params, tie(__VA_ARGS__)); \
    ClassName (const tuple_add_const<decltype(params)>& params_) : ClassName() { copy_tuple(params_, params); } \
    ClassName (PREPEND_(double, APPEND(_, __VA_ARGS__)))         : ClassName() { INIT_STATEMENTS(__VA_ARGS__) }




// copy_tuple(params_, params);   

// Empty struct, that holds parameter names in its template arguments
template <string_consteval... parameter_names>
struct ParamNames {};

// base template
template <string_consteval s, typename P>
struct param_index { static constexpr int value = -1; };

// template specialization
template <string_consteval target_name, string_consteval name, string_consteval... rest_names>
struct param_index<target_name, ParamNames<name, rest_names...>> {
    static constexpr int value = [](){
            if constexpr (name == target_name) {return 0;} 
            else if constexpr (param_index<target_name, ParamNames<rest_names...>>::value == -1) {return -1;}
            else {  return 1 + param_index<target_name, ParamNames<rest_names...>>::value; }
        }();
};


template <typename T>
concept has_param_names = requires (T t) {
    typename T::param_names;
};

template <string_consteval parameter_name, typename T>
requires has_param_names<T>
auto& get_param(T& obj) { 
    using param_names   = T::param_names;
    return std::get<param_index<parameter_name, param_names>::value>(obj.param_tuple); 
}

// template <typename T, typename Arg>
// concept has_params_json_constructor = requires (T t) {
//     T::param_names;
// };

template <typename T>
requires has_param_names<T>
T from_json(json json_params) {
    using param_names   = T::param_names;
    using params_t = tuple_decay<decltype(T::params)>;
    return T(json_unpack<params_t, param_names>(json_params));
}


