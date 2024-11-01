#pragma once

// tool for params

#include "json.hpp"
#include "macros.hpp"
#include "string_consteval.hpp"
#include "tuple.hpp"


// auto [t_finish, h] = json_unpack<Real, "t_finish", "h">(json_params);

#ifndef JSON_PARAMS
#define JSON_PARAMS json_params
#endif

#define JSON_UNPACK(type, ...) \
    auto [__VA_ARGS__] = json_unpack<type, TO_STRINGS(__VA_ARGS__)>(JSON_PARAMS)

using json = nlohmann::json;
using namespace std;

template <typename T, string_consteval... names>
auto json_unpack(json json_obj) {
    return std::tuple((T)json_obj[names.data]...);
}

template <typename T, typename U>
struct json_unpack_s {};

template <template <typename...> class ContainerTypes, typename... Types, template <string_consteval...> class ContainerNames, string_consteval... names>
struct json_unpack_s<ContainerTypes<Types...>, ContainerNames<names...>> {
    inline static auto function(json json_obj) {
        return std::tuple((Types)json_obj[names.data]...);
    }
};

template<typename T, typename U>
inline auto json_unpack(json json_obj) {
    return json_unpack_s<T, U>::function(json_obj);
}

