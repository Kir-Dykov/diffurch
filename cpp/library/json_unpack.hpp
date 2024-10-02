#pragma once

#include "json.hpp"
#include "string_consteval.hpp"
#include "tuple.hpp"

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

