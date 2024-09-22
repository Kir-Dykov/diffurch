#pragma once

#include <tuple>
#include <type_traits>

#define TRUE_AUTO_EQ(variable_name, expression) decltype(expression) variable_name = expression



template <typename T, typename... Args>
std::function<T(const Args& ...)> FuncTupleToArgsTransform(std::function<T(std::tuple<Args...>)>& func) {
    return [func](Args... args) {
        return func(std::tuple<Args...>{args...});
    };
}





template <auto CaseValue, typename CaseType>
struct case_pair {
    static const auto value = CaseValue;
    using type = CaseType;
};

// Forward declaration
template <auto Value, typename Default, typename... Cases>
struct switch_case;

template <auto Value, typename Default, typename Case, typename... Rest>
struct switch_case<Value, Default, Case, Rest...> {
    using type = std::conditional_t<Value == Case::value, typename Case::type, typename switch_case<Value, Default, Rest...>::type>;
};

template <auto Value, typename Default>
struct switch_case<Value, Default> {
    using type = Default;
};


// Helper alias to simplify usage
template <auto Value, typename... CasesAndTypes>
using switch_case_t = typename switch_case<Value, CasesAndTypes...>::type;





// template <typename T>
// struct func_spec;

// template <typename RetType, typename... ArgTypes>
// struct func_spec<function<RetType(ArgTypes...)>> {
//     using return_type = RetType;
//     using argument_types = tuple<ArgTypes...>;
// };





template <typename T, typename Tuple>
struct tuple_prepend;

template <typename T, typename... Types>
struct tuple_prepend<T, tuple<Types...>> {
    using type = tuple<T, Types...>;
};



// Primary template for ArgTuple which collects types and counts void/non-void types
template <typename... Types>
struct ArgTuple;


// Specialization for the case when the first type is not void
template <typename T, typename... Rest>
struct ArgTuple<T, Rest...> {
    using filtered_type = typename tuple_prepend<T, typename ArgTuple<Rest...>::filtered_type>::type;

    static constexpr std::size_t void_count    =      ArgTuple<Rest...>::void_count;
    static constexpr std::size_t non_void_count = 1 + ArgTuple<Rest...>::non_void_count;

    static constexpr std::array<bool, 1 + sizeof...(Rest)> type_is_non_void = [] {
        std::array<bool, 1 + sizeof...(Rest)> arr = { true };
        std::copy(ArgTuple<Rest...>::type_is_non_void.begin(), ArgTuple<Rest...>::type_is_non_void.end(), arr.begin() + 1);
        return arr;
    }();
    
    static constexpr std::array<size_t, 1 + sizeof...(Rest)> index = [] {
        std::array<size_t, 1 + sizeof...(Rest)> arr = {}; // initialized with zeros
        
        size_t i = 0;
        
        // make `i` index the first non-void type
        while (!type_is_non_void[i] && i <= sizeof...(Rest)) { i++; }; 
        
        for ( i = i + 1 ; i <= sizeof...(Rest); i++) {
            arr[i] = arr[i-1] + type_is_non_void[i];
        }
        return arr;
    }();
    
    
    // helper template function
    template <std::size_t... I>
    static auto unpack_by_sequence(std::index_sequence<I...>, filtered_type& args) {
        return std::tie(std::get<index[I]>(args)...);
    }
    
    static auto unpack(filtered_type& args) {
        return unpack_by_sequence(std::make_index_sequence<1 + sizeof...(Rest)>{}, args);
    }
    
};

// Specialization for the case when the first type is void
template <typename... Rest>
struct ArgTuple<void, Rest...> {
    using filtered_type = typename ArgTuple<Rest...>::filtered_type;

    static constexpr std::size_t void_count     = 1 + ArgTuple<Rest...>::void_count;
    static constexpr std::size_t non_void_count =     ArgTuple<Rest...>::non_void_count;
    
    static constexpr std::array<bool, 1 + sizeof...(Rest)> type_is_non_void = [] {
        std::array<bool, 1 + sizeof...(Rest)> arr = { false };
        std::copy(ArgTuple<Rest...>::type_is_non_void.begin(), ArgTuple<Rest...>::type_is_non_void.end(), arr.begin() + 1);
        return arr;
    }();
    
    static constexpr std::array<size_t, 1 + sizeof...(Rest)> index = [] {
        std::array<size_t, 1 + sizeof...(Rest)> arr = {}; // initialized with zeros
        
        size_t i = 0;
        
        // make `i` index the first non-void type
        while (!type_is_non_void[i] && i <= sizeof...(Rest)) { i++; }; 
        
        for ( i = i + 1 ; i <= sizeof...(Rest); i++) {
            arr[i] = arr[i-1] + type_is_non_void[i];
        }
        return arr;
    }();
    
    // helper template function
    template <std::size_t... I>
    static auto unpack_by_sequence(std::index_sequence<I...>, filtered_type& args) {
        return std::tie(std::get<index[I]>(args)...);
    }
    
    static auto unpack(filtered_type& args) {
        return unpack_by_sequence(std::make_index_sequence<1 + sizeof...(Rest)>{}, args);
    }
    
};

// Specialization for the base case (empty pack)
template <>
struct ArgTuple<> {
    using filtered_type = std::tuple<>;

    static constexpr std::size_t void_count = 0;
    static constexpr std::size_t non_void_count = 0;

    static constexpr std::array<bool, 0> type_is_non_void = {};
    
    static constexpr std::array<size_t, 0> index = {};
    
//     // helper template function
//     template <std::size_t... I>
//     static auto unpack_by_sequence(std::index_sequence<I...>, const filtered_type& args) {
//         return std::tie(std::get<index[I]>(args)...);
//     }
    
//     static auto unpack(const filtered_type& args) {
//         return unpack_by_sequence(std::make_index_sequence<1 + sizeof...(Rest)>{}, args);
//     }
    
};