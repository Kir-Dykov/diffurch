#pragma once

#include "time_series.hpp"


template <typename T, typename Arg>
concept has_Return = requires (T t, Arg arg) {
    t.Return(arg);
};

template <typename T, typename Arg>
concept has_event_function = requires (T t, Arg arg) {
    { t.event_function(arg) } -> std::same_as<void>;
};

template <typename T, typename Arg>
concept has_discontinuity_event_function = requires (T t, Arg arg) {
    { t.discontinuity_event_function(arg) } -> std::same_as<void>;   
};



// template <typename ReturnHandler1, typename ReturnHandler2>
// struct ReturnHandlerUnion {
//     ReturnHandler1 rh1;
//     ReturnHandler2 rh2;
    
//     ReturnHandlerUnion(const ReturnHandler1& rh1, const ReturnHandler2& rh2) : rh1(rh1), rh2(rh2) {};
    
//     template <size_t n, size_t phi_derivatives_n, 
//     typename Enable = enable_if_t<
//            has_Return<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>> 
//         || has_Return<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>>>
//     auto Return(RK_TimeSeries<n, phi_derivatives_n>& TS) { 
//         if constexpr (has_Return<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>> &&
//                       has_Return<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>) {
//             r1.event_function(TS);
//             return make_tuple(r1.Return(TS), r2.Return(TS));
//         } else if constexpr (has_Return<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>>) {
//             return r1.Return(TS);
//         } else if constexpr (has_Return<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>) {
//             return r2.Return(TS);
//         }
//     };
    
    
//     template <size_t n, size_t phi_derivatives_n, 
//     typename = typename enable_if_t<has_event_function<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>> || has_event_function<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>>>
//     auto event_function(RK_TimeSeries<n, phi_derivatives_n>& TS) {
//         if constexpr (has_event_function<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>>) {
//             r1.event_function(TS);
//         }
//         if constexpr (has_event_function<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>) {
//             r2.event_function(TS);
//         }
//     }
    
//     template <size_t n, size_t phi_derivatives_n, 
//     typename = typename enable_if_t<has_discontinuity_event_function<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>> 
//                                  || has_discontinuity_event_function<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>>>
//     auto discontinuity_event_function(RK_TimeSeries<n, phi_derivatives_n>& TS) {
//         if constexpr (has_discontinuity_event_function<ReturnHandler1, RK_TimeSeries<n, phi_derivatives_n>>) {
//             r1.discontinuity_event_function(TS);
//         }
//         if constexpr (has_discontinuity_event_function<ReturnHandler2, RK_TimeSeries<n, phi_derivatives_n>>) {
//             r2.discontinuity_event_function(TS);
//         }
//     }
    
// }

// TODO, define UNION for Return handlers:
// if one doesn't have of the methods, type 



