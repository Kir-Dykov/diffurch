#pragma once

#include <type_traits>
#include <tuple>

namespace Events {
    
    template <typename... EventType> 
    struct Events {
    
        tuple<EventType...> events;
        // and filtered 
        
        Events(EventType... events) {
            events = make_tuple(events...);
            // and filtered
        }
    };
    
    template <typename... EventFeature>
    struct Event : std::remove_cvref_t<EventFeature>... {
        Event(EventFeature... event_feature) : std::remove_cvref_t<EventFeature>(event_feature)... {};
    };
    
    template <typename LambdaT>
    struct Save {
        
        // using ValueType = decltype(WhatToSaveLambdaType()); 
        
        vector<???> container;
        
        LambdaT what_to_save;
    
        Save(const LambdaT& lambda) : what_to_save(lambda) {};
        Save(LambdaT&& lambda) : what_to_save(lambda) {}; // Move Constructor
        
        void save(...) {
            container.push_back(what_to_save(...));
        }
        
        auto get() {
            return container;
        }
    };
    
    template <typename LambdaDetectT, typename LambdaFilterT> 
    struct When {
        //
    };
}
// inline auto events() { 
//         using namespace Events;
//         return 
//         Events(
//             Event(
//                 WhenStep{}, 
//                 Save(
//                     [this](const auto& state){
//                         auto& [x,y,z] = state.X;
//                         auto& t = state t;
//                         return make_tuple(t, x, y, z); }
//                 )
//             ),
//             Event(
//                 When([this] STATE_EVAL((x,y,z), (x_<1>(t) - a))),
//                 Save([this] STATE_EVAL((x,y,z), (t)))
//             )
//         )
//     }

// // Base filter implementation: Recurse through variadic types
// template <typename A, typename... Ts>
// struct FilterBaseOf;

// // Specialization for when Ts is empty: return an empty tuple
// template <typename A>
// struct FilterBaseOf<A> {
//     using type = std::tuple<>;
// };

// // Recursive case: filter types inheriting from A
// template <typename A, typename Head, typename... Tail>
// struct FilterBaseOf<A, Head, Tail...> {
// private:
//     // Check if Head inherits from A
//     using TailFiltered = typename FilterBaseOf<A, Tail...>::type;

// public:
//     using type = std::conditional_t<
//         std::is_base_of_v<A, Head>,         // If Head is derived from A
//         decltype(std::tuple_cat(std::tuple<Head>{}, TailFiltered)),  // Include Head
//         TailFiltered                        // Otherwise, skip Head
//     >;
// };



// // concept isEventFeature

// namespace Event {
//     template <typename... EventFeatures>
//     struct Event : EventFeatures... {};

//     template <typename... EventTypes>
//     struct Events {
//         static tuple<EventTypes...> all_events;
//         static filter_types<inherits<WhenStep>, tuple<EventTypes...>> step_events;
//         static filter_types<inherits<WhenRejectedStep>, tuple<EventTypes...>> rejected_step_events;
//         static filter_types<inherits<WhenStopIngetgration>, tuple<EventTypes...>> stop_integration_events;
//         // static filter_types<inherits<>, tuple<EventTypes...>> step_events;
//     }




    struct WhenStep {};
    
    struct WhenRejectedStep {};

    struct WhenStopIngetgration {};

    template <typename LambdaDetect, typename LambdaFilter = void>
    // requires requires 
    struct When {
        bool detect(const auto& state);
        bool locate(const auto& state);
    };
    
    


    template <typename LambdaWhatToSave>
    struct Save {

        vector<???> container;

        auto save(const auto& state);

        auto get();
    }

    template <typename LambdaChange, typename LambdaChangeAbove = void, typename LambdaChangeBelow = void>
    struct Change {
        void change(auto& state); // in state class, make everything, except X and Z, read-only or completely hidden
    }
    
    template <typename LambdaCall>
    struct Call {
        void call(const auto& state); // in state class, make everything, except X and Z, read-only or completely hidden
    }
    
    struct CallStopIntegration {};
    
    
    
    
    
    struct SaveSteps {};
}

// Events<
//        Event<
//             Event::When<[](){}, [](){}>,
//             Event::Save<[](){}>,
//             Event::Change<[](){}>,
//             Event::Call
//        >,
//        Event<
//             Event::When<[](){}, [](){}>,
//             Event::Save<[](){}>,
//             Event::Change<[](){}>,
//             Event::Call
//        >,
// >