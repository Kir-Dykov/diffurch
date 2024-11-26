#pragma once

// tool for tools

#include <boost/preprocessor.hpp>

#define STR2(x) #x
#define STR(x) STR2(x)

#define CAT2(a, b) a##b
#define CAT(a, b) CAT2(a,b)


// macro that replaces auto keyword for non-static member variables
// roughly equivalent to `auto name = value`
#define AUTO(name, value) decltype(value) name = value 

#define TOSTR(elem) #elem
#define TO_STRINGS_OP(r, data, elem) TOSTR(elem)
#define TO_STRINGS_SEQ(...) BOOST_PP_SEQ_TRANSFORM(TO_STRINGS_OP, ~, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
#define TO_STRINGS(...)  BOOST_PP_TUPLE_REM_CTOR(BOOST_PP_SEQ_TO_TUPLE(TO_STRINGS_SEQ(__VA_ARGS__)))

#define PREPEND_OP(r, data, elem) BOOST_PP_CAT(data, elem)
#define PREPEND_OP_(r, data, elem) data elem
#define PREPEND_SEQ(prefix, ...)  BOOST_PP_SEQ_TRANSFORM(PREPEND_OP, prefix, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
#define PREPEND_SEQ_(prefix, ...) BOOST_PP_SEQ_TRANSFORM(PREPEND_OP_, prefix, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
#define PREPEND(...)   BOOST_PP_TUPLE_REM_CTOR(BOOST_PP_SEQ_TO_TUPLE(PREPEND_SEQ (__VA_ARGS__)))
#define PREPEND_(...)  BOOST_PP_TUPLE_REM_CTOR(BOOST_PP_SEQ_TO_TUPLE(PREPEND_SEQ_(__VA_ARGS__)))

#define APPEND_OP(r, data, elem)  BOOST_PP_CAT(elem, data)
#define APPEND_OP_(r, data, elem) elem data
#define APPEND_SEQ(prefix, ...)  BOOST_PP_SEQ_TRANSFORM(APPEND_OP , prefix, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
#define APPEND_SEQ_(prefix, ...) BOOST_PP_SEQ_TRANSFORM(APPEND_OP_, prefix, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
#define APPEND(...)   BOOST_PP_TUPLE_REM_CTOR(BOOST_PP_SEQ_TO_TUPLE(APPEND_SEQ (__VA_ARGS__)))
#define APPEND_(...)  BOOST_PP_TUPLE_REM_CTOR(BOOST_PP_SEQ_TO_TUPLE(APPEND_SEQ_(__VA_ARGS__)))

#define INITIALIZER_LIST_OP(r, data, elem) elem(elem ## _)
#define INITIALIZER_LIST_SEQ(...) BOOST_PP_SEQ_TRANSFORM(INITIALIZER_LIST_OP, ~, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
#define INITIALIZER_LIST(...)   BOOST_PP_TUPLE_REM_CTOR(BOOST_PP_SEQ_TO_TUPLE(INITIALIZER_LIST_SEQ(__VA_ARGS__)))

#define INIT_STATEMENTS_OP(r, data, elem) elem = BOOST_PP_CAT(elem, _);
#define INIT_STATEMENTS(...) BOOST_PP_SEQ_FOR_EACH(INIT_STATEMENTS_OP, ~, BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))
