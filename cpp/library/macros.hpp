#pragma once

#define STR2(x) #x
#define STR(x) STR2(x)

#define CAT2(a, b) a##b
#define CAT(a, b) CAT2(a,b)


// macro that replaces auto keyword for non-static member variables
// roughly equivalent to `auto name = value`
#define AUTO(name, value) decltype(value) name = value 
