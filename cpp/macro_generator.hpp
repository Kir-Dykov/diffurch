// macro_generator.hpp

#include "library/macros.hpp"

#ifndef MACRO_GENERATOR_H
#define MACRO_GENERATOR_H

// Define the macro expansion logic here
#define CONCAT(MACRO_, 1) 4

// Control the iteration using macro redefinition
#if !defined(ITERATION_COUNT)
    #define ITERATION_COUNT 1
#elif ITERATION_COUNT == 1
    #undef ITERATION_COUNT
    #define ITERATION_COUNT 2
#elif ITERATION_COUNT == 2
    #undef ITERATION_COUNT
    #define ITERATION_COUNT 3
#elif ITERATION_COUNT == 3
    #undef ITERATION_COUNT
    #define ITERATION_COUNT 4
#elif ITERATION_COUNT == 4
    #undef ITERATION_COUNT
    #define ITERATION_COUNT 5
#elif ITERATION_COUNT == 5
    #undef ITERATION_COUNT
    #define ITERATION_COUNT 6
#endif

// Generate macros based on the iteration count
#if ITERATION_COUNT <= 5
    DEFINE_MACRO(ITERATION_COUNT);
    #include "macro_generator.hpp"  // Re-include the header to continue the iteration
#endif

#endif  // MACRO_GENERATOR_H
