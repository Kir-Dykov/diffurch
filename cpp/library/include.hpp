#pragma once

#include "de/differential_equation.hpp"
#include "de/return_handler_list.hpp"
#include "de/step_size_controller.hpp"

#include "io/json_unpack.hpp"
#include "io/json.hpp"
#include "io/param_names.hpp"
#include "io/print.hpp"
#include "io/progress_bar.hpp"
#include "io/save.hpp"

#include "utils/find_root.hpp"
#include "utils/macros.hpp"
#include "utils/math.hpp"
#include "utils/polynomial.hpp"
#include "utils/real.hpp"
#include "utils/string_consteval.hpp"
#include "utils/tuple.hpp"
#include "utils/vec.hpp"


#include "equations/ode.hpp"
#include "equations/dde.hpp"
#include "equations/relay_dde.hpp"





