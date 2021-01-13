#ifndef SW_C_LIN_HPP
#define SW_C_LIN_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW{
namespace C {
namespace LIN{
    //Version lineal del algoritmo de alineamiento global Smith-Waterman
    void SW (Alignment& alignment, bool debug);
}
}
}

#endif // SW_C_LIN_HPP