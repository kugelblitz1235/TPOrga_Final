#ifndef SW_C_AVX_HPP
#define SW_C_AVX_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW{
namespace C {
namespace AVX{
    //Version paralelizada, con el set de instrucciones AVX, del algoritmo de alineamiento global Smith-Waterman
    void SW (Alignment& alignment, bool debug);
}
}
}

#endif // SW_C_AVX_HPP