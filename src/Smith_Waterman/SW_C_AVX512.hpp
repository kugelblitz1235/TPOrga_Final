#ifndef SW_C_AVX512_HPP
#define SW_C_AVX512_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW{
namespace C {
namespace AVX512{
    //Version paralelizada, con el set de instrucciones AVX512, del algoritmo de alineamiento global Smith-Waterman
    void SW (Alignment& alignment, bool debug);
}
}
}

#endif // SW_C_AVX512_HPP