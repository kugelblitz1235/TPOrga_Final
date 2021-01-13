#ifndef NW_C_AVX512_HPP
#define NW_C_AVX512_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW{
namespace C {
namespace AVX512{
    //Version paralelizada, con el set de instrucciones AVX512, del algoritmo de alineamiento global Needleman-Wunsch
    void NW (Alignment& alignment, bool debug);
}
}
}

#endif // NW_C_AVX512_HPP