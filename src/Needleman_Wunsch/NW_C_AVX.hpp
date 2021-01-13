#ifndef NW_C_AVX_HPP
#define NW_C_AVX_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW{
namespace C {
namespace AVX{
    //Version paralelizada, con el set de instrucciones AVX, del algoritmo de alineamiento global Needleman-Wunsch
    void NW (Alignment& alignment, bool debug);
}
}
}

#endif // NW_C_AVX_HPP