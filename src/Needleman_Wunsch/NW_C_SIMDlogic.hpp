#ifndef NW_C_SIMDlogic_HPP
#define NW_C_SIMDlogic_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW{
namespace C {
namespace SIMDlogic{
    //Version que emula la paralelizacion para un tamaño de registro dado, del algoritmo de alineamiento global Needleman-Wunsch
    void NW (Alignment& alignment, int vector_len, bool debug);
}
}
}

#endif // NW_C_SIMDlogic_HPP