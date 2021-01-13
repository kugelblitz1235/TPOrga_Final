#ifndef SW_C_SIMDlogic_HPP
#define SW_C_SIMDlogic_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW{
namespace C {
namespace SIMDlogic{
    //Version que emula la paralelizacion para un tamaño de registro dado, del algoritmo de alineamiento global Smith-Waterman
    void SW (Alignment& alignment, int vector_len, bool debug);
}
}
}

#endif // SW_C_SIMDlogic_HPP