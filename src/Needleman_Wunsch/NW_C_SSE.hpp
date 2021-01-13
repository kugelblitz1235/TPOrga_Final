#ifndef NW_C_SSE_HPP
#define NW_C_SSE_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW{
namespace C {
namespace SSE{
    //Version paralelizada, con el set de instrucciones SSE, del algoritmo de alineamiento global Needleman-Wunsch
    void NW (Alignment& alignment, bool debug);
}
}
}

#endif // NW_C_SSE_HPP