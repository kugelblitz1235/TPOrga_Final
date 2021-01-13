#ifndef SW_C_SSE_HPP
#define SW_C_SSE_HPP
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW{
namespace C {
namespace SSE{
    //Version paralelizada, con el set de instrucciones SSE, del algoritmo de alineamiento global Smith-Waterman
    void SW (Alignment& alignment, bool debug);
}
}
}

#endif // SW_C_SSE_HPP