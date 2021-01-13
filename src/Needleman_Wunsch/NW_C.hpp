#ifndef NW_C_H
#define NW_C_H

#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include "NW_C_LIN.hpp"
#include "NW_C_SIMDlogic.hpp"
#include "NW_C_SSE.hpp"
#include "NW_C_AVX.hpp"
#include "NW_C_AVX512.hpp"
#include <string>
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW{
    //Método utilizado para seleccionar la versión del algoritmo a ejecutar
    Alignment* get_alignment(std::string implementation, char * sequence_1, char* sequence_2, short gap, short missmatch, short match);
}

#endif // NW_C_Ha