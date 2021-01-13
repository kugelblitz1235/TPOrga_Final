#ifndef SW_C_H
#define SW_C_H

#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include "SW_C_LIN.hpp"
#include "SW_C_SIMDlogic.hpp"
#include "SW_C_SSE.hpp"
#include "SW_C_AVX.hpp"
#include "SW_C_AVX512.hpp"
#include <string>
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW{
    //se le pasa por parametro la secuencia a alinear y los scores
    Alignment* get_alignment(std::string implementation,char * sequence_1, char* sequence_2, short gap, short missmatch, short match);
}


#endif // SW_C_H