#ifndef SW_C_H
#define SW_C_H

#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include <string>
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace SW {

    void SW_C_LIN(Alignment& alignment,bool debug);
    void SW_C_withLogicSSE(Alignment& alignment, int vector_len, bool debug);
    void SW_C_SSE(Alignment& alignment, bool debug);
    void SW_C_AVX(Alignment& alignment, bool debug);
    namespace AVX512{
        void SW_C (Alignment& alignment, bool debug);
    }
    //se le pasa por parametro la secuencia a alinear y los scores
    Alignment* alignment_by_SW(std::string implementation,char * sequence_1, char* sequence_2, short gap, short missmatch, short match);

}

#endif // SW_C_H