#ifndef NW_C_H
#define NW_C_H
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include <string>
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW {

    void NW_C_LIN (Alignment& alignment, bool debug);
    void NW_C_withLogicSSE (Alignment& alignment, int vector_len, bool debug);
    void NW_C_SSE (Alignment& alignment, bool debug);
    void NW_C_AVX (Alignment& alignment, bool debug);
    namespace AVX512{
        void NW_C (Alignment& alignment, bool debug);
    }
    //se le pasa por parametro la secuencia a alinear y los scores
    Alignment* alignment_by_NW(std::string implementation, char * sequence_1, char* sequence_2, short gap, short missmatch, short match);

}

#endif // NW_C_Ha