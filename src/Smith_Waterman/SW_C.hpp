#ifndef SW_C_H
#define SW_C_H

#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include <string>
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

void SW(Alignment& alignment,bool debug);

//se le pasa por parametro la secuencia a alinear y los scores
Alignment* alignment_by_SW(std::string implementation,char * sequence_1, char* sequence_2, short gap, short missmatch, short match);

#endif // SW_C_H