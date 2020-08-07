#ifndef NW_C_H
#define NW_C_H
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include <string>
#include <limits.h>

void NW_C_LIN (Alignment& alignment);
void NW_C_SSE (Alignment& alignment);

//se le pasa por parametro la secuencia a alinear y los scores
Alignment* alignment_by_NW(std::string implementation, char * sequence_1, char* sequence_2, short gap, short missmatch, short match);

#endif // NW_C_H