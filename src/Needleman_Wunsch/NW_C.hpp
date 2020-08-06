#ifndef NW_C_H
#define NW_C_H
#include "../Misc/Types.hpp"
#include "../Misc/AlignAlgo.hpp"

void NW (Alignment& alignment);

//se le pasa por parametro la secuencia a alinear y los scores
Alignment* alignment_by_NW(char * sequence_1, char* sequence_2, short gap, short missmatch, short match);

#endif // NW_C_H