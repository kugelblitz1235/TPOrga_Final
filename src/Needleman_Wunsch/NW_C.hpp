#ifndef NW_C_H
#define NW_C_H
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include <string>
#include <limits.h>

void NW_C_LIN (Alignment& alignment);
void NW_C_SSE (Alignment& alignment);
typedef short (*score_fun_t)(short* score_matrix, int seq_row_len, int y,int x, int vector_len);

//se le pasa por parametro la secuencia a alinear y los scores
Alignment* alignment_by_NW(std::string implementation, char * sequence_1, char* sequence_2, short gap, short missmatch, short match);

short get_score_LIN(short* score_matrix, int seq_row_len, int y,int x, int vector_len);
short get_score_SSE(short* score_matrix, int seq_row_len, int y,int x, int vector_len);
void backtrack_solution_NW_C(
	short *score_matrix,
	Alignment& alignment,
	int vector_len,
	unsigned int x,unsigned int y,
	bool SW,
	score_fun_t score_fun,
	bool debug
);

#endif // NW_C_H