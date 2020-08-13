#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <string>
#include <string.h>
#include <emmintrin.h>
#include <immintrin.h>

#include "Types.hpp"

using namespace std;
typedef short (*score_fun_t)(short* score_matrix, int seq_row_len, int y,int x, int vector_len);

void printSpaceLine(int dataSize,int columns);

void printDivisionLine(int dataSize,int columns);

void printScoreMatrix(short* matrix,Alignment* alignment,int vec);

void backtracking_C(
	short *score_matrix,
	Alignment& alignment,
	int vector_len,
	unsigned int x,unsigned int y,
	bool SW,
	score_fun_t score_fun,
	bool debug
);

short get_score_SSE(short* score_matrix, int seq_row_len, int y,int x, int vector_len=4);
short get_score_LIN(short* score_matrix, int seq_row_len, int y,int x, int vector_len=4);

void print128_hex(__m128i var);

void word_to_arr8(short* p,__m128i reg);

void word_to_char8(char* p,__m128i reg);
__m128i char_to_word8(char* p);

#endif