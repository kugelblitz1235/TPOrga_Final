#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <emmintrin.h>
#include <immintrin.h>

#include "Types.hpp"

using namespace std;
typedef short (*score_fun_t)(short* score_matrix, unsigned int seq_row_len, int y,int x, int vector_len);

void printStartLine(int row,int vec);
void printEndLine(int row,int vec);

void printScoreMatrix(short* matrix, Alignment* alignment, int vector_len, ofstream &ofs);

extern "C" void backtracking_C(
	short *score_matrix,
	Alignment& alignment,
	int vector_len,
	unsigned int x,unsigned int y,
	bool SW,
	score_fun_t score_fun,
	bool debug
);

extern "C" short get_score_SSE(short* score_matrix, unsigned int seq_row_len, int y,int x, int vector_len=8);
short get_score_LIN(short* score_matrix, unsigned int seq_row_len, int y,int x, int vector_len=8);

void print_arr(short* p, int len);

void print128_hex(__m128i var);

void reg_to_arr(short* p,__m128i reg);

void word_to_arr8(short* p,__m128i reg);

void word_to_char8(char* p,__m128i reg);
__m128i char_to_word8(char* p);

bool valid_alignment(Alignment& alignment);

bool check_scr_matrix_manual(short ** valid_matrix, Alignment* alignment, score_fun_t score_fun, int vector_len=8, bool verbose=true);

#endif