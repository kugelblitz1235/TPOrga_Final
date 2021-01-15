#ifndef UTILITY_H
#define UTILITY_H

#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <string.h>
#include <emmintrin.h>
#include <immintrin.h>

#include "Types.hpp"

#define _mm256_set_m128i(v0, v1)  _mm256_insertf128_si256(_mm256_castsi128_si256(v1), (v0), 1)

using namespace std;
using namespace chrono;

#define NOW high_resolution_clock::now()
#define Now() chrono::high_resolution_clock::now()
static struct Stopwatch {
	chrono::high_resolution_clock::time_point c_time, c_timeout;
	void setTimeout(int us) { c_timeout = c_time + chrono::microseconds(us); }
	void Start() { c_time = Now();}
	inline bool Timeout() { return Now() > c_timeout; }
	long long EllapsedMicroseconds() { return chrono::duration_cast<chrono::microseconds>(Now() - c_time).count(); }
	long long EllapsedMilliseconds() { return chrono::duration_cast<chrono::milliseconds>(Now() - c_time).count(); }
} stopwatch;//} Stopwatch


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

extern "C" void print_registers();

extern "C" void print_xmm(short *stack, unsigned int n);


#endif