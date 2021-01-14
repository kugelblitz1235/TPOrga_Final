#ifndef TYPES_H
#define TYPES_H
#include <iostream>
#include <string.h>

struct Sequence{
  unsigned int length;
  char* sequence; //data
};

extern "C" Sequence* new_Sequence_from_string(char* s);
extern "C" Sequence* new_Sequence(unsigned int length);
void destroy_Sequence(Sequence* seq);
void reverse_Sequence(Sequence* seq);

struct Parameters{
  std::string* algorithm;
  short match;
  short missmatch;
  short gap;
};

struct Result{
  Sequence* sequence_1;
  Sequence* sequence_2;
  short score;
};

Result* new_result();
void destroy_result(Result* result);

struct Metrics{
// long long cycles;
};

struct Alignment{
  Sequence* sequence_1;
  Sequence* sequence_2;
  Parameters* parameters;
  Result* result;
  short* matrix;
};

Alignment* new_alignment();
void destroy_alignment(Alignment* a);

#endif