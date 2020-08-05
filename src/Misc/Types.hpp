#ifndef TYPES_H
#define TYPES_H

struct Sequence{
  unsigned int length;
  char* sequence; //data
};

extern "C" Sequence* new_Sequence_from_string(char* s);
extern "C" Sequence* new_Sequence(unsigned int length);
void destroy_Sequence(Sequence* seq);
void reverse_Sequence(Sequence* seq);

struct Parameters{
  char* algorithm;
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
};

Alignment* new_alignment();
void destroy_alignment(Alignment* a);

struct AlignmentMatrix{
  short* matrix;
};

AlignmentMatrix* new_alignment_matrix(unsigned int vector_length,unsigned int seq_length_1,unsigned int seq_length_2);

#endif