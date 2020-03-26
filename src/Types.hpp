#ifndef TYPES_H
#define TYPES_H

struct String{
  unsigned int length;
  char* sequence;
};

struct Parameters{
  String* algorithm;
  short match;
  short missmatch;
  short gap;
};

struct Result{
  String* sequence_1;
  String* sequence_2;
  short score;
};

struct Metrics{
  // long long cycles;
  
};

struct Alignment{
  String* sequence_1;
  String* sequence_2;
  Parameters* parameters;
  Result* result;
};

extern "C" String* new_string(char *s);
void destroy_string(String* seq);
Result* new_result();
void destroy_result(Result* result);
Alignment* new_alignment();
void destroy_alignment(Alignment* a);

struct AlignmentMatrix{
  short* matrix;
};

AlignmentMatrix* new_alignment_matrix(unsigned int vector_length,unsigned int seq_length_1,unsigned int seq_length_2);

#endif