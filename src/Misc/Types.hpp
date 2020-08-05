#ifndef TYPES_H
#define TYPES_H

struct RNA_Sequence{
  unsigned int length;
  char* sequence;
};

extern "C" RNA_Sequence* new_RNA_Sequence(char* s);
extern "C" RNA_Sequence* new_RNA_Sequence_length(unsigned int length);
void destroy_RNA_Sequence(RNA_Sequence* seq);
void reverse_RNA_Sequence(RNA_Sequence* seq);

struct Parameters{
  RNA_Sequence* algorithm;
  short match;
  short missmatch;
  short gap;
};

struct Result{
  RNA_Sequence* sequence_1;
  RNA_Sequence* sequence_2;
  short score;
};

Result* new_result();
void destroy_result(Result* result);

struct Metrics{
// long long cycles;
};

struct Alignment{
  RNA_Sequence* sequence_1;
  RNA_Sequence* sequence_2;
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