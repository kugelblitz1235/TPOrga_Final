#ifndef TYPES_H
#define TYPES_H

struct Sequence{
  unsigned int length;
  char* sequence;
};

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

struct Metrics{
  // long long cycles;
  
};

struct Alignment{
  Sequence* sequence_1;
  Sequence* sequence_2;
  Parameters* parameters;
  Result* result;
};


#endif