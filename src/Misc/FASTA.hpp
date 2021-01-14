#ifndef FASTA_H
#define FASTA_H

#include <fstream>
#include <iostream>
#include "Types.hpp"

using namespace std;

void FASTA_to_alignment(Alignment* alignment, const char* file1, const char* file2);

#endif
