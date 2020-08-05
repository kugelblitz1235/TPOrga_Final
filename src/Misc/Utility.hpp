#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <string>
#include <string.h>

#include "Types.hpp"

using namespace std;

void printSpaceLine(int dataSize,int columns);
void printDivisionLine(int dataSize,int columns);
void printScoreMatrix(short* matrix,Alignment* alignment,int vec);

#endif