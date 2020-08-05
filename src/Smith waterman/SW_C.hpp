#ifndef SW_C_H
#define SW_C_H

#include "Types.hpp"

void SW(unsigned int str_cnt_1,char* str1,
			   unsigned int str_cnt_2,char* str2,
			   short match_score,short missmatch_pen,short gap_pen,
			   Alignment& alignment);

#endif // SW_C_H