#ifndef NW_C_H
#define NW_C_H

#include "../Misc/Types.hpp"

void NW(unsigned int str_cnt_1,char* str1,
			   unsigned int str_cnt_2,char* str2,
			   short match_score,short missmatch_pen,short gap_pen,
			   Alignment& alignment);

#endif // NW_C_H