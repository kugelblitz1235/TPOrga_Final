#include <iostream>
#include "Types.h"

using namespace std;


void SW(unsigned int str_cnt_1,char* str1,
			   unsigned int str_cnt_2,char* str2,
			   short match_score,short missmatch_pen,short gap_pen,
			   Alignment& alignment){
	
	unsigned int best_y = 0;
	unsigned int best_x = 0;
	
	short** scores = (short**)malloc((str_cnt_1+1)*sizeof(short*));
	
	for(unsigned int y = 0;y < str_cnt_1+1;y++)
		scores[y] = (short*)malloc((str_cnt_2+1)*sizeof(short));
	
	scores[0][0] = 0;
	
	for(unsigned int y = 1;y < str_cnt_1+1;y++){
		scores[y][0] = scores[y-1][0] + gap_pen;
	}
	
	for(unsigned int x = 1;x < str_cnt_2+1;x++){
		scores[0][x] = scores[0][x-1] + gap_pen;
	}
	
	for(unsigned int y = 1;y < str_cnt_1+1;y++){
			
		for(unsigned int x = 1;x < str_cnt_2+1;x++){
			
			short score_diag = scores[y-1][x-1] + 
						  missmatch_pen*(str1[y-1]!= str2[x-1]) + 
						  match_score*(str1[y-1] == str2[x-1]);
			
			short score_left = scores[y][x-1] + gap_pen;
			
			short score_up = scores[y-1][x] + gap_pen;
			
			short best_score = 0;
			
			best_score = max(best_score,score_diag);
			
			best_score = max(best_score,score_left);
			
			best_score = max(best_score,score_up);
			
			scores[y][x] = best_score;
			
			if(scores[y][x] > scores[best_y][best_x]){
				best_y = y;
				best_x = x;
			}
		}
	}
	
	for(unsigned int y = 0;y < str_cnt_1+1;y++){
		for(unsigned int x = 0;x < str_cnt_2+1;x++){
			cerr << (int)scores[y][x] << " ";
		}
		cerr << endl << endl;
	}	
	
	char* best_sequence_1 = (char*)malloc(str_cnt_1+1+str_cnt_2+1);
	char* best_sequence_2 = (char*)malloc(str_cnt_1+1+str_cnt_2+1);
	
	unsigned int length = 0;
	
	unsigned int y = best_y;
	unsigned int x = best_x;
	
	while(y != 0 || x != 0){
		
		//cerr << x << " " << y << " " << scores[y][x] << endl;
		
		if(y > 0 && x > 0){
			short score_diag = scores[y-1][x-1] + 
						  missmatch_pen*(str1[y-1]!= str2[x-1]) + 
						  match_score*(str1[y-1] == str2[x-1]);
			
			short score_left = scores[y][x-1] + gap_pen;
			
			short score_up = scores[y-1][x] + gap_pen;
			
			if(scores[y][x] == 0){
				break;
			}else if(score_diag == scores[y][x]){
				best_sequence_1[length] = str1[y-1];
				best_sequence_2[length] = str2[x-1];
				x--;
				y--;
			}else if(score_up == scores[y][x]){
				best_sequence_1[length] = str1[y-1];
				best_sequence_2[length] = '_';
				y--;
			}else {
				best_sequence_1[length] = '_';
				best_sequence_2[length] = str2[x-1];
				x--;
			}
			
		}else if(x > 0){
			best_sequence_1[length] = '_';
			best_sequence_2[length] = str2[x-1];
			x--;
		}else if(y > 0){
			best_sequence_1[length] = str1[y-1];
			best_sequence_2[length] = '_';
			y--;
		}
		
		length++;
	}
	
	for(int i = 0;i < length/2;i++){
		char swap = best_sequence_1[i];
		best_sequence_1[i] = best_sequence_1[length-1-i];
		best_sequence_1[length-1-i] = swap;
		swap = best_sequence_2[i];
		best_sequence_2[i] = best_sequence_2[length-1-i];
		best_sequence_2[length-1-i] = swap;
	}
	
	for(int i = 0;i < length;i++){
		cerr << best_sequence_1[i];
	}cerr << endl;
	for(int i = 0;i < length;i++){
		cerr << best_sequence_2[i];
	}cerr << endl;
		
	for(unsigned int y = 0;y < str_cnt_1+1;y++)
		free(scores[y]);
	
	free(scores);
	
	alignment.length = length;
	alignment.sequence_1 = best_sequence_1;
	alignment.sequence_2 = best_sequence_2;
}

