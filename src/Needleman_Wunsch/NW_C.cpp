#include <iostream>
#include "../Misc/Types.cpp"
#include "NW_C.hpp"

using namespace std;

void NW(Alignment& alignment){
	
	short** scores = (short**)malloc((alignment.sequence_1->length+1)*sizeof(short*));
	
	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++)
		scores[y] = (short*)malloc((alignment.sequence_2->length+1)*sizeof(short));
	
	scores[0][0] = 0;
	
	for(unsigned int y = 1;y < alignment.sequence_1->length+1;y++){
		scores[y][0] = scores[y-1][0] + alignment.parameters->gap;
	}
	
	for(unsigned int x = 1;x < alignment.sequence_2->length+1;x++){
		scores[0][x] = scores[0][x-1] + alignment.parameters->gap;
	}
	
	for(unsigned int y = 1;y < alignment.sequence_1->length+1;y++){
			
		for(unsigned int x = 1;x < alignment.sequence_2->length+1;x++){
			
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(alignment.sequence_1->sequence[y-1]!= alignment.sequence_2->sequence[x-1]) + 
						  alignment.parameters->match*(alignment.sequence_1->sequence[y-1] == alignment.sequence_2->sequence[x-1]);
			
			short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			
			short best_score = max(score_diag,score_left);
			
			best_score = max(best_score,score_up);
			
			scores[y][x] = best_score;
			
		}
	}
	
	cerr << "Score Matrix" << endl;
	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++){
		for(unsigned int x = 0;x < alignment.sequence_2->length+1;x++){
			cerr << (int)scores[y][x] << " ";
		}
		cerr << endl << endl;
	}	
	
	char* best_sequence_1 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	char* best_sequence_2 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	
	unsigned int length = 0;
	
	unsigned int y = alignment.sequence_1->length;
	unsigned int x = alignment.sequence_2->length;
	
	while(y != 0 || x != 0){
		//cerr << x << " " << y << " " << scores[y][x] << endl;
		if(y > 0 && x > 0){
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(alignment.sequence_1->sequence[y-1]!= alignment.sequence_2->sequence[x-1]) + 
						  alignment.parameters->match*(alignment.sequence_1->sequence[y-1] == alignment.sequence_2->sequence[x-1]);
			
			short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			
			if(score_diag == scores[y][x]){
				best_sequence_1[length] = alignment.sequence_1->sequence[y-1];
				best_sequence_2[length] = alignment.sequence_2->sequence[x-1];
				x--;
				y--;
			}else if(score_up == scores[y][x]){
				best_sequence_1[length] = alignment.sequence_1->sequence[y-1];
				best_sequence_2[length] = '_';
				y--;
			}else{
				best_sequence_1[length] = '_';
				best_sequence_2[length] = alignment.sequence_2->sequence[x-1];
				x--;
			}
			
		}else if(x > 0){
			best_sequence_1[length] = '_';
			best_sequence_2[length] = alignment.sequence_2->sequence[x-1];
			x--;
		}else if(y > 0){
			best_sequence_1[length] = alignment.sequence_1->sequence[y-1];
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
	
	cerr << "Best sequences" << endl;
	for(int i = 0;i < length;i++){
		cerr << best_sequence_1[i];
	}cerr << endl;
	for(int i = 0;i < length;i++){
		cerr << best_sequence_2[i];
	}cerr << endl;

	alignment.result->sequence_1 = new_RNA_Sequence(best_sequence_1);
	alignment.result->sequence_2 = new_RNA_Sequence(best_sequence_2);
	alignment.result->score = scores[alignment.sequence_1->length][alignment.sequence_2->length];

	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++)
		free(scores[y]);
	
	free(scores);
	
	// alignment.length = length;
	// alignment.sequence_1 = best_sequence_1;
	// alignment.sequence_2 = best_sequence_2;


}
