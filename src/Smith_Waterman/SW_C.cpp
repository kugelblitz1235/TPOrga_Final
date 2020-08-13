#include <iostream>
#include "SW_C.hpp"

using namespace std;

void SW(
	Alignment& alignment,
	bool debug
){
	unsigned int seq1_len = alignment.sequence_1->length;
	unsigned int seq2_len = alignment.sequence_2->length;
	char* seq1 = alignment.sequence_1->sequence;
	char* seq2 = alignment.sequence_2->sequence;
	
	short** scores = (short**)malloc((seq2_len)*sizeof(short*));
	
	for(unsigned int y = 0;y < seq2_len;y++)
		scores[y] = (short*)malloc((seq1_len)*sizeof(short));
	
	for(unsigned int y = 0;y < seq2_len;y++)
		scores[y][0] = 0;
	
	for(unsigned int x = 0;x < seq1_len;x++)
		scores[0][x] = 0;
	
	short best_global = 0;
	unsigned int best_x=0;
	unsigned int best_y=0;

	for(unsigned int y = 1;y < seq2_len;y++){
			
		for(unsigned int x = 1;x < seq1_len;x++){
			
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(seq2[y]!= seq1[x]) + 
						  alignment.parameters->match*(seq2[y] == seq1[x]);

			short score_left = scores[y][x-1] + alignment.parameters->gap;
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			
			short best_score = 0;
			best_score = max(best_score,score_left);
			best_score = max(best_score,score_up);
			best_score = max(best_score,score_diag);
			scores[y][x] = best_score;

			if(best_global < best_score){
				best_global = best_score;
				best_y=y;
				best_x=x;
			}
		}
	}
	
	if(debug){
		cerr << "Score Matrix" << endl;
		for(unsigned int y = 0;y < seq2_len;y++){
			for(unsigned int x = 0;x < seq1_len;x++){
				cerr << (int)scores[y][x] << " ";
			}
			cerr << endl << endl;
		}	
	}
	
	backtracking_C(
		(short*)scores,
		alignment,
		1,
		best_x,best_y,
		true,
		(score_fun_t)get_score_LIN,
		false
	);

	
	for(unsigned int y = 0;y < seq2_len;y++)
		free(scores[y]);
	
	free(scores);
}

Alignment* alignment_by_SW(std::string implementation, char * sequence_1, char* sequence_2, short gap, short missmatch, short match){
	//creo la estructura vacia
	Alignment* alignment = new_alignment();
	
	//relleno la estructura con los valores correspondientes
	alignment->sequence_1=new_Sequence_from_string(sequence_1);
	alignment->sequence_2=new_Sequence_from_string(sequence_2);
	(alignment->parameters)->match=match;
	(alignment->parameters)->missmatch=missmatch;
	(alignment->parameters)->gap=gap;

	if(implementation.compare("C") == 0){
		//ejecuto la implementación en c con debug
		SW(*alignment, false);
	
	}else if(implementation.compare("LIN") == 0){
		
		//ejecuto el algoritmo en asm lineal
		SWLIN(alignment);
	}
	else{
		throw "No existe la implementación ingresada.";
	}


	//devuelvo la estructura modificada
	return alignment;
}
