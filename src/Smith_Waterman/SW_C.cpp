#include <iostream>
#include "../Misc/Types.cpp"

using namespace std;

void find_scores(
	short **scores,
	Alignment& alignment,
	unsigned int &best_y,unsigned int &best_x,
	bool debug = false
){
	RNA_Sequence* seq1 = alignment.sequence_1;
	RNA_Sequence* seq2 = alignment.sequence_2;
	
	for(unsigned int y = 1;y < alignment.sequence_1->length+1;y++){
		for(unsigned int x = 1;x < alignment.sequence_2->length+1;x++){
			
			short score_left = scores[y][x-1] + alignment.parameters->gap;
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(seq1->sequence[y-1] != seq2->sequence[x-1]) + 
						  alignment.parameters->match*(seq1->sequence[y-1] == seq2->sequence[x-1]);
			
			short best_score = 0;
			
			best_score = max(best_score,score_left);
			best_score = max(best_score,score_up);
			best_score = max(best_score,score_diag);
			
			scores[y][x] = best_score;
			
			if(scores[y][x] > scores[best_y][best_x]){
				best_y = y;
				best_x = x;
			}
		}
	}

	if(debug){
		for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++){
			for(unsigned int x = 0;x < alignment.sequence_2->length+1;x++){
				cerr << (int)scores[y][x] << " ";
			}
			cerr << endl << endl;
		}	
	}

}

void backtrack_solution(
	short **scores,
	Alignment& alignment,
	unsigned int best_x,unsigned int best_y,
	bool debug = false
){
	
	RNA_Sequence* seq1 = alignment.sequence_1;
	RNA_Sequence* seq2 = alignment.sequence_2;
	
	RNA_Sequence* best_seq_1 = alignment.result->sequence_1;
	RNA_Sequence* best_seq_2 = alignment.result->sequence_2;

	unsigned int length = 0;
	
	unsigned int y = best_y;
	unsigned int x = best_x;
	
	while(y != 0 || x != 0){
		
		if(debug)
			cerr << x << " " << y << " " << scores[y][x] << endl;
		
		if(y > 0 && x > 0){

			short score_left = scores[y][x-1] + alignment.parameters->gap;
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(seq1->sequence[y-1] != seq2->sequence[x-1]) + 
						  alignment.parameters->match*(seq1->sequence[y-1] == seq2->sequence[x-1]);
			
			if(scores[y][x] == 0){
				break;
			}else if(score_diag == scores[y][x]){
				best_seq_1->sequence[length] = seq1->sequence[y-1];
				best_seq_2->sequence[length] = seq2->sequence[x-1];
				x--;
				y--;
			}else if(score_up == scores[y][x]){
				best_seq_1->sequence[length] = seq1->sequence[y-1];
				best_seq_2->sequence[length] = '_';
				y--;
			}else {
				best_seq_1->sequence[length] = '_';
				best_seq_2->sequence[length] = seq2->sequence[x-1];
				x--;
			}
			
		}else if(x > 0){
			best_seq_1->sequence[length] = '_';
			best_seq_2->sequence[length] = seq2->sequence[x-1];
			x--;
		}else if(y > 0){
			best_seq_1->sequence[length] = seq1->sequence[y-1];
			best_seq_2->sequence[length] = '_';
			y--;
		}
		
		length++;
	}

	//creamos los nuevos strings y destruimos los anteriores
	best_seq_1 = new_RNA_Sequence(best_seq_1->sequence);
	best_seq_2 = new_RNA_Sequence(best_seq_2->sequence);
	destroy_RNA_Sequence(alignment.result->sequence_1);
	destroy_RNA_Sequence(alignment.result->sequence_2);
	alignment.result->sequence_1 = best_seq_1;
	alignment.result->sequence_2 = best_seq_2;
}

void SW(
	Alignment& alignment,
	bool debug
){
	unsigned int best_y = 0;
	unsigned int best_x = 0;
	
	//Inicializar matriz de scores
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
	
	find_scores(
		scores,
		alignment,
		best_y,best_x,
		debug
	);

	alignment.result->sequence_1 = new_RNA_Sequence_length(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	alignment.result->sequence_2 = new_RNA_Sequence_length(alignment.sequence_1->length+1+alignment.sequence_2->length+1);

	backtrack_solution(
		scores,
		alignment,
		best_x,best_y
	);

	if(debug){
		for(int i = 0;i < alignment.result->sequence_1->length;i++)
			cerr << alignment.result->sequence_1->sequence[i];cerr << endl;
		for(int i = 0;i < alignment.result->sequence_2->length;i++)
			cerr << alignment.result->sequence_1->sequence[i];cerr << endl;
	}
		
	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++)
		free(scores[y]);
	
	free(scores);
}

