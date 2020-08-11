#include <iostream>
#include "NW_C.hpp"

#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;


void NW_C_LIN(Alignment& alignment){
	
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
			
			//short score_left = scores[y][x-1] + alignment.parameters->gap;
			
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
	
	for(unsigned int i = 0;i < length/2;i++){
		char swap = best_sequence_1[i];
		best_sequence_1[i] = best_sequence_1[length-1-i];
		best_sequence_1[length-1-i] = swap;
		swap = best_sequence_2[i];
		best_sequence_2[i] = best_sequence_2[length-1-i];
		best_sequence_2[length-1-i] = swap;
	}
	
	cerr << "Best sequences" << endl;
	for(unsigned int i = 0;i < length;i++){
		cerr << best_sequence_1[i];
	}cerr << endl;
	for(unsigned int i = 0;i < length;i++){
		cerr << best_sequence_2[i];
	}cerr << endl;

	alignment.result->sequence_1 = new_Sequence_from_string(best_sequence_1);
	alignment.result->sequence_2 = new_Sequence_from_string(best_sequence_2);
	alignment.result->score = scores[alignment.sequence_1->length][alignment.sequence_2->length];

	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++)
		free(scores[y]);
	
	free(scores);
	
	// alignment.length = length;
	// alignment.sequence_1 = best_sequence_1;
	// alignment.sequence_2 = best_sequence_2;


}

Alignment* alignment_by_NW(std::string implementation, char* sequence_1, char* sequence_2, short gap, short missmatch, short match){
	//creo la estructura vacia
	Alignment* alignment = new_alignment();
	
	//relleno la estructura con los valores correspondientes
	alignment->sequence_1 = new_Sequence_from_string(sequence_1);
	alignment->sequence_2 = new_Sequence_from_string(sequence_2);
	(alignment->parameters)->match = match;
	(alignment->parameters)->missmatch = missmatch;
	(alignment->parameters)->gap = gap;

	if(implementation.compare("C") == 0){
		//ejecuto la implementación en c
		NW_C_LIN(*alignment);
	
	}else if(implementation.compare("LIN") == 0){
		
		//ejecuto el algoritmo en asm lineal
		NWLIN(alignment);

	}
	else{
		throw "No existe la implementación ingresada.";
	}

	//devuelvo la estructura modificada
	return alignment;
}



void NW_C_SSE (Alignment& alignment){
	char* seq1 = alignment.sequence_1->sequence;	
	char* seq2 = alignment.sequence_2->sequence;
	unsigned int seq1_len = alignment.sequence_1->length;
	unsigned int seq2_len = alignment.sequence_2->length;
	
	
	//en este caso hardcodeamos el tamaño del vector
	int vector_len = 4;
	
	int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
	int width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
	int score_matrix_sz = height * width * vector_len; // largo de diagonal
		
	short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
	
	short* v_aux = (short*)malloc((width-1)*sizeof(short));
	//llenamos el vector auxiliar
	for(int i = 0;i < width-1;i++){
		v_aux[i] = SHRT_MIN/2;
	}

	int count=0;
	for(int i = 0 ; i < height ; i++){
		unsigned int offset_y = i * width * vector_len;
		for( int j = 0; j < 2 ; j++){
			unsigned int offset_x = j * vector_len;
			//emulamos simd
			for( int k = 0;k < vector_len;k++){
				if( j==1 && k == vector_len-1)
					score_matrix[offset_y + offset_x + k] = i * vector_len * alignment.parameters->gap;
				else
					score_matrix[offset_y + offset_x + k] = SHRT_MIN/2;
				count++;
			}			
		}
	}
	/******************************************************************************************************/
	//arrays auxiliares para el calculo de los scores
	char* str_row = (char*)malloc(vector_len * sizeof(char)); 
	char* str_col = (char*)malloc(vector_len * sizeof(char));
	short* left_score = (short*)malloc(vector_len * sizeof(short));
	short* up_score =(short*)malloc(vector_len * sizeof(short));
	short* diag_score =(short*)malloc(vector_len * sizeof(short));
	short* constant_gap = (short*)malloc(vector_len * sizeof(short));
	short* cmp_missmatch = (short*)malloc(vector_len * sizeof(short));
	short* cmp_match = (short*)malloc(vector_len * sizeof(short));

	for( int k = 0;k < vector_len;k++){
		constant_gap[k] = alignment.parameters->gap;
	}

	for( int i = 0 ; i < height ; i++){
		int offset_y = i * width * vector_len;

		if((i+1)*vector_len >= seq2_len){
			int offset_col = (i+1)*vector_len - seq2_len;
			//simd : leer de memoria (movdqu)
			//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más 
			for( int k = 0;k < vector_len;k++){
				str_col[vector_len - 1 - k] = seq2[i * vector_len - offset_col + k];
			}

			//simd : shift
			for(int s = 0;s < offset_col; s++){
				for(int k = vector_len-2;k >= 0;k--){
					str_col[k+1] = str_col[k];
				}
			}


		}else{
			//simd : leer de memoria (movdqu)
			//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
			for( int k = 0;k < vector_len;k++){
				str_col[vector_len - 1 - k] = seq2[i * vector_len + k];
			}
		}

		//simd : shuffle
		for (int i = 0 ; i < vector_len / 2; i++){
			char temp = str_col[i];
			str_col[i] = str_col[vector_len - i - 1];
			str_col[vector_len - i - 1] = temp;
		}

		for( int j = 2; j < width ; j++){
			int offset_x = j * vector_len;
			//emulamos simd
			if(j-vector_len < 0){ //desborde por izquierda
				//simd : desplazamiento de puntero y levantar datos de memoria
				//levantamos el string con el puntero offseteado para no acceder afuera
				
				int offset_str_row = vector_len - j;
				//simd : leer de memoria (movdqu)
				//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
				for(int k = 0;k < vector_len;k++){
			//		cerr<<k<<endl;
					//Esto se simplifica haciendo seq1[k], porque siempre queremos
					//empezar desde el principio del string por izquierda
					str_row[vector_len - 1 - k] = seq1[j - vector_len + offset_str_row + k ];
				}
				//simd : shift
				for(int s = 0;s < offset_str_row; s++){
					for(int k = 1;k < vector_len ;k++){
						str_row[k-1] = str_row[k];
					}
			}
			}else if(j > width-vector_len){ // desborde por derecha
				//simd : desplazamiento de puntero y levantar datos de memoria
				//levantamos el string con el puntero offseteado para no acceder afuera
				int offset_str_row = j - (width-vector_len);
				
				//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
				for(int k = 0;k < vector_len;k++){
					//Es lo mismo levantar siempre en este caso los ultimos vector_len caracteres
					str_row[vector_len - 1 - k] = seq1[j - vector_len - offset_str_row + k];
				}
				//simd : shift
				//shifteamos los chars para que quede bien (el sentido es contrario a cuando lo hagamos en simd)
				for(int s = 0;s < offset_str_row; s++){
				for(int k = vector_len-2;k >= 0;k--){
					str_row[k+1] = str_row[k];
				}
			}	

			}else{ //caso feliz
			//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
				for(int k = 0;k < vector_len;k++){
					str_row[vector_len - 1 - k] = seq1[j - vector_len + k];
				}
			}

			//vemos como los strings que levantamos se comparan en cada iteracion
			cerr<<"str_row = ";
			for(int i=0;i<vector_len;i++)
				cerr<<str_row[i];cerr<<endl;
			
			cerr<<"str_col = ";
			for(int i=0;i<vector_len;i++)
				cerr<<str_col[i];cerr<<endl;
			
			//left score
			//simd : leer de memoria (movdqu)
			for( int k = 0;k < vector_len;k++){
				left_score[vector_len - 1 - k] = score_matrix[offset_y + offset_x - vector_len + k];
			}
			
			//simd : padddw
			for( int k = 0;k < vector_len;k++){
				left_score[k] += constant_gap[k];
			}

			//up score
			//simd : copiar registro
			for( int k = 0;k < vector_len;k++){
				up_score[vector_len - 1 - k] = score_matrix[offset_y + offset_x - vector_len + k];
			}
			//simd : shift
			for(int k = vector_len-2;k >= 0;k--){
				up_score[k+1] = up_score[k];
			}
			//simd : insert
			up_score[0] = v_aux[j-1];

			//simd : padddw
			for( int k = 0;k < vector_len;k++){
				up_score[k] += constant_gap[k];
			}

			//diag score
			//simd : leer de memoria (movdqu)
			for( int k = 0;k < vector_len;k++){
				diag_score[vector_len - 1 - k] = score_matrix[offset_y + offset_x - 2*vector_len + k];
			}

			//simd : shift
			for(int k = vector_len-2;k >= 0;k--){
				diag_score[k+1] = diag_score[k];
			}

			//simd : insert
			diag_score[0] = v_aux[j-2];
		
			//simd : PUNPCKLBW
			for( int k = 0;k < vector_len;k++){
				cmp_missmatch[k] = str_row[k];
			}

			//simd : PUNPCKLBW
			for( int k = 0;k < vector_len;k++){
				cmp_match[k] = str_col[k];
			}

			//simd : compare
			for( int k = 0;k < vector_len;k++){
				cmp_match[k] = (cmp_match[k] == cmp_missmatch[k]);
			}

			//simd : pandn
			for( int k = 0;k < vector_len;k++){
				cmp_missmatch[k] = 1-cmp_match[k];
			}

			//simd : pmult por match
			for( int k = 0;k < vector_len;k++){
				cmp_match[k] *= alignment.parameters->match;
			}

			//simd : pmult por missmatch
			for( int k = 0;k < vector_len;k++){
				cmp_missmatch[k] *= alignment.parameters->missmatch;
			}

			//simd : padddw
			for( int k = 0;k < vector_len;k++){
				diag_score[k] += cmp_match[k] + cmp_missmatch[k];
			}
			
			cerr<<"Left_score : ";
			for( int k = 0;k < vector_len;k++){
				cerr << left_score[k] << " ";
			}cerr<<endl;

			cerr<<"up score : ";
			for( int k = 0;k < vector_len;k++){
				cerr << up_score[k] << " ";
			}cerr<<endl;
			cerr<<"diag score : ";
			for( int k = 0;k < vector_len;k++){
				cerr << diag_score[k] << " ";
			}cerr<<endl;

			//PMAXSW
			for( int k = 0;k < vector_len;k++){
				diag_score[k] = max(diag_score[k],up_score[k]);
			}
			
			//PMAXSW
			for( int k = 0;k < vector_len;k++){
				diag_score[k] = max(diag_score[k],left_score[k]);
			}

			//simd : mover a memoria
			for( int k = 0;k < vector_len;k++){
				score_matrix[offset_y + offset_x + k] = diag_score[vector_len - 1 - k];
			}

			//simd : PEXTRW
			if(j-vector_len >= 0)
				v_aux[j-vector_len] = diag_score[vector_len-1];
			
		}	
	}

	printScoreMatrix(score_matrix, &alignment, (int)vector_len);

	free(str_row);
	free(str_col);
	free(left_score);
	free(up_score);
	free(diag_score);
	free(constant_gap);
	free(cmp_missmatch);
	free(cmp_match);
	
	backtrack_solution_NW_C(
		score_matrix,
		alignment,
		vector_len,
		alignment.sequence_1->length-1,alignment.sequence_2->length-1,
		false,
		(score_fun_t)get_score_SSE,
		false
	);

	free(score_matrix);
}


void backtrack_solution_NW_C(
	short *score_matrix,
	Alignment& alignment,
	int vector_len,
	unsigned int x,unsigned int y,
	bool SW,
	score_fun_t score_fun,
	bool debug = false
){
	
	char* seq1 = alignment.sequence_1->sequence;
	char* seq2 = alignment.sequence_2->sequence;

	int seq1_len = alignment.sequence_1->length;
	int seq2_len = alignment.sequence_2->length;
	
	Sequence* best_seq_1 = alignment.result->sequence_1;
	Sequence* best_seq_2 = alignment.result->sequence_2;

	char* best_sequence_1 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	char* best_sequence_2 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	
	unsigned int length = 0;
	
	while(y != 0 || x != 0){

		if(y > 0 && x > 0){
			if(SW && score_fun(score_matrix, seq1_len,y,x,vector_len) == 0)
				break;

			short score_diag = score_fun(score_matrix, seq1_len,y-1,x-1,vector_len) + 
						  alignment.parameters->missmatch*(seq2[y] != seq1[x]) + 
						  alignment.parameters->match*(seq2[y] == seq1[x]);
			
			//short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = score_fun(score_matrix, seq1_len,y-1,x,vector_len) + alignment.parameters->gap;
			
			if(score_diag ==  score_fun(score_matrix, seq1_len,y,x,vector_len)){
				best_sequence_1[length] = seq2[y];
				best_sequence_2[length] = seq1[x];
				x--;
				y--;
			}else if(score_up == score_fun(score_matrix, seq1_len,y,x,vector_len)){
				best_sequence_1[length] = seq2[y];
				best_sequence_2[length] = '_';
				y--;
			}else{
				best_sequence_1[length] = '_';
				best_sequence_2[length] = seq1[x];
				x--;
			}
		}else if(x > 0){
			best_sequence_1[length] = '_';
			best_sequence_2[length] = seq2[x-1];
			x--;
		}else if(y > 0){
			best_sequence_1[length] = seq1[y-1];
			best_sequence_2[length] = '_';
			y--;
		}
		
		length++;
	}
	
	for(unsigned int i = 0;i < length/2;i++){
		char swap = best_sequence_1[i];
		best_sequence_1[i] = best_sequence_1[length-1-i];
		best_sequence_1[length-1-i] = swap;
		swap = best_sequence_2[i];
		best_sequence_2[i] = best_sequence_2[length-1-i];
		best_sequence_2[length-1-i] = swap;
	}
	
	cerr << "Best sequences" << endl;
	best_sequence_1[length]=0;
	best_sequence_2[length]=0;
	DBG(best_sequence_2);
	DBG(best_sequence_1);
	alignment.result->sequence_1 = new_Sequence_from_string(best_sequence_1);
	alignment.result->sequence_2 = new_Sequence_from_string(best_sequence_2);
	alignment.result->score = score_fun(score_matrix, seq1_len,alignment.sequence_2->length-1,alignment.sequence_1->length-1,vector_len);

	//creamos los nuevos strings y destruimos los anteriores
	best_seq_1 = new_Sequence_from_string(best_sequence_1);
	best_seq_2 = new_Sequence_from_string(best_sequence_2);

	destroy_Sequence(alignment.result->sequence_1);
	destroy_Sequence(alignment.result->sequence_2);
	alignment.result->sequence_1 = best_seq_1;
	alignment.result->sequence_2 = best_seq_2;

}

short get_score_SSE(short* score_matrix, int seq_row_len, int y,int x, int vector_len=4){
	int width = (1 + seq_row_len + vector_len - 1); //cantidad de diagonales por franja

	int franja = y / vector_len;	
	int offset_y = franja * width * vector_len;

	int diagonal = y % vector_len + x + 1;
	int cell = vector_len - 1 - y % vector_len;
	int offset_x = diagonal * vector_len + cell;

	return score_matrix[offset_y + offset_x];
}

short get_score_LIN(short* score_matrix, int seq_row_len, int y,int x, int vector_len=4){
	short** matrix_ptr = (short**) score_matrix;
	return matrix_ptr[y][x];
}