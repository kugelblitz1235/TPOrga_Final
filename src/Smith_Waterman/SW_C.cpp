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


void SW_C_withLogicSSE (Alignment& alignment){
	char* seq1 = alignment.sequence_1->sequence;	
	char* seq2 = alignment.sequence_2->sequence;
	unsigned int seq1_len = alignment.sequence_1->length;
	unsigned int seq2_len = alignment.sequence_2->length;
	
	
	//en este caso hardcodeamos el tamaño del vector
	int vector_len = 16;
	
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
					score_matrix[offset_y + offset_x + k] = 0;
				else
					score_matrix[offset_y + offset_x + k] = SHRT_MIN/2;
				count++;
			}			
		}
	}
	/******************************************************************************************************/
	short best_global = 0;
	unsigned int best_x=0;
	unsigned int best_y=0;

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

		if((i+1)*vector_len >= (int)seq2_len){
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
			for( int k = 0 ; k < vector_len ; k++){
				cmp_match[k] = str_col[k];
			}

			//simd : compare
			for( int k = 0 ; k < vector_len ; k++){
				cmp_match[k] = (cmp_match[k] == cmp_missmatch[k]);
			}

			//simd : pandn
			for( int k = 0 ; k < vector_len ; k++){
				cmp_missmatch[k] = 1-cmp_match[k];
			}

			//simd : pmult por match
			for( int k = 0 ; k < vector_len ; k++){
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
			
			//PMAXSW
			for( int k = 0;k < vector_len;k++){
				diag_score[k] = max(diag_score[k],up_score[k]);
			}
			
			//PMAXSW
			for( int k = 0;k < vector_len;k++){
				diag_score[k] = max(diag_score[k],left_score[k]);
			}
				
			//PMAXSW
			for( int k = 0;k < vector_len;k++){
				diag_score[k] = max(diag_score[k],(short)0);
			}

			int best_score = 0;
			if(best_global < best_score){
				best_global = best_score;
				best_y=y;
				best_x=x;
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

	free(str_row);
	free(str_col);
	free(left_score);
	free(up_score);
	free(diag_score);
	free(constant_gap);
	free(cmp_missmatch);
	free(cmp_match);
	
	backtracking_C(
		score_matrix,
		alignment,
		vector_len,
		best_x,best_y,
		true,
		(score_fun_t)get_score_SSE,
		false
	);

	free(score_matrix);
}

