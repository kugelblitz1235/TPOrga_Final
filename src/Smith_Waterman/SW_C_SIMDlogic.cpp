#include "SW_C_SIMDlogic.hpp"

namespace SW{
namespace C{
namespace SIMDlogic{
    void SW (Alignment& alignment, int vector_len, bool debug){
		// tamaños y punteros a las secuencias
		char* seq1 = alignment.sequence_1->sequence;	
		char* seq2 = alignment.sequence_2->sequence;
		unsigned int seq1_len = alignment.sequence_1->length;
		unsigned int seq2_len = alignment.sequence_2->length;
		
		int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
		int width = (1 + seq1_len + vector_len - 1); 			//cantidad de diagonales por franja
		int score_matrix_sz = height * width * vector_len; 		// Cantidad de celdas de la matriz de puntajes

		// Matriz de puntajes	
		short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
		
		// Llenar el vector auxiliar
		short* v_aux = (short*)malloc((width-1)*sizeof(short));
		for(int i = 0;i < width-1;i++){
			v_aux[i] = SHRT_MIN/2;
		}
		
		// Inicializamos la matriz 
		int count=0;
		for(int i = 0 ; i < height ; i++){
			unsigned int offset_y = i * width * vector_len;
			for( int j = 0; j < 2 ; j++){
				unsigned int offset_x = j * vector_len;
				// Emular simd
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
		short best_global = 0;								// Máximo puntaje obtenido
		unsigned int best_x = 0;							// Coordenada x del máximo puntaje
		unsigned int best_y = 0;							// Coordenada y del máximo puntaje

		// Arrays auxiliares para el calculo de los scores
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
			// Levantar strings -------------------------------------------------------------------
			// String vertical --------------------------------------------------------------------
			if((i+1)*vector_len >= (int)seq2_len){
				int offset_col = (i+1)*vector_len - seq2_len;
				// Caso de desborde por abajo
				// Evita levantar de mas en la secuencia vertical
				// lee el tamanio del vector sin pasarse
				// y corrije shifteando
				for( int k = 0;k < vector_len;k++){
					str_col[vector_len - 1 - k] = seq2[i * vector_len - offset_col + k];
				}

				//simd : shift
				for(int s = 0;s < offset_col; s++){
					for(int k = vector_len-2;k >= 0;k--){
						str_col[k+1] = str_col[k];
					}
				}
				//simd : or para evitar basura
				for(int s = 0; s < offset_col; s++){
					str_col[s] = 1;
				}

			}else{
				//simd : leer de memoria (movdqu)
				//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
				for( int k = 0;k < vector_len;k++){
					str_col[vector_len - 1 - k] = seq2[i * vector_len + k];
				}
			}

			// hace reverse de la secuencia vertical
			// para poder comparar directamente con la horizontal
			//simd : shuffle
			for (int i = 0 ; i < vector_len / 2; i++){
				char temp = str_col[i];
				str_col[i] = str_col[vector_len - i - 1];
				str_col[vector_len - i - 1] = temp;
			}

			for( int j = 2; j < width ; j++){
				// Procesamiento de las diagonales dentro de la franja
				int offset_x = j * vector_len;
				// Emular simd
				// String horizontal ------------------------------------------------------------------
				if(j-vector_len < 0){ // Caso de desborde por izquierda
					// parte del vector cae fuera de la matriz
					// levantamos el string con el puntero offseteado para no acceder afuera
					// shifteamos para corregir.
					//simd : desplazamiento de puntero y levantar datos de memoria
					
					int offset_str_row = vector_len - j;
					//simd : leer de memoria (movdqu)
					//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
					for(int k = 0;k < vector_len;k++){
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
					
					//simd: agrega 0 para evitar basura
					for(int s = 0; s < offset_str_row; s++){
						str_row[vector_len-1-s] = 0;
					}
				}else if(j > width-vector_len){ // Caso de desborde por derecha
					// levantamos el string con el puntero offseteado para no acceder afuera
					// shifteamos para corregir.
					//simd : desplazamiento de puntero y levantar datos de memoria
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

					//simd: agrega 0 para evitar basura
					for(int s = 0; s < offset_str_row; s++){
						str_row[s] = 0;
					}

				}else{ // Caso sin desborde
				//aclaracion hay que levantar la cantidad de caracteres que nos indica vector_len, no más
					for(int k = 0;k < vector_len;k++){
						str_row[vector_len - 1 - k] = seq1[j - vector_len + k];
					}
				}

				// Calculo scores --------------------------------------------------------------------
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
				int index = 0;
				for( int k = 0;k < vector_len;k++){
					if(best_score < diag_score[k]){
						best_score = diag_score[k];
						index=vector_len-1-k;
					}
				}

				// Actualizamos la posición máxima junto con sus posiciones
				if(best_global < best_score){
					best_global = best_score;
					best_y = vector_len * i + (vector_len-1) - index;
					best_x = j - vector_len + index;
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
	
		if(debug){// Utilizar para debuggear la matriz de puntajes obtenida
			alignment.matrix = score_matrix;

		}

		// Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion del puntaje máximo obtenido
		backtracking_C(
			score_matrix,
			alignment,
			vector_len,
			best_x,best_y,
			true,
			(score_fun_t)get_score_SSE,
			true
		);

		if(!debug)free(score_matrix);
	}
}
}
}
