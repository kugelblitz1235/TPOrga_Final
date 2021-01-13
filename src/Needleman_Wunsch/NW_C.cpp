#include <iostream>
#include <bitset>

#include "NW_C.hpp"

#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;

namespace NW{

	namespace C{

		namespace LIN{
			void NW(Alignment& alignment, bool debug){
				// Tamaños y punteros a las secuencias
				unsigned int seq1_len = alignment.sequence_1->length;
				unsigned int seq2_len = alignment.sequence_2->length;
				char* seq1 = alignment.sequence_1->sequence;
				char* seq2 = alignment.sequence_2->sequence;
				
				// Matriz de los puntajes que genera el algoritmo NW
				short** scores = (short**)malloc((seq2_len)*sizeof(short*));
				short* tmp_scores = (short*)malloc((seq2_len*seq1_len)*sizeof(short));
				
				// Inicializamos la matriz
				for(unsigned int y = 0;y < seq2_len;y++)
					scores[y] = &tmp_scores[y*seq1_len];
				
				scores[0][0] = 0;
				for(unsigned int y = 1;y < seq2_len;y++){
					scores[y][0] = scores[y-1][0] + alignment.parameters->gap;
				}
				for(unsigned int x = 1;x < seq1_len;x++){
					scores[0][x] = scores[0][x-1] + alignment.parameters->gap;
				}
				
				for(unsigned int y = 1;y < seq2_len;y++){
					for(unsigned int x = 1;x < seq1_len;x++){
						// Moverse por la diagonal, chequeando si se genera un matching
						short score_diag = scores[y-1][x-1] + 
									alignment.parameters->missmatch*(seq2[y]!= seq1[x]) + 
									alignment.parameters->match*(seq2[y] == seq1[x]);
						
						// Moverse por la izquierda sumando el puntaje de gap
						short score_left = scores[y][x-1] + alignment.parameters->gap;
						// Moverse por arriba sumando el puntaje de gap
						short score_up = scores[y-1][x] + alignment.parameters->gap;
						// Guardar en la matriz el maximo entre las 3 posibles direcciones
						short best_score = max(score_diag,score_left);
						best_score = max(best_score,score_up);
						scores[y][x] = best_score;
					}
				}
				
				if(debug){// Matriz utilizada para debuggear la matriz de puntajes
					alignment.matrix = (short *) scores;
				}
				
				// Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion mas inferior derecha
				backtracking_C(
					(short*)scores,
					alignment,
					1,
					seq1_len-1,seq2_len-1,
					false,
					(score_fun_t)get_score_LIN,
					false
				);

				
				if(!debug){
					free(tmp_scores);
					free(scores);
				}
			}
		}
		
		void with_logic_SSE (Alignment& alignment, int vector_len, bool debug){
			// tamaños y punteros a las secuencias
			char* seq1 = alignment.sequence_1->sequence;	
			char* seq2 = alignment.sequence_2->sequence;
			unsigned int seq1_len = alignment.sequence_1->length;
			unsigned int seq2_len = alignment.sequence_2->length;
			
			int height = ((seq2_len + vector_len - 1)/ vector_len); // Cantidad de "franjas" de diagonales
			int width = (1 + seq1_len + vector_len - 1); 			// Cantidad de diagonales por franja
			int score_matrix_sz = height * width * vector_len; 		// Cantidad de celdas de la matriz de puntajes

			// Matriz de puntajes				
			short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
			
			short* v_aux = (short*)malloc((width-1)*sizeof(short));
			// Llenar el vector auxiliar
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
							score_matrix[offset_y + offset_x + k] = i * vector_len * alignment.parameters->gap;
						else
							score_matrix[offset_y + offset_x + k] = SHRT_MIN/2;
						count++;
					}			
				}
			}
			/******************************************************************************************************/
			
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
					// Caso de desborde por abajo
					// Evita levantar de mas en la secuencia vertical
					// lee el tamanio del vector sin pasarse
					// y corrije shifteando
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
					//emulamos simd
					// String horizontal ------------------------------------------------------------------
					if(j-vector_len < 0){ 
						// desborde por izquierda
						// parte del vector cae fuera de la matriz
						// levantamos el string con el puntero offseteado para no acceder afuera
						// shifteamos para corregir.
						
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

					} else if(j > width-vector_len){ // Caso de desborde por derecha
						// levantamos el string con el puntero offseteado para no acceder afuera
						// shifteamos para corregir.

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
				
					// Comparacion de secuencias
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
			
			if(debug){
				alignment.matrix = score_matrix;
			}

			backtracking_C(
				score_matrix,
				alignment,
				vector_len,
				alignment.sequence_1->length-1,alignment.sequence_2->length-1,
				false,
				(score_fun_t)get_score_SSE,
				false
			);

			if(!debug)free(score_matrix);
		}

		namespace SSE{
			// Inicializar los valores del vector auxiliar y la matriz de puntajes
			void inicializar_casos_base(int width, int height, int vector_len, short* v_aux, short* score_matrix, Alignment& alignment){
				// Llenar vector auxiliar con un valor inicial negativo grande para no afectar los calculos
				for(int i = 0;i < width-1;i++){
					v_aux[i] = SHRT_MIN/2; 
				}

				// Inicializar por cada franja las primeras 2 diagonales. 
				// Se pone en cada posicion un valor negativo grande para no afectar los calculos posteriores
				__m128i diag;
				for(int i = 0 ; i < height ; i++){
					unsigned int offset_y = i * width * vector_len;
					// Broadcastear el valor SHRT_MIN/2 a nivel word en el registro diag
					diag = _mm_insert_epi16(diag,SHRT_MIN/2,0);
					diag = _mm_shufflelo_epi16(diag,0b0);
					diag = _mm_shuffle_epi32 (diag,0b0);
					_mm_storeu_si128((__m128i*)(score_matrix + offset_y), diag);
					// Insertar la penalidad adecuada para la franja en la posicion mas alta de la diagonal
					diag = _mm_insert_epi16(diag,i * vector_len * alignment.parameters->gap, 7); 			// 7 = vector_len - 1
					_mm_storeu_si128((__m128i*)(score_matrix + offset_y + vector_len), diag);
				}
			}
			
			// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
			__m128i leer_secuencia_columna(
				int i,
				int vector_len,
				int seq2_len,
				char* seq2,
				__m128i zeroes_xmm,
				__m128i reverse_mask_xmm){

				__m128i str_col_xmm;
				if((i+1)*vector_len >= (int)seq2_len){// Caso de desborde por abajo

					int offset_col = (i+1)*vector_len - seq2_len; 										// Indica cuanto hay que shiftear el string columna luego de levantarlo
					str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + seq2_len - vector_len) ); 

					__m128i shift_count_xmm = zeroes_xmm; 												// Se utiliza junto con offset_col para shiftear str_col_xmm
					__m128i shift_mask_xmm = zeroes_xmm;
					// Acomodar los caracteres en str_col_xmm para que queden en las posiciones correctas para la comparación mas adelante
					shift_count_xmm = _mm_insert_epi8(shift_count_xmm, offset_col * 8, 0);
					str_col_xmm = _mm_srl_epi64(str_col_xmm, shift_count_xmm); 							// str_col_xmm = |0...0|str_col|
					
					// shift_mask_xmm va a tener 1's donde haya caracteres inválidos
					shift_mask_xmm = _mm_cmpeq_epi8(shift_mask_xmm, shift_mask_xmm);
					shift_count_xmm = _mm_insert_epi8(shift_count_xmm, (char)(8 - offset_col) * 8, 0);
					shift_mask_xmm = _mm_sll_epi64(shift_mask_xmm, shift_count_xmm); 					// shift_mask_xmm = |1...1|0...0|
					// Combinar la mascara de caracteres invalidos con los caracteres validos
					str_col_xmm = _mm_or_si128(str_col_xmm, shift_mask_xmm); 							// str_col_xmm = |1...1|str_col|
				}else{ // Caso sin desborde
					// Levantar directamente de memoria, no hay desborde
					str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len) );
				}
				// Desempaquetar los caracteres en str_col_xmm para trabajar con words
				str_col_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
				// Invertir la secuencia de caracteres para compararlos correctamente mas adelante
				str_col_xmm = _mm_shuffle_epi8(str_col_xmm, reverse_mask_xmm);
				
				return str_col_xmm;
			}

			// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
			__m128i leer_secuencia_fila(
				int j,
				int vector_len,
				int width,
				char* seq1,
				__m128i zeroes_xmm
			) {
				__m128i str_row_xmm;
				__m128i shift_count = zeroes_xmm; // Se utiliza junto con offset_row para shiftear str_row_xmm en los casos de desborde
				
				if(j-vector_len < 0){ // Caso de desborde por izquierda
				
					int offset_str_row = vector_len - j; // Indica cuanto hay que shiftear a la izquierda el string fila luego de levantarlo
					str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1));
					
					// Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
					shift_count = _mm_insert_epi8(shift_count, offset_str_row * 8, 0);
					str_row_xmm = _mm_sll_epi64(str_row_xmm, shift_count); 				// str_row_xmm = |str_row|0...0|
						
				} else if(j > width-vector_len){ // Caso de desborde por derecha
					//Desplazamiento de puntero a derecha y levantar datos de memoria
					int offset_str_row = j - (width-vector_len); // Indica cuanto hay que shiftear a la derecha el string fila luego de levantarlo
					
					str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len - offset_str_row) ); 
					// Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
					shift_count = _mm_insert_epi8(shift_count, offset_str_row * 8, 0);
					str_row_xmm = _mm_srl_epi64(str_row_xmm, shift_count); 				// str_row_xmm = |0...0|str_row|
						
				}else{ // Caso sin desborde
					str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len) );
				}
				// Desempaquetar los caracteres en str_row_xmm para trabajar con words
				str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);
				return str_row_xmm;
			}

			//Calcula los puntajes resultantes de las comparaciones entre caracteres
			void calcular_scores(
				__m128i& left_score_xmm,
				__m128i& up_score_xmm,
				__m128i& diag_score_xmm,
				short* score_matrix,
				short* v_aux,
				int j,
				int offset_y,
				int offset_x,
				int vector_len,
				__m128i constant_gap_xmm,
				__m128i str_col_xmm,
				__m128i str_row_xmm,
				__m128i constant_missmatch_xmm,
				__m128i constant_match_xmm
				){

				// Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
				left_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - vector_len));
				left_score_xmm = _mm_add_epi16(left_score_xmm, constant_gap_xmm);
				
				// Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
				up_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - vector_len));
				up_score_xmm = _mm_srli_si128(up_score_xmm, 2);
				up_score_xmm = _mm_insert_epi16(up_score_xmm,v_aux[j-1],0b111);
				up_score_xmm = _mm_add_epi16(up_score_xmm, constant_gap_xmm);
				
				// Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
				// si coinciden o no los caracteres de la fila y columna correspondientes
				diag_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - 2*vector_len));
				diag_score_xmm = _mm_srli_si128(diag_score_xmm, 2);
				diag_score_xmm = _mm_insert_epi16(diag_score_xmm,v_aux[j-2],0b111);
				// Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
				__m128i cmp_match_xmm;
				cmp_match_xmm = _mm_cmpeq_epi16(str_col_xmm,str_row_xmm); 									// Mascara con unos en las posiciones donde coinciden los caracteres
				cmp_match_xmm = _mm_blendv_epi8(constant_missmatch_xmm,constant_match_xmm,cmp_match_xmm); 	// Seleccionar para cada posicion el puntaje correcto basado en la mascara previa
				diag_score_xmm = _mm_add_epi16(diag_score_xmm, cmp_match_xmm);
			}
		

			void NW (Alignment& alignment, bool debug){
				// Tamaños y punteros a las secuencias
				char* seq1 = alignment.sequence_1->sequence;	
				char* seq2 = alignment.sequence_2->sequence;
				unsigned int seq1_len = alignment.sequence_1->length;
				unsigned int seq2_len = alignment.sequence_2->length;
				
				// Equivalentes a registros nombrados:
				__m128i constant_gap_xmm, constant_missmatch_xmm, constant_match_xmm, zeroes_xmm;
				__m128i str_row_xmm, str_col_xmm, left_score_xmm, up_score_xmm, diag_score_xmm;
				__m128i reverse_mask_xmm;
				
				// Broadcastear el valor de gap, a nivel word, en el registro
				constant_gap_xmm = _mm_insert_epi16(constant_gap_xmm,alignment.parameters->gap,0);
				constant_gap_xmm = _mm_shufflelo_epi16(constant_gap_xmm,0b0);
				constant_gap_xmm = _mm_shuffle_epi32 (constant_gap_xmm,0b0);
				// Broadcastear el valor de missmatch, a nivel word, en el registro
				constant_missmatch_xmm = _mm_insert_epi16(constant_missmatch_xmm,alignment.parameters->missmatch,0);
				constant_missmatch_xmm = _mm_shufflelo_epi16(constant_missmatch_xmm,0b0);
				constant_missmatch_xmm = _mm_shuffle_epi32 (constant_missmatch_xmm,0b0);
				// Broadcastear el valor de match, a nivel word, en el registro
				constant_match_xmm = _mm_insert_epi16(constant_match_xmm,alignment.parameters->match,0);
				constant_match_xmm = _mm_shufflelo_epi16(constant_match_xmm,0b0);
				constant_match_xmm = _mm_shuffle_epi32 (constant_match_xmm,0b0);
				// Máscara de ceros
				zeroes_xmm = _mm_setzero_si128();
				
				// Máscara utilizada para invertir el string almacenado en un registro
				char reverse_mask[16] = {0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1};
				reverse_mask_xmm = _mm_loadu_si128((__m128i*)reverse_mask);
				
				// El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
				int vector_len = 8;
				
				int height = ((seq2_len + vector_len - 1)/ vector_len); 			// Cantidad de "franjas" de diagonales
				int width = (1 + seq1_len + vector_len - 1); 						// Cantidad de diagonales por franja
				int score_matrix_sz = height * width * vector_len; 					// Cantidad de celdas de la matriz
				
				// Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores
				short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
				short* v_aux = (short*)malloc((width-1)*sizeof(short));
					
				inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, alignment);
				
				/******************************************************************************************************/

				for( int i = 0 ; i < height ; i++){
					int offset_y = i * width * vector_len;

					str_col_xmm = leer_secuencia_columna(i, vector_len, seq2_len, seq2, zeroes_xmm, reverse_mask_xmm);

					for( int j = 2; j < width ; j++){
						int offset_x = j * vector_len;
						// String horizontal ------------------------------------------------------------------
						str_row_xmm = leer_secuencia_fila(j, vector_len, width, seq1, zeroes_xmm);
						
						// Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
						calcular_scores(
							left_score_xmm,
							up_score_xmm,
							diag_score_xmm,
							score_matrix,
							v_aux,
							j,
							offset_y,
							offset_x,
							vector_len,
							constant_gap_xmm,
							str_col_xmm,
							str_row_xmm,
							constant_missmatch_xmm,
							constant_match_xmm
						);
						
						// Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
						diag_score_xmm = _mm_max_epi16(diag_score_xmm,up_score_xmm);
						diag_score_xmm = _mm_max_epi16(diag_score_xmm,left_score_xmm);

						// Almacenamos el puntaje máximo en la posición correcta de la matriz
						_mm_storeu_si128((__m128i*)(score_matrix + offset_y + offset_x), diag_score_xmm);

						if(j>=vector_len){
							v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_xmm, 0b0000);
						}

					}	
				}

				if(debug){// Utilizar para debuggear los valores en la matriz de puntajes
					alignment.matrix = score_matrix;
				}

				// Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion mas inferior derecha
				backtracking_C(
					score_matrix,
					alignment,
					vector_len,
					alignment.sequence_1->length-1,alignment.sequence_2->length-1,
					false,
					(score_fun_t)get_score_SSE,
					false
				);

				if(!debug) free(score_matrix);
			}
		}
		
		namespace AVX{
			// Inicializar los valores del vector auxiliar y la matriz de puntajes
			void inicializar_casos_base(int width, int height, int vector_len, short* v_aux, short* score_matrix, Alignment& alignment){
				// Llenar vector auxiliar con un valor inicial negativo grande para no afectar los calculos
				for(int i = 0;i < width-1;i++){
					v_aux[i] = SHRT_MIN/2;
				}

				__m128i temp_xmm;
				__m256i diag_ymm;
				// Inicializar por cada franja las primeras 2 diagonales. 
				// Se pone en cada posicion un valor negativo grande para no afectar los calculos posteriores
				for(int i = 0 ; i < height ; i++){
					unsigned int offset_y = i * width * vector_len;
					// Broadcastear el valor SHRT_MIN/2 a nivel word en el registro diag
					temp_xmm = _mm_insert_epi16(temp_xmm,SHRT_MIN/2,0);
					diag_ymm = _mm256_broadcastw_epi16(temp_xmm);
					_mm256_storeu_si256((__m256i*)(score_matrix + offset_y), diag_ymm);
					// Insertar la penalidad adecuada para la franja en la posicion mas alta de la diagonal
					diag_ymm = _mm256_insert_epi16(diag_ymm,i * vector_len * alignment.parameters->gap, 15); // 15 = vector_len - 1
					_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + vector_len), diag_ymm);
				}
			}

			// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
			__m256i leer_secuencia_columna(
				int i,
				int vector_len,
				int seq2_len,
				char* seq2,
				__m256i zeroes_ymm,
				__m128i reverse_mask_xmm,
				__m128i shift_mask_right_xmm){

				__m128i str_col_xmm;
				__m256i str_col_ymm;

				if((i+1)*vector_len >= (int)seq2_len){// Caso de desborde por abajo
					int offset_str_col = (i+1)*vector_len - seq2_len;										// Indica cuanto hay que shiftear el string columna luego de levantarlo
																								
					str_col_xmm = _mm_loadu_si128((__m128i*)(seq2 + seq2_len - vector_len));				// str_col_xmm = |0|str_col|	

					// Poner el valor offset_str_col en todos los bytes de offset_str_col_xmm, 
					// el cual indica cuantos bytes hay que shiftear los caracteres para dejarlos en la posición correcta
					__m128i offset_str_col_xmm = _mm_insert_epi8(offset_str_col_xmm,offset_str_col,0);
					offset_str_col_xmm = _mm_broadcastb_epi8(offset_str_col_xmm);						// offset_str_col_xmm = |offset|...|offset|
					// Sumar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
					offset_str_col_xmm = _mm_add_epi8(shift_mask_right_xmm,offset_str_col_xmm);			// offset_str_col_xmm = |0x7F + offset|0x7E + offset|...|0x70 + offset|
					// Shiftear utilizando la mascara offset_str_col_xmm 
					str_col_xmm = _mm_shuffle_epi8(str_col_xmm,offset_str_col_xmm);
					// Todos los elementos que sean basura van a convertirse en el valor 0xFFFF, haciendo que nunca matcheen mas adelante
					str_col_xmm = _mm_blendv_epi8(str_col_xmm,offset_str_col_xmm, offset_str_col_xmm);	// str_col_xmm = |1...1|str_col|
				}else{ // Caso sin desborde
					str_col_xmm = _mm_loadu_si128((__m128i*)(seq2 + i * vector_len));
				}

				// Invertir el string almacenado en str_col_xmm
				str_col_xmm = _mm_shuffle_epi8 (str_col_xmm, reverse_mask_xmm);
				__m128i zeroes_xmm = _mm_setzero_si128();
				// Desempaquetar los caracteres almacenados en str_col_xmm para trabajar con words
				__m128i str_col_hi_xmm = _mm_unpackhi_epi8(str_col_xmm, zeroes_xmm);
				__m128i str_col_lo_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
				str_col_ymm = _mm256_set_m128i(str_col_hi_xmm,str_col_lo_xmm);
				
				return str_col_ymm;
			}

			// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
			__m256i leer_secuencia_fila(
				int j,
				int vector_len,
				int width,
				char* seq1,
				__m256i zeroes_ymm,
				__m128i shift_mask_right_xmm,
				__m128i shift_mask_left_xmm
			) {
				__m128i str_row_xmm;
				__m256i str_row_ymm;

				if(j-vector_len < 0){ // Caso de desborde por izquierda
					int offset_str_row = vector_len - j;													// Indica cuanto hay que shiftear a la izquierda el string fila luego de levantarlo
					str_row_xmm = _mm_loadu_si128((__m128i*)(seq1) );										
					
					// Broadcastear el offset_str_row en offset_str_row_xmm 
					__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,offset_str_row,0);
					offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);							// offset_str_col_xmm = |offset|...|offset|
					// Sumar con shift_mask_left_xmm para indicar cuanto hay que shiftear a la izquierda los caracteres para que queden en las posiciones correctas
					offset_str_row_xmm = _mm_sub_epi8(shift_mask_left_xmm,offset_str_row_xmm);				// offset_str_col_xmm = |0xF - offset|0xE - offset|...|0x0 - offset|
					
					// Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
					str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);							// str_row_xmm = |str_row|0...0|
						
				}else if(j > width-vector_len){ // Caso de desborde por derecha
					//Desplazamiento de puntero a derecha y levantar datos de memoria
					int offset_str_row = j - (width-vector_len);											// Indica cuanto hay que shiftear a la derecha el string fila luego de levantarlo
					str_row_xmm = _mm_loadu_si128((__m128i*)(seq1 + j - vector_len - offset_str_row) );
					
					// Broadcastear el offset_str_row en offset_str_row_xmm 
					__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,offset_str_row,0);
					offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);							// offset_str_col_xmm = |offset|...|offset|
					// Sumar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
					offset_str_row_xmm = _mm_add_epi8(shift_mask_right_xmm,offset_str_row_xmm);				// offset_str_col_xmm = |0xF + offset|0xE + offset|...|0x0 + offset|
					// Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
					str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);							// str_row_xmm = |0...0|str_row|
				
				}else{ // Caso sin desborde
					str_row_xmm = _mm_loadu_si128((__m128i*)(seq1 + j - vector_len));
				}
				
				// Desempaquetar los caracteres en str_row_xmm para trabajar con words
				__m128i zeroes_xmm = _mm_setzero_si128();
				__m128i str_row_hi_xmm = _mm_unpackhi_epi8(str_row_xmm, zeroes_xmm);
				__m128i str_row_lo_xmm = _mm_unpacklo_epi8(str_row_xmm, zeroes_xmm);
				str_row_ymm = _mm256_set_m128i(str_row_hi_xmm,str_row_lo_xmm);
				return str_row_ymm;
			}
			
			//Calcula los puntajes resultantes de las comparaciones entre caracteres
			void calcular_scores(
				__m256i& left_score_ymm,
				__m256i& up_score_ymm,
				__m256i& diag_score_ymm,
				short* score_matrix,
				short* v_aux,
				int j,
				int vector_len,
				__m256i constant_gap_ymm,
				__m256i str_col_ymm,
				__m256i str_row_ymm,
				__m256i constant_missmatch_ymm,
				__m256i constant_match_ymm,
				__m256i diag1_ymm,
				__m256i diag2_ymm
				){

					

				// Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
				left_score_ymm = diag2_ymm;
				left_score_ymm = _mm256_add_epi16 (left_score_ymm, constant_gap_ymm);
				
				// Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
				up_score_ymm = diag2_ymm;
				// El shift de a bytes se hace en cada linea de 128 bits,
				// por lo que hay que mover manualmente el word en la posicion mas baja en la linea de 128 bits
				// mas alta, e insertarlo en la posicion mas alta de la linea de 128 bits mas baja
				short medium_score = _mm256_extract_epi16( up_score_ymm, 0b1000);					// |0xF...0x8|0x7...0x0| -> 0x8
				up_score_ymm = _mm256_bsrli_epi128(up_score_ymm, 2); 								// |0xF...0x8|0x7...0x0| -> |0x0...0x9|0x0...0x1|
				up_score_ymm = _mm256_insert_epi16(up_score_ymm,medium_score,0b0111); 				// |0x0...0x9|0x0...0x1| -> |0x0...0x9|0x8...0x1|
				up_score_ymm = _mm256_insert_epi16(up_score_ymm,v_aux[j-1],15); 					// 15 = vector_len - 1
				up_score_ymm = _mm256_add_epi16(up_score_ymm, constant_gap_ymm);
				
				// Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
				// si coinciden o no los caracteres de la fila y columna correspondientes
				diag_score_ymm = diag1_ymm;
				medium_score = _mm256_extract_epi16( diag_score_ymm, 0b1000);
				diag_score_ymm = _mm256_bsrli_epi128(diag_score_ymm, 2);
				diag_score_ymm = _mm256_insert_epi16(diag_score_ymm,medium_score,0b0111);
				
				diag_score_ymm = _mm256_insert_epi16(diag_score_ymm,v_aux[j-2],15); 				// 15 = vector_len - 1
				
				// Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
				__m256i cmp_match_ymm = _mm256_cmpeq_epi16(str_col_ymm,str_row_ymm);
				cmp_match_ymm = _mm256_blendv_epi8(constant_missmatch_ymm,constant_match_ymm,cmp_match_ymm); 
				diag_score_ymm = _mm256_add_epi16(diag_score_ymm, cmp_match_ymm);

			}

			void NW (Alignment& alignment, bool debug){
				// Tamaños y punteros a las secuencias
				char* seq1 = alignment.sequence_1->sequence;	
				char* seq2 = alignment.sequence_2->sequence;
				unsigned int seq1_len = alignment.sequence_1->length;
				unsigned int seq2_len = alignment.sequence_2->length;
				
				// Equivalentes a registros nombrados:
				__m256i constant_gap_ymm, constant_missmatch_ymm, constant_match_ymm, zeroes_ymm;
				__m256i str_row_ymm, str_col_ymm, left_score_ymm, up_score_ymm, diag_score_ymm;
				__m128i reverse_mask_xmm, shift_mask_right_xmm, shift_mask_left_xmm;
				__m256i diag1_ymm, diag2_ymm;

				__m128i temp_xmm;
				// Broadcastear el valor de gap, a nivel word, en el registro
				temp_xmm = _mm_insert_epi16(temp_xmm,alignment.parameters->gap,0);
				constant_gap_ymm = _mm256_broadcastw_epi16(temp_xmm);
				// Broadcastear el valor de missmatch, a nivel word, en el registro
				temp_xmm = _mm_insert_epi16(temp_xmm,alignment.parameters->missmatch,0);
				constant_missmatch_ymm = _mm256_broadcastw_epi16(temp_xmm);
				// Broadcastear el valor de match, a nivel word, en el registro
				temp_xmm = _mm_insert_epi16(temp_xmm,alignment.parameters->match,0);
				constant_match_ymm = _mm256_broadcastw_epi16(temp_xmm);
				// Máscara de ceros
				zeroes_ymm = _mm256_setzero_si256();

				// Cada posicion en la máscara tiene 0x70 sumado, para que luego de la suma el bit más significativo esté activado para las posiciones que caen fuera de los límites
				char shift_mask_right[16] = {0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F};
				char shift_mask_left[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF};
				char reverse_mask[16] = {0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0};
				// Máscara utilizada para invertir el string almacenado en un registro
				reverse_mask_xmm = _mm_loadu_si128 ((__m128i*)reverse_mask);
				// Máscara utilizada para shiftear a derecha los caracteres 
				shift_mask_right_xmm =  _mm_loadu_si128((__m128i*)shift_mask_right);
				// Máscara utilizada para shiftear a izquierda los caracteres
				shift_mask_left_xmm =  _mm_loadu_si128((__m128i*)shift_mask_left);
				
				// El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
				int vector_len = 16;
				
				int height = ((seq2_len + vector_len - 1)/ vector_len); // Cantidad de "franjas" de diagonales
				int width = (1 + seq1_len + vector_len - 1); 			// Cantidad de diagonales por franja
				int score_matrix_sz = height * width * vector_len; 		// Cantidad de celdas de la matriz

				// Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores	
				short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
				short* v_aux = (short*)malloc((width-1)*sizeof(short));
				
				inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, alignment);
				/******************************************************************************************************/
				
				for( int i = 0 ; i < height ; i++){
					int offset_y = i * width * vector_len;
					// String vertical --------------------------------------------------------------------
					str_col_ymm = leer_secuencia_columna(i, vector_len, seq2_len, seq2, zeroes_ymm, reverse_mask_xmm, shift_mask_right_xmm);
					// Cargar las primeras 2 diagonales de la franja actual necesarios para el primer calculo dentro de la franja
					diag1_ymm = _mm256_loadu_si256((__m256i const*) (score_matrix + offset_y));
					diag2_ymm = _mm256_loadu_si256((__m256i const*) (score_matrix + offset_y + vector_len));
					for( int j = 2; j < width ; j++){
						int offset_x = j * vector_len;
						// String horizontal ------------------------------------------------------------------
						str_row_ymm = leer_secuencia_fila(j, vector_len, width, seq1, zeroes_ymm, shift_mask_right_xmm, shift_mask_left_xmm);
						// Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
						calcular_scores(
							left_score_ymm,
							up_score_ymm,
							diag_score_ymm,
							score_matrix,
							v_aux,
							j,
							vector_len,
							constant_gap_ymm,
							str_col_ymm,
							str_row_ymm,
							constant_missmatch_ymm,
							constant_match_ymm,
							diag1_ymm,
							diag2_ymm
						);
						// Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
						diag_score_ymm = _mm256_max_epi16(diag_score_ymm,up_score_ymm);
						diag_score_ymm = _mm256_max_epi16(diag_score_ymm,left_score_ymm);

						// Almacenamos el puntaje máximo en la posición correcta de la matriz
						_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + offset_x), diag_score_ymm);

						if(j>=vector_len){
							v_aux[j - vector_len] =  _mm256_extract_epi16 (diag_score_ymm, 0x0);
						}

						// Actualizar las 2 diagonales anteriores para la siguiente iteracion, actualizando diag2 con el valor de la actual
						diag1_ymm = diag2_ymm;
						diag2_ymm = diag_score_ymm;
					}	
				}
				if(debug){// Utilizar para debuggear los valores en la matriz de puntajes
					alignment.matrix = score_matrix;
				}
				
				//Ejecutamos backtracking partiendo de la posición inferior derecha de la matriz de puntajes
				backtracking_C(
					score_matrix,
					alignment,
					vector_len,
					alignment.sequence_1->length-1,alignment.sequence_2->length-1,
					false,
					(score_fun_t)get_score_SSE,
					false
				);

				if(!debug) free(score_matrix);
			}
		}
			
		namespace AVX512{
			// Directiva utilizada para emular la representación de un zmm con sus parte ymm y xmm
			// |          zmm          |
			// |    256    |    ymm    |
			// | 128 | 128 | 128 | xmm |
			union SIMDreg{
				__m128i x;
				__m256i y;
				__m512i z;
			};

			// Registros globales utilizados
			SIMDreg constant_gap_mm, constant_missmatch_mm, constant_match_mm, zeroes_mm, ones_mm;
			SIMDreg str_row_mm, str_col_mm, left_score_mm, up_score_mm, diag_score_mm;
			SIMDreg str_reverse_mask_mm, str_shift_right_mask_mm, str_shift_left_mask_mm;
			SIMDreg str_512_unpacklo_epi8_mask_mm;
			SIMDreg score_512_rot_right_word_mask_mm;
			SIMDreg diag1_mm, diag2_mm;
			
			//tamaños y punteros a las secuencias
			char *seq1, *seq2;
			unsigned int seq1_len, seq2_len;
			
			short str_reverse_mask[32] = {
				0x1F,0x1E,0x1D,0x1C,0x1B,0x1A,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0
			};
			uint64_t str_512_unpacklo_epi8_mask[8] = {
				0x0,0xFF,0x1,0xFF,0x2,0xFF,0x3,0xFF
			};
			// Mascara para rotar a derecha a nivel word un zmm
			short score_512_rot_right_word_mask[32] = {
				0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x0
			};
			
			// El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
			const int vector_len = 32;
			
			int height; 				//cantidad de "franjas" de diagonales
			int width; 					//cantidad de diagonales por franja

			short* score_matrix;
			short* v_aux;

			// Inicializamos los valores del vector auxiliar y la matriz de puntajes
			void inicializar_casos_base(Alignment& alignment){
				// Llenar vector auxiliar con un valor inicial negativo grande para no afectar los calculos
				for(int i = 0;i < width-1;i++){
					v_aux[i] = SHRT_MIN/2;
				}

				// Inicializar por cada franja las primeras 2 diagonales. 
				// Se pone en cada posicion un valor negativo grande para no afectar los calculos posteriores
				for(int i = 0 ; i < height ; i++){
					unsigned int offset_y = i * width * vector_len;
					SIMDreg temp_mm;
					SIMDreg diag_mm;
					// Broadcastear el valor SHRT_MIN/2 a nivel word en el registro diag_zmm
					diag_mm.z = _mm512_maskz_set1_epi16(0xFFFFFFFF, SHRT_MIN/2);
					_mm512_storeu_si512((__m512i*)(score_matrix + offset_y), diag_mm.z);
					temp_mm.x = diag_mm.x;
					// Insertar la penalidad adecuada para la franja en la posicion mas alta de la diagonal
					temp_mm.x = _mm_insert_epi16(temp_mm.x, i * vector_len * alignment.parameters->gap, 7);
					diag_mm.z = _mm512_inserti64x2(diag_mm.z, temp_mm.x, 3);
					_mm512_storeu_si512((__m512i*)(score_matrix + offset_y + vector_len), diag_mm.z);
				}
			}

			// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
			void leer_secuencia_columna(int i){

				if((i+1)*vector_len >= (int)seq2_len){// Caso de desborde por abajo
					int offset_str_col = (i+1)*vector_len - seq2_len; // Indica cuanto hay que shiftear a la derecha el shift_right_mask

					// Shiftear a derecha la cantidad de posiciones equivalente a caracteres invalidos,
					// para evitar levantar memoria invalida
					__mmask32 shift_right_mask = 0xFFFFFFFF;
					shift_right_mask = (shift_right_mask >> offset_str_col); 
					// Seleccionar los caracteres validos a derecha con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
					// A su vez poner en los lugares de posiciones invalidas todos 0s, para evitar que coincida con algun caracter de la secuencia columna
					str_col_mm.y = _mm256_mask_loadu_epi8(ones_mm.y, shift_right_mask, (__m256i*)(seq2 + seq2_len - vector_len + offset_str_col));  	// str_row_mm = |0...0|str_row|
				}else{ // Caso sin desborde
					str_col_mm.y = _mm256_loadu_si256((__m256i*)(seq2 + i * vector_len));
				}
				// Desempaquetar el string en str_col_ymm para trabajar a nivel words
				str_col_mm.z = _mm512_permutexvar_epi64(str_512_unpacklo_epi8_mask_mm.z, str_col_mm.z);
				str_col_mm.z = _mm512_unpacklo_epi8(str_col_mm.z, zeroes_mm.z); 
				
				// Invertir el string en str_col_zmm
				str_col_mm.z = _mm512_permutexvar_epi16 (str_reverse_mask_mm.z, str_col_mm.z);
			}
	
			// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
			void leer_secuencia_fila(int j) {
				if(j-vector_len < 0){ // Caso de desborde por izquierda
					int offset_str_row = vector_len - j;
					// Shiftear a izquierda la cantidad de posiciones equivalente a caracteres invalidos,
					// para evitar levantar memoria invalida
					__mmask32 shift_left_mask = 0xFFFFFFFF;
					shift_left_mask <<= offset_str_row; 
					// Seleccionar los caracteres validos a derecha con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
					// A su vez poner en los lugares de posiciones invalidas todos 0s, para evitar que coincida con algun caracter de la secuencia columna
					str_row_mm.y = _mm256_maskz_loadu_epi8(shift_left_mask, (__m256i*)(seq1 - offset_str_row)); 	// str_row_mm = |str_row|0...0|

				}else if(j > width-vector_len){ // Caso de desborde por derecha
					int offset_str_row = j - (width-vector_len);
					// Shiftear a derecha la cantidad de posiciones equivalente a caracteres invalidos,
					// para evitar levantar memoria invalida
					__mmask32 shift_right_mask = 0xFFFFFFFF;
					shift_right_mask >>= offset_str_row; 
					// Seleccionar los caracteres validos a derecha con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
					// A su vez poner en los lugares de posiciones invalidas todos 0s, para evitar que coincida con algun caracter de la secuencia columna
					str_row_mm.y = _mm256_maskz_loadu_epi8(shift_right_mask, (__m256i*)(seq1 + j - vector_len)); 	// str_row_mm = |0...0|str_row|

				}else{ // Caso sin desborde
					str_row_mm.y = _mm256_loadu_si256((__m256i*)(seq1 + j - vector_len));
				}
				
				// Desempaquetamos los caracteres en str_row_ymm para trabajr a nivel word
				str_row_mm.z = _mm512_permutexvar_epi64(str_512_unpacklo_epi8_mask_mm.z, str_row_mm.z);
				str_row_mm.z = _mm512_unpacklo_epi8(str_row_mm.z, zeroes_mm.z); 
			}

			//Calcula los puntajes resultantes de las comparaciones entre caracteres
			void calcular_scores(int j){
				// Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
				left_score_mm.z = diag2_mm.z;
				left_score_mm.z = _mm512_add_epi16 (left_score_mm.z, constant_gap_mm.z);
				
				// Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
				up_score_mm.z = diag2_mm.z;
				up_score_mm.x = _mm_insert_epi16(up_score_mm.x,v_aux[j-1],0b0);
				up_score_mm.z = _mm512_permutexvar_epi16(score_512_rot_right_word_mask_mm.z, up_score_mm.z);
				up_score_mm.z = _mm512_add_epi16(up_score_mm.z, constant_gap_mm.z);
				
				// Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
				// si coinciden o no los caracteres de la fila y columna correspondientes
				diag_score_mm.z = diag1_mm.z;
				diag_score_mm.x = _mm_insert_epi16(diag_score_mm.x,v_aux[j-2],0b0);
				diag_score_mm.z = _mm512_permutexvar_epi16 (score_512_rot_right_word_mask_mm.z, diag_score_mm.z);
				
				// Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
				SIMDreg cmp_match_mm;
				__mmask32 cmp_mask = _mm512_cmpeq_epi16_mask(str_col_mm.z,str_row_mm.z);
				cmp_match_mm.z = _mm512_mask_blend_epi16(cmp_mask, constant_missmatch_mm.z,constant_match_mm.z); 
				diag_score_mm.z = _mm512_add_epi16(diag_score_mm.z, cmp_match_mm.z);
			}


			void NW (Alignment& alignment, bool debug){
				// Tamaños y punteros de las secuencias
				seq1 = alignment.sequence_1->sequence;	
				seq2 = alignment.sequence_2->sequence;
				seq1_len = alignment.sequence_1->length;
				seq2_len = alignment.sequence_2->length;
				
				// Broadcastear el valor de gap, a nivel word, en el registro
				constant_gap_mm.z = _mm512_maskz_set1_epi16(0xFFFFFFFF, alignment.parameters->gap);
				// Broadcastear el valor de missmatch, a nivel word, en el registro
				constant_missmatch_mm.z = _mm512_maskz_set1_epi16(0xFFFFFFFF, alignment.parameters->missmatch);
				// Broadcastear el valor de match, a nivel word, en el registro
				constant_match_mm.z = _mm512_maskz_set1_epi16(0xFFFFFFFF, alignment.parameters->match);
				// Máscara de ceros
				zeroes_mm.z = _mm512_setzero_si512();
				// Máscara de unos
				ones_mm.z = _mm512_maskz_set1_epi8(0xFFFFFFFFFFFFFFFFULL, (uint8_t)0xFF);
				
				// Máscara utilizada para invertir el string almacenado en un registro
				str_reverse_mask_mm.z = _mm512_loadu_si512 ((__m512i*)str_reverse_mask);
				// Máscara utilizada para desempaquetar los caracteres a nivel word 
				str_512_unpacklo_epi8_mask_mm.z =  _mm512_loadu_si512((__m512i*)str_512_unpacklo_epi8_mask);
				// Máscara utilizada para rotar el string a la derecha
				score_512_rot_right_word_mask_mm.z =  _mm512_loadu_si512((__m512i*)score_512_rot_right_word_mask);
				
				height = ((seq2_len + vector_len - 1)/ vector_len);					//cantidad de "franjas" de diagonales
				width = (1 + seq1_len + vector_len - 1); 							//cantidad de diagonales por franja
				int score_matrix_sz = height * width * vector_len; 					// largo de diagonal
				
				// Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores	
				score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
				v_aux = (short*)malloc((width-1)*sizeof(short));

				inicializar_casos_base(alignment);
				/******************************************************************************************************/

				for( int i = 0 ; i < height ; i++){
					int offset_y = i * width * vector_len;

					// String vertical --------------------------------------------------------------------
					leer_secuencia_columna(i);
					// Cargar las primeras 2 diagonales de la franja actual necesarios para el primer calculo dentro de la franja
					diag1_mm.z = _mm512_loadu_si512((__m512i const*) (score_matrix + offset_y));
					diag2_mm.z = _mm512_loadu_si512((__m512i const*) (score_matrix + offset_y + vector_len));
					
					for( int j = 2; j < width ; j++){
						int offset_x = j * vector_len;
						// String horizontal ------------------------------------------------------------------
						leer_secuencia_fila(j);
						
						//Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
						calcular_scores(j);
						
						// Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
						diag_score_mm.z = _mm512_max_epi16(diag_score_mm.z,up_score_mm.z);
						diag_score_mm.z = _mm512_max_epi16(diag_score_mm.z,left_score_mm.z);

						// Almacenamos el puntaje máximo en la posición correcta de la matriz
						_mm512_storeu_si512((__m512i*)(score_matrix + offset_y + offset_x), diag_score_mm.z);

						if(j>=vector_len){
							v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_mm.x, 0x0);
						}
						// Actualizar las 2 diagonales anteriores para la siguiente iteracion, actualizando diag2 con el valor de la actual
						diag1_mm.z = diag2_mm.z;
						diag2_mm.z = diag_score_mm.z;
					}	
				}
				if(debug){	// Utilizar para debuggear los valores en la matriz de puntajes
					alignment.matrix = score_matrix;
				}
							
				//Ejecutamos backtracking partiendo de la posición inferior derecha de la matriz de puntajes
				backtracking_C(
					score_matrix,
					alignment,
					vector_len,
					alignment.sequence_1->length-1,alignment.sequence_2->length-1,
					false,
					(score_fun_t)get_score_SSE,
					false
				);

				if(!debug) free(score_matrix);
			}
		}
	}

	Alignment* get_alignment(std::string implementation, char* sequence_1, char* sequence_2, short gap, short missmatch, short match){
		// Creo la estructura vacia
		Alignment* alignment = new_alignment();
		
		// Relleno la estructura con los valores correspondientes
		alignment->sequence_1 = new_Sequence_from_string(sequence_1);
		alignment->sequence_2 = new_Sequence_from_string(sequence_2);
		(alignment->parameters)->match = match;
		(alignment->parameters)->missmatch = missmatch;
		(alignment->parameters)->gap = gap;

		if(implementation.compare("C_LIN") == 0) NW::C::LIN::NW(*alignment, true);
		else if(implementation.compare("C_SSE") == 0)NW::C::SSE::NW(*alignment, true);
		else if(implementation.compare("C_AVX512") == 0)NW::C::AVX512::NW(*alignment, true);
		else if(implementation.compare("C_AVX") == 0)NW::C::AVX::NW(*alignment, true);
		else if(implementation.compare("ASM_LIN") == 0)NW_ASM_LIN(alignment);
		else if(implementation.compare("ASM_SSE") == 0)NW_ASM_SSE(alignment, true);
		else if(implementation.compare("ASM_AVX") == 0)NW_ASM_AVX(alignment, true);
		else if(implementation.compare("ASM_AVX512") == 0)NW_ASM_AVX512(alignment, true);
		else throw "No existe la implementación ingresada.";
		

		//devuelvo la estructura modificada
		return alignment;
	}
}
