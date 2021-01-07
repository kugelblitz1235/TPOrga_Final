#include <iostream>
#include "NW_C.hpp"

#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;

namespace NW{

namespace LIN{
	void NW_C_LIN(Alignment& alignment, bool debug){
		
		unsigned int seq1_len = alignment.sequence_1->length;
		unsigned int seq2_len = alignment.sequence_2->length;
		char* seq1 = alignment.sequence_1->sequence;
		char* seq2 = alignment.sequence_2->sequence;
		
			
		short** scores = (short**)malloc((seq2_len)*sizeof(short*));
		short* tmp_scores = (short*)malloc((seq2_len*seq1_len)*sizeof(short));
		
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
				
				short score_diag = scores[y-1][x-1] + 
							alignment.parameters->missmatch*(seq2[y]!= seq1[x]) + 
							alignment.parameters->match*(seq2[y] == seq1[x]);
				
				short score_left = scores[y][x-1] + alignment.parameters->gap;
				
				short score_up = scores[y-1][x] + alignment.parameters->gap;
				
				short best_score = max(score_diag,score_left);
				
				best_score = max(best_score,score_up);
				scores[y][x] = best_score;
			}
		}
		
		if(debug){
			alignment.matrix = (short *) scores;
			// cerr << "Score Matrix" << endl;
			// for(unsigned int y = 0;y < seq2_len;y++){
			// 	for(unsigned int x = 0;x < seq1_len;x++){
			// 		cerr << (int)scores[y][x] << " ";
			// 	}
			// 	cerr << endl << endl;
			// }	
		}
		
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
namespace SSE{

	void NW_C_withLogicSSE (Alignment& alignment, bool debug){
		char* seq1 = alignment.sequence_1->sequence;	
		char* seq2 = alignment.sequence_2->sequence;
		unsigned int seq1_len = alignment.sequence_1->length;
		unsigned int seq2_len = alignment.sequence_2->length;
		
		
		//en este caso hardcodeamos el tamaño del vector
		int vector_len = 8;
		
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

			// Levantar strings -------------------------------------------------------------------
			// String vertical --------------------------------------------------------------------
			if((i+1)*vector_len >= (int)seq2_len){
				// Desborde por abajo
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
					
					//simd: agrega 0 para evitar basura
					for(int s = 0; s < offset_str_row; s++){
						str_row[vector_len-1-s] = 0;
					}

				} else if(j > width-vector_len){ 
					// desborde por derecha
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

				}else{ //caso feliz
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

	void inicializar_casos_base(int width, int height, int vector_len, short* v_aux, short* score_matrix, Alignment& alignment){
		// llenamos el vector auxiliar
		for(int i = 0;i < width-1;i++){
			v_aux[i] = SHRT_MIN/2;
		}

		// inicializar casos base en matriz
		for(int i = 0 ; i < height ; i++){
			unsigned int offset_y = i * width * vector_len;
			__m128i diag;
			diag = _mm_insert_epi16(diag,SHRT_MIN/2,0);
			diag = _mm_shufflelo_epi16(diag,0b0);
			diag = _mm_shuffle_epi32 (diag,0b0);
			_mm_storeu_si128((__m128i*)(score_matrix + offset_y), diag);
			diag = _mm_insert_epi16(diag,i * vector_len * alignment.parameters->gap, 7); //7 = vector_len - 1
			_mm_storeu_si128((__m128i*)(score_matrix + offset_y + vector_len), diag);

			//printf("Width: %d\nHeight: %d\noffset_y: %d\n", width, height, offset_y);
		}
	}
	
	__m128i leer_secuencia_columna(
		int i,
		int vector_len,
		int seq2_len,
		char* seq2,
		__m128i zeroes_xmm,
		__m128i reverse_mask_xmm){

		__m128i str_col_xmm;
		if((i+1)*vector_len >= (int)seq2_len){
			// Desborde por abajo
			// Evita levantar de mas en la secuencia vertical
			// lee el tamanio del vector sin pasarse
			// y corrije shifteando
			int offset_col = (i+1)*vector_len - seq2_len;
		
			//simd : leer de memoria (movdqu)
			str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + seq2_len - vector_len) );

			__m128i shift_count = zeroes_xmm;
			__m128i shift_mask;
			shift_count = _mm_insert_epi8(shift_count, offset_col*8, 0);
			str_col_xmm = _mm_srl_epi64(str_col_xmm, shift_count);

			shift_mask = _mm_cmpeq_epi8(shift_mask, shift_mask);
			shift_count = _mm_insert_epi8(shift_count, (char)(8-offset_col)*8, 0);
			shift_mask = _mm_sll_epi64(shift_mask, shift_count);

			str_col_xmm = _mm_or_si128(str_col_xmm, shift_mask);

			
		}else{
			//simd : leer de memoria (movdqu)
			str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len) );
			
		}
		str_col_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
		// Reverse de la secuencia vertical
		str_col_xmm = _mm_shuffle_epi8(str_col_xmm, reverse_mask_xmm);
		
		return str_col_xmm;
	}


	__m128i leer_secuencia_fila(
		int j,
		int vector_len,
		int width,
		char* seq1,
		__m128i zeroes_xmm
	) {
		__m128i str_row_xmm;
		
		if(j-vector_len < 0){ //desborde por izquierda
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = vector_len - j;
				//simd : leer de memoria (movdqu)
				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1));
				
				__m128i shift_count = zeroes_xmm;
				shift_count = _mm_insert_epi8(shift_count, offset_str_row*8, 0);
				str_row_xmm = _mm_sll_epi64(str_row_xmm, shift_count);
				
		} else if(j > width-vector_len){ // desborde por derecha
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = j - (width-vector_len);

				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len - offset_str_row) );
				__m128i shift_count = zeroes_xmm;
				shift_count = _mm_insert_epi8(shift_count, offset_str_row*8, 0);
				str_row_xmm = _mm_srl_epi64(str_row_xmm, shift_count);
				
		}else{ //caso feliz
				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len) );
		}
		
		str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);
		return str_row_xmm;
	}

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

			
		//left score
		left_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - vector_len));
		left_score_xmm = _mm_add_epi16(left_score_xmm, constant_gap_xmm);
		
		//up score
		up_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - vector_len));
		up_score_xmm = _mm_srli_si128(up_score_xmm, 2);
		up_score_xmm = _mm_insert_epi16(up_score_xmm,v_aux[j-1],0b111);
		up_score_xmm = _mm_add_epi16(up_score_xmm, constant_gap_xmm);
		
		//diag score
		diag_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - 2*vector_len));
		diag_score_xmm = _mm_srli_si128(diag_score_xmm, 2);
		diag_score_xmm = _mm_insert_epi16(diag_score_xmm,v_aux[j-2],0b111);
		
		//compare the 2 strings and put the right penalty (match or missmatch) on each position
		__m128i cmp_match_xmm;
		cmp_match_xmm = _mm_cmpeq_epi16(str_col_xmm,str_row_xmm);
		cmp_match_xmm = _mm_blendv_epi8(constant_missmatch_xmm,constant_match_xmm,cmp_match_xmm); 
		diag_score_xmm = _mm_add_epi16(diag_score_xmm, cmp_match_xmm);
	}


	void NW_C_SSE (Alignment& alignment, bool debug){
		char* seq1 = alignment.sequence_1->sequence;	
		char* seq2 = alignment.sequence_2->sequence;
		unsigned int seq1_len = alignment.sequence_1->length;
		unsigned int seq2_len = alignment.sequence_2->length;
		
		// Equivalentes a registros nombrados:
		__m128i constant_gap_xmm, constant_missmatch_xmm, constant_match_xmm, zeroes_xmm;
		__m128i str_row_xmm, str_col_xmm, left_score_xmm, up_score_xmm, diag_score_xmm;
		__m128i reverse_mask_xmm;
		
		constant_gap_xmm = _mm_insert_epi16(constant_gap_xmm,alignment.parameters->gap,0);
		constant_gap_xmm = _mm_shufflelo_epi16(constant_gap_xmm,0b0);
		constant_gap_xmm = _mm_shuffle_epi32 (constant_gap_xmm,0b0);
		constant_missmatch_xmm = _mm_insert_epi16(constant_missmatch_xmm,alignment.parameters->missmatch,0);
		constant_missmatch_xmm = _mm_shufflelo_epi16(constant_missmatch_xmm,0b0);
		constant_missmatch_xmm = _mm_shuffle_epi32 (constant_missmatch_xmm,0b0);
		constant_match_xmm = _mm_insert_epi16(constant_match_xmm,alignment.parameters->match,0);
		constant_match_xmm = _mm_shufflelo_epi16(constant_match_xmm,0b0);
		constant_match_xmm = _mm_shuffle_epi32 (constant_match_xmm,0b0);
		zeroes_xmm = _mm_setzero_si128();

		//each element has a 0x70 added, so after addition the most significative bit is activated for the trash characters
		//|0001|0000|0003|0002|0005|0004|0007|0006|0009|0008|000B|000A|000D|000C|000F|000E|
		char reverse_mask[16] = {0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1};
		
		reverse_mask_xmm = _mm_loadu_si128((__m128i*)reverse_mask);
		
		//en este caso hardcodeamos el tamaño del vector
		int vector_len = 8;
		
		int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
		int width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
		int score_matrix_sz = height * width * vector_len; // largo de diagonal
			
		short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
		
		short* v_aux = (short*)malloc((width-1)*sizeof(short));
			
		inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, alignment);
		
		/******************************************************************************************************/

		for( int i = 0 ; i < height ; i++){
			int offset_y = i * width * vector_len;

			// Levantar strings -------------------------------------------------------------------
			// String vertical --------------------------------------------------------------------
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
				
				diag_score_xmm = _mm_max_epi16(diag_score_xmm,up_score_xmm);
				diag_score_xmm = _mm_max_epi16(diag_score_xmm,left_score_xmm);

				//save the max score in the right position of score matrix
				_mm_storeu_si128((__m128i*)(score_matrix + offset_y + offset_x), diag_score_xmm);

				if(j>=vector_len){
					v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_xmm, 0b0000);
				}

			}	
		}
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

		if(!debug) free(score_matrix);
	}
}
namespace AVX{
	void inicializar_casos_base(int width, int height, int vector_len, short* v_aux, short* score_matrix, Alignment& alignment){
			// llenamos el vector auxiliar
			for(int i = 0;i < width-1;i++){
				v_aux[i] = SHRT_MIN/2;
			}
		
			// inicializar casos base en matriz
			for(int i = 0 ; i < height ; i++){
				unsigned int offset_y = i * width * vector_len;
				__m128i temp_xmm;
				__m256i diag_ymm;
				temp_xmm = _mm_insert_epi16(temp_xmm,SHRT_MIN/2,0);
				diag_ymm = _mm256_broadcastw_epi16(temp_xmm);
				_mm256_storeu_si256((__m256i*)(score_matrix + offset_y), diag_ymm);
				diag_ymm = _mm256_insert_epi16(diag_ymm,i * vector_len * alignment.parameters->gap, 7);
				_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + vector_len), diag_ymm);
			}
		}
		
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

			if((i+1)*vector_len >= (int)seq2_len){
				// Desborde por abajo
				// Evita levantar de mas en la secuencia vertical
				// lee el tamanio del vector sin pasarse
				// y corrije shifteando
				int offset_str_col = (i+1)*vector_len - seq2_len;
			
				//simd : leer de memoria (movdqu)
				str_col_xmm = _mm_load_si128((__m128i*)(seq2 + seq2_len - vector_len));

				__m128i offset_str_col_xmm = _mm_insert_epi8(offset_str_col_xmm,offset_str_col,0);
				offset_str_col_xmm = _mm_broadcastb_epi8(offset_str_col_xmm);
				offset_str_col_xmm = _mm_add_epi8(shift_mask_right_xmm,offset_str_col_xmm);
				
				str_col_xmm = _mm_shuffle_epi8(str_col_xmm,offset_str_col_xmm);
				//todos los elementos que sean basura van a convertirse en el valor 0xFFFF, haciendo que nunca matcheen mas adelante ni de casualidad
				str_col_xmm = _mm_blendv_epi8(str_col_xmm,offset_str_col_xmm, offset_str_col_xmm);
			}else{
				//simd : leer de memoria (movdqu)
				str_col_xmm = _mm_load_si128((__m128i*)(seq2 + i * vector_len));
			}
			str_col_xmm = _mm_shuffle_epi8 (str_col_xmm, reverse_mask_xmm);
			str_col_ymm = _mm256_set_m128i(str_col_xmm,str_col_xmm);
			str_col_ymm = _mm256_unpacklo_epi8 (str_col_ymm, zeroes_ymm);
			// Reverse de la secuencia vertical
			
			return str_col_ymm;
		}


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
			
			if(j-vector_len < 0){ //desborde por izquierda
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = vector_len - j;
				//simd : leer de memoria (movdqu)
				str_row_xmm = _mm_load_si128((__m128i*)(seq1 + j - vector_len + offset_str_row) );

				__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,offset_str_row,0);
				offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
				offset_str_row_xmm = _mm_sub_epi8(shift_mask_left_xmm,offset_str_row_xmm);
		
				str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
			
			}else if(j > width-vector_len){ // desborde por derecha
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = j - (width-vector_len);
				
				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len - offset_str_row) );

				__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,offset_str_row,0);
				offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
				offset_str_row_xmm = _mm_add_epi8(shift_mask_right_xmm,offset_str_row_xmm);
				str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
			
			}else{ //caso feliz
				str_row_xmm = _mm_load_si128((__m128i*)(seq1 + j - vector_len));
			}
			str_row_ymm = _mm256_set_m128i(str_row_xmm,str_row_xmm);
			str_row_ymm = _mm256_unpacklo_epi8(str_row_ymm,zeroes_ymm);
			return str_row_ymm;
		}

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

				
			//left score
			left_score_ymm = diag2_ymm;
			left_score_ymm = _mm256_add_epi16 (left_score_ymm, constant_gap_ymm);
			
			//up score
			up_score_ymm = diag2_ymm;
			up_score_ymm = _mm256_bsrli_epi128(up_score_ymm, 2);
			up_score_ymm = _mm256_insert_epi16(up_score_ymm,v_aux[j-1],7);
			up_score_ymm = _mm256_add_epi16(up_score_ymm, constant_gap_ymm);
			
			//diag score
			diag_score_ymm = diag1_ymm;
			diag_score_ymm = _mm256_bsrli_epi128(diag_score_ymm, 2);
			diag_score_ymm = _mm256_insert_epi16(diag_score_ymm,v_aux[j-2],7);
			
			//compare the 2 strings and put the right penalty (match or missmatch) on each position
			__m256i cmp_match_ymm = _mm256_cmpeq_epi16(str_col_ymm,str_row_ymm);
			cmp_match_ymm = _mm256_blendv_epi8(constant_missmatch_ymm,constant_match_ymm,cmp_match_ymm); 
			diag_score_ymm = _mm256_add_epi16(diag_score_ymm, cmp_match_ymm);
		}


		void NW_C_AVX (Alignment& alignment, bool debug){
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
			temp_xmm = _mm_insert_epi16(temp_xmm,alignment.parameters->gap,0);
			constant_gap_ymm = _mm256_broadcastw_epi16(temp_xmm);
			temp_xmm = _mm_insert_epi16(temp_xmm,alignment.parameters->missmatch,0);
			constant_missmatch_ymm = _mm256_broadcastw_epi16(temp_xmm);
			temp_xmm = _mm_insert_epi16(temp_xmm,alignment.parameters->match,0);
			constant_match_ymm = _mm256_broadcastw_epi16(temp_xmm);
			zeroes_ymm = _mm256_setzero_si256();

			//each element has a 0x70 added, so after addition the most significative bit is activated for the trash characters
			char shift_mask_right[16] = {0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F};
			char shift_mask_left[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF};
			char reverse_mask[16] = {0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0};


			reverse_mask_xmm = _mm_loadu_si128 ((__m128i*)reverse_mask);
			shift_mask_right_xmm =  _mm_loadu_si128((__m128i*)shift_mask_right);
			shift_mask_left_xmm =  _mm_loadu_si128((__m128i*)shift_mask_left);
			
			int vector_len = 16;
			
			int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
			int width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
			int score_matrix_sz = height * width * vector_len; // largo de diagonal
				
			short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
			
			short* v_aux = (short*)malloc((width-1)*sizeof(short));
				
			inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, alignment);
			
			/******************************************************************************************************/
			
			for( int i = 0 ; i < height ; i++){
				int offset_y = i * width * vector_len;

				// Levantar strings -------------------------------------------------------------------
				// String vertical --------------------------------------------------------------------
				str_col_ymm = leer_secuencia_columna(i, vector_len, seq2_len, seq2, zeroes_ymm, reverse_mask_xmm, shift_mask_right_xmm);

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
					
					diag_score_ymm = _mm256_max_epi16(diag_score_ymm,up_score_ymm);
					diag_score_ymm = _mm256_max_epi16(diag_score_ymm,left_score_ymm);

					//save the max score in the right position of score matrix
					_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + offset_x), diag_score_ymm);

					if(j>=vector_len){
						v_aux[j - vector_len] =  _mm256_extract_epi16 (diag_score_ymm, 0x0);
					}

					diag1_ymm = diag2_ymm;
					diag2_ymm = diag_score_ymm;
				}	
			}
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

			if(!debug) free(score_matrix);
		}
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

		if(implementation.compare("C_LIN") == 0) LIN::NW_C_LIN(*alignment, true);
		else if(implementation.compare("C_SSE") == 0)SSE::NW_C_SSE(*alignment, true);
		else if(implementation.compare("C_SSE") == 0)AVX::NW_C_AVX(*alignment, true);
		else if(implementation.compare("ASM_LIN") == 0)NW_ASM_LIN(alignment);
		else if(implementation.compare("ASM_SSE") == 0)NW_ASM_SSE(alignment, true);
		//else if(implementation.compare("ASM_AVX") == 0)NW_ASM_AVX(alignment, true);
		else throw "No existe la implementación ingresada.";
		

		//devuelvo la estructura modificada
		return alignment;
	}
}