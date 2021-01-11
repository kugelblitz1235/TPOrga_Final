#include <iostream>
#include <bitset>
#include "SW_C.hpp"

#define DBG(x) cerr << #x << " = " << (x) <<"\n"


using namespace std;

namespace SW {

void SW_C_LIN(
	Alignment& alignment,
	bool debug
){
	unsigned int seq1_len = alignment.sequence_1->length;
	unsigned int seq2_len = alignment.sequence_2->length;
	char* seq1 = alignment.sequence_1->sequence;
	char* seq2 = alignment.sequence_2->sequence;
	
	short** scores = (short**)malloc((seq2_len)*sizeof(short*));
	short* tmp_scores = (short*)malloc((seq2_len*seq1_len)*sizeof(short));
	
	for(unsigned int y = 0;y < seq2_len;y++)
		scores[y] = &tmp_scores[y*seq1_len];
	
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
	
	/*if(debug){
		cerr << "Score Matrix" << endl;
		for(unsigned int y = 0;y < seq2_len;y++){
			for(unsigned int x = 0;x < seq1_len;x++){
				cerr << (int)scores[y][x] << " ";
			}
			cerr << endl << endl;
		}	
	}*/
	

	if(debug){
		alignment.matrix = (short *)scores;
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

	if(!debug){
		free(tmp_scores);
		free(scores);
	}
}


	void SW_C_withLogicSSE (Alignment& alignment, int vector_len, bool debug){
		char* seq1 = alignment.sequence_1->sequence;	
		char* seq2 = alignment.sequence_2->sequence;
		unsigned int seq1_len = alignment.sequence_1->length;
		unsigned int seq2_len = alignment.sequence_2->length;
		
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
		unsigned int best_x = 0;
		unsigned int best_y = 0;

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
					
					//simd: agrega 0 para evitar basura
					for(int s = 0; s < offset_str_row; s++){
						str_row[vector_len-1-s] = 0;
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
				int index;
				for( int k = 0;k < vector_len;k++){
					if(best_score < diag_score[k]){
						best_score = diag_score[k];
						index=vector_len-1-k;
					}
				}

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
		
		// DBG(best_global);

		if(debug){
			alignment.matrix = score_matrix;

		}


		backtracking_C(
			score_matrix,
			alignment,
			vector_len,
			best_x,best_y,
			true,
			(score_fun_t)get_score_SSE,
			false
		);

		if(!debug)free(score_matrix);
	}
	namespace SSE{

		void inicializar_casos_base(int width, int height, int vector_len, short* v_aux, short* score_matrix, __m128i zeroes_xmm) {
			//llenamos el vector auxiliar
			for(int i = 0;i < width-1;i++){
				v_aux[i] = SHRT_MIN/2;
			}

			for(int i = 0 ; i < height ; i++){
				unsigned int offset_y = i * width * vector_len;
				_mm_storeu_si128((__m128i*)(score_matrix + offset_y), zeroes_xmm);
				_mm_storeu_si128((__m128i*)(score_matrix + offset_y + vector_len), zeroes_xmm);
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
				str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len - offset_col) );

				__m128i shift_count = zeroes_xmm;
				__m128i shift_mask = zeroes_xmm;
				shift_count = _mm_insert_epi8(shift_count, offset_col*8, 0);
				str_col_xmm = _mm_srl_epi64(str_col_xmm, shift_count);

				shift_mask = _mm_cmpeq_epi8(shift_mask, shift_mask);
				shift_count = _mm_insert_epi8(shift_count, (char)(8-offset_col)*8, 0);
				shift_mask = _mm_sll_epi64(shift_mask, shift_count);

				str_col_xmm = _mm_or_si128(str_col_xmm, shift_mask);

				// str_col_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
				
			}else{
				//simd : leer de memoria (movdqu)
				str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len) );
				// str_col_xmm = _mm_unpacklo_epi8(str_col_xmm,zeroes_xmm);
				
			}
			str_col_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
			// Reverse de la secuencia vertical
			str_col_xmm = _mm_shuffle_epi8(str_col_xmm,reverse_mask_xmm);
			
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
					str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len + offset_str_row) );
					
					__m128i shift_count = zeroes_xmm;
					shift_count = _mm_insert_epi8(shift_count, offset_str_row*8, 0);
					str_row_xmm = _mm_sll_epi64(str_row_xmm, shift_count);
					// str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);
					
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

		void actualizar_posicion_maxima(
			int &best_global,
			int &best_x,
			int &best_y,
			int vector_len,
			int i,
			int j,
			__m128i diag_score_xmm
		){
			//find the index of the maximum word in the 128bit register
			__m128i nums_xmm =  diag_score_xmm;
			__m128i nums_copy_xmm = nums_xmm;
			__m128i nums_s_xmm;
			
			nums_s_xmm = _mm_srli_si128(nums_xmm,1*2);	
			nums_xmm = _mm_max_epi16(nums_xmm,nums_s_xmm);	
			nums_s_xmm = _mm_srli_si128 (nums_xmm,2*2);
			nums_xmm = _mm_max_epi16(nums_xmm,nums_s_xmm);
			nums_s_xmm = _mm_srli_si128 (nums_xmm,4*2);
			nums_xmm = _mm_max_epi16(nums_xmm,nums_s_xmm);
			
			nums_xmm = _mm_shufflelo_epi16(nums_xmm,0b0);
			nums_xmm = _mm_shuffle_epi32 (nums_xmm,0b0);

			__m128i index_xmm = _mm_cmpeq_epi16(nums_xmm,nums_copy_xmm);
			index_xmm = _mm_packs_epi16(index_xmm,index_xmm);
			int64_t index_mask = _mm_extract_epi64(index_xmm,0);

			int max_index = (__builtin_ffsll(index_mask)-1)/8;
			short max_local_score =  _mm_extract_epi16 (nums_xmm, 0b0000);
			if(best_global < max_local_score){
				
				best_global = max_local_score;
				best_y = vector_len * i + (vector_len-1) - max_index;
				best_x = j - vector_len + max_index;
			}
		}

	}

	void SW_C_SSE(Alignment& alignment, bool debug){
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
		
		SSE::inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, zeroes_xmm);

	/******************************************************************************************************/

		int best_global = 0;
		int best_y = 0;
		int best_x = 0;

		for( int i = 0 ; i < height ; i++){
			int offset_y = i * width * vector_len;

			// Levantar strings -------------------------------------------------------------------
			// String vertical --------------------------------------------------------------------
			str_col_xmm = SSE::leer_secuencia_columna(i, vector_len, seq2_len, seq2, zeroes_xmm, reverse_mask_xmm);

			for( int j = 2; j < width ; j++){
				int offset_x = j * vector_len;
				//emulamos simd
				str_row_xmm = SSE::leer_secuencia_fila(j, vector_len, width, seq1, zeroes_xmm);
				
				// Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
				SSE::calcular_scores(
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
				diag_score_xmm = _mm_max_epi16(diag_score_xmm,zeroes_xmm);

				//save the max score in the right position of score matrix
				_mm_storeu_si128((__m128i*)(score_matrix + offset_y + offset_x), diag_score_xmm);
				
				if(j>=vector_len){
					v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_xmm, 0b0000);
				}
				
				SSE::actualizar_posicion_maxima(best_global,best_x,best_y,vector_len,i,j,diag_score_xmm);
			}	
		}

		// printScoreMatrix2(score_matrix,&alignment,vector_len);
		if(debug){
			alignment.matrix = score_matrix;
		}

		// printf("%s\n",seq1);
		// printf("%s\n",seq2);
		
		backtracking_C(
			score_matrix,
			alignment,
			vector_len,
			best_x,best_y,
			true,
			(score_fun_t)get_score_SSE,
			false
		);

		if(!debug) free(score_matrix);
	}

	
	namespace AVX{
		void inicializar_casos_base(
			int width, 
			int height, 
			int vector_len, 
			short* v_aux, 
			short* score_matrix,
			__m256i zeroes_ymm, 
			Alignment& alignment){
			// llenamos el vector auxiliar
			for(int i = 0;i < width-1;i++){
				v_aux[i] = SHRT_MIN/2;
			}

			// inicializar casos base en matriz
			for(int i = 0 ; i < height ; i++){
				unsigned int offset_y = i * width * vector_len;
				_mm256_storeu_si256((__m256i*)(score_matrix + offset_y), zeroes_ymm);
				_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + vector_len), zeroes_ymm);
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
				str_col_xmm = _mm_loadu_si128((__m128i*)(seq2 + seq2_len - vector_len));

				__m128i offset_str_col_xmm = _mm_insert_epi8(offset_str_col_xmm,offset_str_col,0);
				offset_str_col_xmm = _mm_broadcastb_epi8(offset_str_col_xmm);
				offset_str_col_xmm = _mm_add_epi8(shift_mask_right_xmm,offset_str_col_xmm);
				
				str_col_xmm = _mm_shuffle_epi8(str_col_xmm,offset_str_col_xmm);
				//todos los elementos que sean basura van a convertirse en el valor 0xFFFF, haciendo que nunca matcheen mas adelante ni de casualidad
				str_col_xmm = _mm_blendv_epi8(str_col_xmm,offset_str_col_xmm, offset_str_col_xmm);
			}else{
				//simd : leer de memoria (movdqu)
				str_col_xmm = _mm_loadu_si128((__m128i*)(seq2 + i * vector_len));
			}
			str_col_xmm = _mm_shuffle_epi8 (str_col_xmm, reverse_mask_xmm);
			
			__m128i zeroes_xmm = _mm_setzero_si128();
			__m128i str_col_hi_xmm = _mm_unpackhi_epi8(str_col_xmm, zeroes_xmm);
			__m128i str_col_lo_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
			str_col_ymm = _mm256_set_m128i(str_col_hi_xmm,str_col_lo_xmm);
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
				str_row_xmm = _mm_loadu_si128((__m128i*)(seq1) );
				
				__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,offset_str_row,0);
				offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
				offset_str_row_xmm = _mm_sub_epi8(shift_mask_left_xmm,offset_str_row_xmm);
				
				str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
					
			}else if(j > width-vector_len){ // desborde por derecha
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = j - (width-vector_len);
				str_row_xmm = _mm_loadu_si128((__m128i*)(seq1 + j - vector_len - offset_str_row) );
				
				__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,offset_str_row,0);
				offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
				offset_str_row_xmm = _mm_add_epi8(shift_mask_right_xmm,offset_str_row_xmm);
				str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
			
			}else{ //caso feliz
				str_row_xmm = _mm_loadu_si128((__m128i*)(seq1 + j - vector_len));
			}
			__m128i zeroes_xmm = _mm_setzero_si128();
			__m128i str_row_hi_xmm = _mm_unpackhi_epi8(str_row_xmm, zeroes_xmm);
			__m128i str_row_lo_xmm = _mm_unpacklo_epi8(str_row_xmm, zeroes_xmm);
			str_row_ymm = _mm256_set_m128i(str_row_hi_xmm,str_row_lo_xmm);
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
			short medium_score = _mm256_extract_epi16( up_score_ymm, 0b1000);
			up_score_ymm = _mm256_bsrli_epi128(up_score_ymm, 2); // TODO
			up_score_ymm = _mm256_insert_epi16(up_score_ymm,medium_score,0b0111);

			up_score_ymm = _mm256_insert_epi16(up_score_ymm,v_aux[j-1],15);
			up_score_ymm = _mm256_add_epi16(up_score_ymm, constant_gap_ymm);
			
			//diag score
			diag_score_ymm = diag1_ymm;
			medium_score = _mm256_extract_epi16( diag_score_ymm, 0b1000);
			diag_score_ymm = _mm256_bsrli_epi128(diag_score_ymm, 2); // TODO
			diag_score_ymm = _mm256_insert_epi16(diag_score_ymm,medium_score,0b0111);
			
			diag_score_ymm = _mm256_insert_epi16(diag_score_ymm,v_aux[j-2],15);
			
			//compare the 2 strings and put the right penalty (match or missmatch) on each position
			__m256i cmp_match_ymm = _mm256_cmpeq_epi16(str_col_ymm,str_row_ymm);
			cmp_match_ymm = _mm256_blendv_epi8(constant_missmatch_ymm,constant_match_ymm,cmp_match_ymm); 
			diag_score_ymm = _mm256_add_epi16(diag_score_ymm, cmp_match_ymm);

		}

		
		void actualizar_posicion_maxima(
			int &best_global,
			int &best_x,
			int &best_y,
			int vector_len,
			int i,
			int j,
			__m256i diag_score_ymm
		){
			//find the index of the maximum word in the 128bit register
			__m256i nums_ymm =  diag_score_ymm;
			__m256i nums_copy_ymm = nums_ymm;
			__m256i nums_s_ymm;
			
			nums_s_ymm = _mm256_srli_si256 (nums_ymm,1*2);	
			nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);	
			nums_s_ymm = _mm256_srli_si256  (nums_ymm,2*2);
			nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);
			nums_s_ymm = _mm256_srli_si256  (nums_ymm,4*2);
			nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);

			__m128i max_hi256 = _mm256_extracti128_si256 (nums_ymm, 0b1);
			nums_s_ymm = _mm256_broadcastw_epi16(max_hi256);
			nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);
			
			nums_ymm = _mm256_shufflelo_epi16(nums_ymm,0b0);
			nums_ymm = _mm256_shuffle_epi32 (nums_ymm,0b0);

			__m256i index_ymm = _mm256_cmpeq_epi16(nums_ymm,nums_copy_ymm);
			int mask = _mm256_movemask_epi8(index_ymm);
			int max_index = (__builtin_ffsll(mask)-1)/2;
			
			short max_local_score =  _mm256_extract_epi16 (nums_ymm, 0b0000);
			if(best_global < max_local_score){
				
				best_global = max_local_score;
				best_y = vector_len * i + (vector_len-1) - max_index;
				best_x = j - vector_len + max_index;
			}
		}
	}

	

	void SW_C_AVX (Alignment& alignment, bool debug){
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
		
		AVX::inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, zeroes_ymm, alignment);
		/******************************************************************************************************/
		
		int best_global = 0;
		int best_y = 0;
		int best_x = 0;

		for( int i = 0 ; i < height; i++){
			int offset_y = i * width * vector_len;

			// Levantar strings -------------------------------------------------------------------
			// String vertical --------------------------------------------------------------------
			str_col_ymm = AVX::leer_secuencia_columna(i, vector_len, seq2_len, seq2, zeroes_ymm, reverse_mask_xmm, shift_mask_right_xmm);
			diag1_ymm = _mm256_loadu_si256((__m256i const*) (score_matrix + offset_y));
			diag2_ymm = _mm256_loadu_si256((__m256i const*) (score_matrix + offset_y + vector_len));
			for( int j = 2; j < width ; j++){
				int offset_x = j * vector_len;
				// String horizontal ------------------------------------------------------------------
				str_row_ymm = AVX::leer_secuencia_fila(j, vector_len, width, seq1, zeroes_ymm, shift_mask_right_xmm, shift_mask_left_xmm);
				// Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
				AVX::calcular_scores(
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
				diag_score_ymm = _mm256_max_epi16(diag_score_ymm,zeroes_ymm);

				//save the max score in the right position of score matrix
				_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + offset_x), diag_score_ymm);

				if(j>=vector_len){
					v_aux[j - vector_len] =  _mm256_extract_epi16 (diag_score_ymm, 0x0);
				}

				diag1_ymm = diag2_ymm;
				diag2_ymm = diag_score_ymm;
					
				AVX::actualizar_posicion_maxima(best_global,best_x,best_y,vector_len,i,j,diag_score_ymm);
			}	
		}
		if(debug){
			alignment.matrix = score_matrix;
		}

		backtracking_C(
			score_matrix,
			alignment,
			vector_len,
			best_x,best_y,
			true,
			(score_fun_t)get_score_SSE,
			false
		);

		if(!debug) free(score_matrix);
	}

	namespace AVX512{
		union SIMDreg{
			__m128i x;
			__m256i y;
			__m512i z;
		};

		void printzmm(SIMDreg reg, bool word){
			short dw[32];
			char db[64];
			if(word){
				_mm512_storeu_si512((__m512i*)dw,reg.z);
				for(short e : dw)cout << e << "|";cout<<endl;
				
			}else{
				_mm512_storeu_si512((__m512i*)db,reg.z);
				for(char e : db){
					cout << (short)e << "|";
				}cout<<endl;
			}
		}

		SIMDreg constant_gap_mm, constant_missmatch_mm, constant_match_mm, zeroes_mm;
		SIMDreg str_row_mm, str_col_mm, left_score_mm, up_score_mm, diag_score_mm;
		SIMDreg str_reverse_mask_mm, str_shift_right_mask_mm, str_shift_left_mask_mm;
		SIMDreg str_512_unpacklo_epi8_mask_mm;
		SIMDreg score_512_rot_right_word_mask_mm;
		SIMDreg diag1_mm, diag2_mm;

		char *seq1, *seq2;
		unsigned int seq1_len, seq2_len;
		
		//each element has a 0x70 added, so after addition the most significative bit is activated for the trash characters
		short str_shift_right_mask[32] = {
			0x7FE0,0x7FE1,0x7FE2,0x7FE3,0x7FE4,0x7FE5,0x7FE6,0x7FE7,0x7FE8,0x7FE9,0x7FEA,0x7FEB,0x7FEC,0x7FED,0x7FEE,0x7FEF,
			0x7FF0,0x7FF1,0x7FF2,0x7FF3,0x7FF4,0x7FF5,0x7FF6,0x7FF7,0x7FF8,0x7FF9,0x7FFA,0x7FFB,0x7FFC,0x7FFD,0x7FFE,0x7FFF
		};
		short str_shift_left_mask[32] = {
			0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F
		};
		short str_reverse_mask[32] = {
			0x1F,0x1E,0x1D,0x1C,0x1B,0x1A,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0
		};
		uint8_t str_512_unpacklo_epi8_mask[64] = {
			0x0,0xFF,0x1,0xFF,0x2,0xFF,0x3,0xFF,0x4,0xFF,0x5,0xFF,0x6,0xFF,0x7,0xFF,0x8,0xFF,0x9,0xFF,0xA,0xFF,0xB,0xFF,0xC,0xFF,0xD,0xFF,0xE,0xFF,0xF,0xFF,
			0x10,0xFF,0x11,0xFF,0x12,0xFF,0x13,0xFF,0x14,0xFF,0x15,0xFF,0x16,0xFF,0x17,0xFF,0x18,0xFF,0x19,0xFF,0x1A,0xFF,0x1B,0xFF,0x1C,0xFF,0x1D,0xFF,0x1E,0xFF,0x1F,0xFF
		};

		short score_512_rot_right_word_mask[32] = {
			0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x0
		};

		const int vector_len = 32;
		
		int height; //cantidad de "franjas" de diagonales
		int width; //cantidad de diagonales por franja

		short* score_matrix;
		short* v_aux;

		int best_global;
		int best_y;
		int best_x;
			
		void inicializar_casos_base(Alignment& alignment){
			// llenamos el vector auxiliar
			for(int i = 0;i < width-1;i++){
				v_aux[i] = SHRT_MIN/2;
			}

			// inicializar casos base en matriz
			for(int i = 0 ; i < height ; i++){
				unsigned int offset_y = i * width * vector_len;
				_mm512_storeu_si512((__m512i*)(score_matrix + offset_y), zeroes_mm.z);
				_mm512_storeu_si512((__m512i*)(score_matrix + offset_y + vector_len), zeroes_mm.z);
			}
		}

		void leer_secuencia_columna(int i){

			if((i+1)*vector_len >= (int)seq2_len){
				// Desborde por abajo
				// Evita levantar de mas en la secuencia vertical
				// lee el tamanio del vector sin pasarse
				// y corrije shifteando
				int offset_str_col = (i+1)*vector_len - seq2_len;
			
				//simd : leer de memoria (movdqu)
				str_col_mm.y = _mm256_loadu_si256((__m256i*)(seq2 + seq2_len - vector_len));

				SIMDreg offset_str_col_mm;
				offset_str_col_mm.x = _mm_insert_epi16(offset_str_col_mm.x, offset_str_col,0);
				offset_str_col_mm.z = _mm512_broadcastw_epi16(offset_str_col_mm.x);
				offset_str_col_mm.z = _mm512_add_epi16(str_shift_right_mask_mm.z,offset_str_col_mm.z);

				str_col_mm.y = _mm256_permute4x64_epi64 (str_col_mm.y, 0b11011000); // 3|1|2|0	
				SIMDreg str_col_lo_mm;
				str_col_lo_mm.y = _mm256_unpacklo_epi8 (str_col_mm.y, zeroes_mm.y); // z|1|z|0
				SIMDreg str_col_hi_mm;
				str_col_hi_mm.y = _mm256_unpackhi_epi8 (str_col_mm.y, zeroes_mm.y); // z|3|z|2
				str_col_mm.z = _mm512_inserti64x4(str_col_lo_mm.z, str_col_hi_mm.y, 0b1);
				
				__mmask32 shift_right_mask = _mm512_movepi16_mask(offset_str_col_mm.z);
				shift_right_mask = _knot_mask32(shift_right_mask);
				str_col_mm.z = _mm512_maskz_permutexvar_epi16 (shift_right_mask, offset_str_col_mm.z, str_col_mm.z);
				//todos los elementos que sean basura van a convertirse en el valor 0xFFFF, haciendo que nunca matcheen mas adelante ni de casualidad
				str_col_mm.z = _mm512_mask_blend_epi16(shift_right_mask, offset_str_col_mm.z,str_col_mm.z);
			}else{
				//simd : leer de memoria (movdqu)
				str_col_mm.y = _mm256_loadu_si256((__m256i*)(seq2 + i * vector_len));
				str_col_mm.y = _mm256_permute4x64_epi64 (str_col_mm.y, 0b11011000); // 3|1|2|0	
				SIMDreg str_col_lo_mm;
				str_col_lo_mm.y = _mm256_unpacklo_epi8 (str_col_mm.y, zeroes_mm.y); // z|1|z|0
				SIMDreg str_col_hi_mm;
				str_col_hi_mm.y = _mm256_unpackhi_epi8 (str_col_mm.y, zeroes_mm.y); // z|3|z|2
				str_col_mm.z = _mm512_inserti64x4(str_col_lo_mm.z, str_col_hi_mm.y, 0b1);
			}
			str_col_mm.z = _mm512_permutexvar_epi16 (str_reverse_mask_mm.z, str_col_mm.z);
		
		}

		void leer_secuencia_fila(int j) {
			if(j-vector_len < 0){ //desborde por izquierda
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = vector_len - j;
				//simd : leer de memoria (movdqu)
				str_row_mm.y = _mm256_loadu_si256((__m256i*)(seq1));
				
				SIMDreg offset_str_row_mm;
				offset_str_row_mm.x = _mm_insert_epi16(offset_str_row_mm.x, offset_str_row,0);
				offset_str_row_mm.z = _mm512_broadcastw_epi16(offset_str_row_mm.x);
				offset_str_row_mm.z = _mm512_sub_epi16(str_shift_left_mask_mm.z,offset_str_row_mm.z);
				
				str_row_mm.y = _mm256_permute4x64_epi64 (str_row_mm.y, 0b11011000); // 3|1|2|0	
				SIMDreg str_row_lo_mm;
				str_row_lo_mm.y = _mm256_unpacklo_epi8 (str_row_mm.y, zeroes_mm.y); // z|1|z|0
				SIMDreg str_row_hi_mm;
				str_row_hi_mm.y = _mm256_unpackhi_epi8 (str_row_mm.y, zeroes_mm.y); // z|3|z|2
				str_row_mm.z = _mm512_inserti64x4(str_row_lo_mm.z, str_row_hi_mm.y, 0b1);

				__mmask32 shift_left_mask = _mm512_movepi16_mask(offset_str_row_mm.z);
				shift_left_mask = _knot_mask32(shift_left_mask);
				str_row_mm.z = _mm512_maskz_permutexvar_epi16 (shift_left_mask, offset_str_row_mm.z, str_row_mm.z);

			}else if(j > width-vector_len){ // desborde por derecha
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = j - (width-vector_len);
				
				str_row_mm.y = _mm256_loadu_si256((__m256i*)(seq1 + j - vector_len - offset_str_row) );
				
				SIMDreg offset_str_row_mm;
				offset_str_row_mm.x = _mm_insert_epi16(offset_str_row_mm.x, offset_str_row,0);
				offset_str_row_mm.z = _mm512_broadcastw_epi16(offset_str_row_mm.x);
				offset_str_row_mm.z = _mm512_add_epi16(str_shift_right_mask_mm.z,offset_str_row_mm.z);
				
				str_row_mm.y = _mm256_permute4x64_epi64 (str_row_mm.y, 0b11011000); // 3|1|2|0	
				SIMDreg str_row_lo_mm;
				str_row_lo_mm.y = _mm256_unpacklo_epi8 (str_row_mm.y, zeroes_mm.y); // z|1|z|0
				SIMDreg str_row_hi_mm;
				str_row_hi_mm.y = _mm256_unpackhi_epi8 (str_row_mm.y, zeroes_mm.y); // z|3|z|2
				str_row_mm.z = _mm512_inserti64x4(str_row_lo_mm.z, str_row_hi_mm.y, 0b1);

				__mmask32 shift_right_mask = _mm512_movepi16_mask(offset_str_row_mm.z);
				shift_right_mask = _knot_mask32(shift_right_mask);
				str_row_mm.z = _mm512_maskz_permutexvar_epi16 (shift_right_mask, offset_str_row_mm.z, str_row_mm.z);
			}else{ //caso feliz
				str_row_mm.y = _mm256_loadu_si256((__m256i*)(seq1 + j - vector_len));

				str_row_mm.y = _mm256_permute4x64_epi64 (str_row_mm.y, 0b11011000); // 3|1|2|0	
				SIMDreg str_row_lo_mm;
				str_row_lo_mm.y = _mm256_unpacklo_epi8 (str_row_mm.y, zeroes_mm.y); // z|1|z|0
				SIMDreg str_row_hi_mm;
				str_row_hi_mm.y = _mm256_unpackhi_epi8 (str_row_mm.y, zeroes_mm.y); // z|3|z|2
				str_row_mm.z = _mm512_inserti64x4(str_row_lo_mm.z, str_row_hi_mm.y, 0b1);
				
			}
		}

		void calcular_scores(int j){
			//left score
			left_score_mm.z = diag2_mm.z;
			left_score_mm.z = _mm512_add_epi16 (left_score_mm.z, constant_gap_mm.z);
			
			//up score
			up_score_mm.z = diag2_mm.z;
			up_score_mm.x = _mm_insert_epi16(up_score_mm.x,v_aux[j-1],0b0);
			up_score_mm.z = _mm512_permutexvar_epi16(score_512_rot_right_word_mask_mm.z, up_score_mm.z);
			up_score_mm.z = _mm512_add_epi16(up_score_mm.z, constant_gap_mm.z);
			
			//diag score

			diag_score_mm.z = diag1_mm.z;
			diag_score_mm.x = _mm_insert_epi16(diag_score_mm.x,v_aux[j-2],0b0);
			diag_score_mm.z = _mm512_permutexvar_epi16 (score_512_rot_right_word_mask_mm.z, diag_score_mm.z);
			SIMDreg cmp_match_mm;
			__mmask32 cmp_mask = _mm512_cmpeq_epi16_mask(str_col_mm.z,str_row_mm.z);
			cmp_match_mm.z = _mm512_mask_blend_epi16(cmp_mask, constant_missmatch_mm.z,constant_match_mm.z); 
			diag_score_mm.z = _mm512_add_epi16(diag_score_mm.z, cmp_match_mm.z);
		}

		void actualizar_posicion_maxima(int i,int j){
			//find the index of the maximum word in the 128bit register
			SIMDreg nums_mm =  diag_score_mm;
			SIMDreg nums_s_mm;
			
			//nums_mm = |WWWWWWWW|WWWWWWWW|WWWWWWWW|WWWWWWWW|
			nums_s_mm.z = _mm512_bsrli_epi128(nums_mm.z,1*2);	
			nums_mm.z = _mm512_max_epi16 (nums_mm.z,nums_s_mm.z);	
			//nums_mm = | W W W W| W W W W| W W W W| W W W W|
			nums_s_mm.z = _mm512_bsrli_epi128(nums_mm.z,2*2);	
			nums_mm.z = _mm512_max_epi16 (nums_mm.z,nums_s_mm.z);	
			//nums_mm = |   W   W|   W   W|   W   W|   W   W|
			nums_s_mm.z = _mm512_bsrli_epi128(nums_mm.z,4*2);	
			nums_mm.z = _mm512_max_epi16 (nums_mm.z,nums_s_mm.z);	
			//nums_mm = |       W|       W|       W|       W|
			__mmask16 mask = 0x1111;
			nums_mm.z = _mm512_mask_compress_epi32(nums_mm.z, mask, nums_mm.z);
			//nums_mm = |        |        |        | W W W W|
			nums_s_mm.x = _mm_bsrli_si128(nums_mm.x,2*2);
			nums_mm.x = _mm_max_epi16 (nums_mm.x,nums_s_mm.x);
			//nums_mm = |        |        |        |   W   W|
			nums_s_mm.x = _mm_bsrli_si128(nums_mm.x,4*2);
			nums_mm.x = _mm_max_epi16 (nums_mm.x,nums_s_mm.x);
			//nums_mm = |        |        |        |       W|

			nums_mm.z = _mm512_broadcastw_epi16 (nums_mm.x);
			//cout<<"i"<<i<<"j"<<j<<"|diag:"<<endl;
			//printzmm(nums_mm,true);	
			__mmask32 index_mask = _mm512_cmpeq_epi16_mask(nums_mm.z,diag_score_mm.z);
			int max_index = __builtin_ffs(index_mask)-1;
			
			short max_local_score =  _mm_extract_epi16(nums_mm.x, 0b0000);
			if(best_global < max_local_score){
				
				best_global = max_local_score;
				best_y = vector_len * i + (vector_len-1) - max_index;
				best_x = j - vector_len + max_index;
				//cerr<<"update: x:"<<best_x<<" y:"<<best_y<<" index_mask:"<<bitset<32>(index_mask)<<" index:"<<max_index<<" max_local_score:"<<max_local_score<<endl;
			}
		}

		void SW_C (Alignment& alignment, bool debug){
			seq1 = alignment.sequence_1->sequence;	
			seq2 = alignment.sequence_2->sequence;
			seq1_len = alignment.sequence_1->length;
			seq2_len = alignment.sequence_2->length;

			constant_gap_mm.x = _mm_insert_epi16(constant_gap_mm.x,alignment.parameters->gap,0);
			constant_gap_mm.z = _mm512_broadcastw_epi16(constant_gap_mm.x);
			constant_missmatch_mm.x = _mm_insert_epi16(constant_missmatch_mm.x,alignment.parameters->missmatch,0);
			constant_missmatch_mm.z = _mm512_broadcastw_epi16(constant_missmatch_mm.x);
			constant_match_mm.x = _mm_insert_epi16(constant_match_mm.x,alignment.parameters->match,0);
			constant_match_mm.z = _mm512_broadcastw_epi16(constant_match_mm.x);
			zeroes_mm.z = _mm512_setzero_si512();

			str_reverse_mask_mm.z = _mm512_loadu_si512 ((__m512i*)str_reverse_mask);
			str_shift_right_mask_mm.z =  _mm512_loadu_si512((__m512i*)str_shift_right_mask);
			str_shift_left_mask_mm.z =  _mm512_loadu_si512((__m512i*)str_shift_left_mask);
			str_512_unpacklo_epi8_mask_mm.z =  _mm512_loadu_si512((__m512i*)str_512_unpacklo_epi8_mask);
			score_512_rot_right_word_mask_mm.z =  _mm512_loadu_si512((__m512i*)score_512_rot_right_word_mask);
			
			height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
			width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
			int score_matrix_sz = height * width * vector_len; // largo de diagonal
				
			score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
			
			v_aux = (short*)malloc((width-1)*sizeof(short));
			//cerr<<"inicializar_casos_base ----------------------------------------------------------------------"<<endl;
			inicializar_casos_base(alignment);
			/******************************************************************************************************/
			
			best_global = 0;
			best_y = 0;
			best_x = 0;

			for( int i = 0 ; i < height ; i++){
				int offset_y = i * width * vector_len;

				// Levantar strings -------------------------------------------------------------------
				// String vertical --------------------------------------------------------------------
				leer_secuencia_columna(i);
				
				diag1_mm.z = _mm512_loadu_si512((__m512i const*) (score_matrix + offset_y));
				diag2_mm.z = _mm512_loadu_si512((__m512i const*) (score_matrix + offset_y + vector_len));
				
				for( int j = 2; j < width ; j++){
					int offset_x = j * vector_len;
					// String horizontal ------------------------------------------------------------------
					leer_secuencia_fila(j);
					
					//Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
					calcular_scores(j);
					
					diag_score_mm.z = _mm512_max_epi16(diag_score_mm.z,up_score_mm.z);
					diag_score_mm.z = _mm512_max_epi16(diag_score_mm.z,left_score_mm.z);
					diag_score_mm.z = _mm512_max_epi16(diag_score_mm.z,zeroes_mm.z);

					//save the max score in the right position of score matrix
					_mm512_storeu_si512((__m512i*)(score_matrix + offset_y + offset_x), diag_score_mm.z);

					if(j>=vector_len){
						v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_mm.x, 0x0);
					}

					actualizar_posicion_maxima(i,j);

					diag1_mm.z = diag2_mm.z;
					diag2_mm.z = diag_score_mm.z;
				}	
			}
			if(debug){
				alignment.matrix = score_matrix;
			}
			//cerr<<best_global<<" "<<best_x<<" "<<best_y<<endl;
			backtracking_C(
				score_matrix,
				alignment,
				vector_len,
				best_x,best_y,
				true,
				(score_fun_t)get_score_SSE,
				false
			);

			if(!debug) free(score_matrix);
		}
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

		if(implementation.compare("C_LIN") == 0) SW_C_LIN(*alignment, true);
		//else if(implementation.compare("C_LogicSSE") == 0)SW_C_withLogicSSE(*alignment, true);
		else if(implementation.compare("C_SSE") == 0)SW_C_SSE(*alignment, true);
		else if(implementation.compare("C_AVX512") == 0)AVX512::SW_C(*alignment, true);
		else if(implementation.compare("C_AVX") == 0)SW_C_AVX(*alignment, true);
		else if(implementation.compare("ASM_LIN") == 0)SW_ASM_LIN(alignment);
		else if(implementation.compare("ASM_SSE") == 0)SW_ASM_SSE(alignment, true);
		//else if(implementation.compare("ASM_AVX") == 0)SW_ASM_AVX(alignment, true);
		else throw "No existe la implementación ingresada.";


		//devuelvo la estructura modificada
		return alignment;
	}
}