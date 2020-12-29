#include <iostream>
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
		alignment.matrix = new_alignment_matrix(1, seq1_len, seq2_len);
		alignment.matrix->matrix = (short *)scores;
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

void SW_C_withLogicSSE (Alignment& alignment, bool debug){
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
					score_matrix[offset_y + offset_x + k] = 0;
				else
					score_matrix[offset_y + offset_x + k] = SHRT_MIN/2;
				count++;
			}			
		}
	}
	/******************************************************************************************************/
	short best_global = 0;
	unsigned int best_x;
	unsigned int best_y;

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
		alignment.matrix = new_alignment_matrix(vector_len, seq1_len, seq2_len);
		alignment.matrix->matrix = score_matrix;

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

void inicializar_casos_base(int width, int height, int vector_len, short* v_aux, short* score_matrix) {
	//llenamos el vector auxiliar
	for(int i = 0;i < width-1;i++){
		v_aux[i] = SHRT_MIN/2;
	}

	for(int i = 0 ; i < height ; i++){
		unsigned int offset_y = i * width * vector_len;
		// for( int j = 0; j < 2 ; j++){
		// 	unsigned int offset_x = j * vector_len;
		// 	//emulamos simd
		// 	for( int k = 0;k < vector_len;k++){
		// 		if( j==1 && k == vector_len-1)
		// 			score_matrix[offset_y + offset_x + k] = 0;
		// 		else
		// 			score_matrix[offset_y + offset_x + k] = SHRT_MIN/2;
		// 	}			
		// }
		__m128i diag;
		diag = _mm_insert_epi16(diag,SHRT_MIN/2,0);
		diag = _mm_broadcastw_epi16(diag);
		_mm_storeu_si128((__m128i*)(score_matrix + offset_y), diag);
		diag = _mm_insert_epi16(diag, 0, 7);
		_mm_storeu_si128((__m128i*)(score_matrix + offset_y + vector_len), diag);
	}
}

__m128i leer_secuencia_columna(
	int i,
	int vector_len,
	int seq2_len,
	char* seq2,
	__m128i zeroes_xmm,
	__m128i shift_mask_col_xmm,
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
		str_col_xmm = _mm_unpacklo_epi8(str_col_xmm,zeroes_xmm);
		//realizamos el shift con shuffle de a bytes, una mascara con que preserva la identidad del registro de byte
		//agarramos el offset, lo broadcasteamos en un registro de 128 enpaquetado de bytes
		//para realizar un shift a derecha se realiza la suma del offset y la mascara y luego con eso realizamos el shuffle b
		//simd : shift right
		__m128i offset_str_col_xmm = _mm_insert_epi8(offset_str_col_xmm,2*offset_col,0);
		offset_str_col_xmm = _mm_broadcastb_epi8(offset_str_col_xmm);
		offset_str_col_xmm = _mm_add_epi8(shift_mask_col_xmm,offset_str_col_xmm);
		//las posiciones que se correspondan con caracteres basura (que no existen) van a tener un 1 en la posicion mas significativa
		__m128i ones_mask = _mm_srai_epi16 (offset_str_col_xmm, 15);
		ones_mask = _mm_slli_epi16(ones_mask,15);
		
		str_col_xmm = _mm_shuffle_epi8(str_col_xmm,offset_str_col_xmm);
		//todos los elementos que sean basura van a convertirse en el valor 0xFFFF, haciendo que nunca matcheen mas adelante ni de casualidad
		str_col_xmm = _mm_or_si128(str_col_xmm,ones_mask);
		

	}else{
		//simd : leer de memoria (movdqu)
		str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len) );
		str_col_xmm = _mm_unpacklo_epi8(str_col_xmm,zeroes_xmm);
		
	}

	// Reverse de la secuencia vertical
	str_col_xmm = _mm_shuffle_epi8(str_col_xmm,reverse_mask_xmm);
	
	return str_col_xmm;
}


__m128i leer_secuencia_fila(
	int j,
	int vector_len,
	int width,
	char* seq1,
	__m128i zeroes_xmm,
	__m128i shift_mask_row_xmm
) {
	__m128i str_row_xmm;
	
	if(j-vector_len < 0){ //desborde por izquierda
			//simd : desplazamiento de puntero y levantar datos de memoria
			int offset_str_row = vector_len - j;
			//simd : leer de memoria (movdqu)
			str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len + offset_str_row) );
			str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);
			//realizamos el shift con shuffle de a bytes, una mascara con que preserva la identidad del registro de byte
			//agarramos el offset, lo broadcasteamos en un registro de 128 enpaquetado de bytes
			//para realizar un shift a izquierda se realiza la resta del offset y la mascara y luego con eso realizamos el shuffle b
			//simd : shift left
			__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,2*offset_str_row,0);
			offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
			offset_str_row_xmm = _mm_sub_epi8(shift_mask_row_xmm,offset_str_row_xmm);
			str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);

	}else if(j > width-vector_len){ // desborde por derecha
			//simd : desplazamiento de puntero y levantar datos de memoria
			int offset_str_row = j - (width-vector_len);
			
			str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len - offset_str_row) );
			str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);
			//realizamos el shift con shuffle de a bytes, una mascara con que preserva la identidad del registro de byte
			//agarramos el offset, lo broadcasteamos en un registro de 128 enpaquetado de bytes
			//para realizar un shift a izquierda se realiza la resta del offset y la mascara y luego con eso realizamos el shuffle b
			//simd : shift right
			__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,2*offset_str_row,0);
			offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
			offset_str_row_xmm = _mm_add_epi8(shift_mask_row_xmm,offset_str_row_xmm);
			str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
			
	}else{ //caso feliz
			str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len) );
			str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);
	}
	
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
	__m128i cmp_match_xmm = str_col_xmm;
	cmp_match_xmm = _mm_cmpeq_epi16(str_col_xmm,str_row_xmm);
	str_row_xmm = _mm_andnot_si128(cmp_match_xmm,constant_missmatch_xmm);
	cmp_match_xmm = _mm_and_si128(cmp_match_xmm,constant_match_xmm);
	
	//get the max score of diag,up,left
	diag_score_xmm = _mm_add_epi16(diag_score_xmm, cmp_match_xmm);
	diag_score_xmm = _mm_add_epi16(diag_score_xmm, str_row_xmm);

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
	
	nums_xmm = _mm_broadcastw_epi16(nums_xmm);
	
	__m128i index_xmm = _mm_cmpeq_epi16(nums_xmm,nums_copy_xmm);
	index_xmm = _mm_packs_epi16(index_xmm,index_xmm);
	int64_t index_mask = _mm_extract_epi64(index_xmm,0);
	
	int max_index = __builtin_ffsll(index_mask)/8;
	short max_local_score =  _mm_extract_epi16 (nums_xmm, 0b0000);
	if(best_global < max_local_score){
		
		best_global = max_local_score;
		best_y = vector_len * i + (vector_len-1) - max_index;
		best_x = j - vector_len + max_index;
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
	__m128i reverse_mask_xmm, shift_mask_col_xmm, shift_mask_row_xmm;

	constant_gap_xmm = _mm_insert_epi16(constant_gap_xmm,alignment.parameters->gap,0);
	constant_gap_xmm = _mm_broadcastw_epi16(constant_gap_xmm);
	constant_missmatch_xmm = _mm_insert_epi16(constant_missmatch_xmm,alignment.parameters->missmatch,0);
	constant_missmatch_xmm = _mm_broadcastw_epi16(constant_missmatch_xmm);
	constant_match_xmm = _mm_insert_epi16(constant_match_xmm,alignment.parameters->match,0);
	constant_match_xmm = _mm_broadcastw_epi16(constant_match_xmm);
	zeroes_xmm = _mm_setzero_si128();

	//each element has a 0x70 added, so after addition the most significative bit is activated for the trash characters
	char shift_mask_col[16] = {0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F};
	char shift_mask_row[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF};

	//|0001|0000|0003|0002|0005|0004|0007|0006|0009|0008|000B|000A|000D|000C|000F|000E|
	char reverse_mask[16] = {0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1};
	reverse_mask_xmm = _mm_loadu_si128((__m128i*)reverse_mask);
	shift_mask_col_xmm =  _mm_loadu_si128((__m128i*)shift_mask_col);
	shift_mask_row_xmm =  _mm_loadu_si128((__m128i*)shift_mask_row);
	
	//en este caso hardcodeamos el tamaño del vector
	int vector_len = 8;
	
	int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
	int width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
	int score_matrix_sz = height * width * vector_len; // largo de diagonal
		
	short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
	
	short* v_aux = (short*)malloc((width-1)*sizeof(short));
	
	inicializar_casos_base(width, height, vector_len, v_aux, score_matrix);

/******************************************************************************************************/

	int best_global=0;
	int best_y;
	int best_x;

	for( int i = 0 ; i < height ; i++){
		int offset_y = i * width * vector_len;

		// Levantar strings -------------------------------------------------------------------
		// String vertical --------------------------------------------------------------------
		str_col_xmm = leer_secuencia_columna(i, vector_len, seq2_len, seq2, zeroes_xmm, shift_mask_col_xmm, reverse_mask_xmm);

		for( int j = 2; j < width ; j++){
			int offset_x = j * vector_len;
			//emulamos simd
			str_row_xmm = leer_secuencia_fila(j, vector_len, width, seq1, zeroes_xmm, shift_mask_row_xmm);
			
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
			diag_score_xmm = _mm_max_epi16(diag_score_xmm,zeroes_xmm);

			//save the max score in the right position of score matrix
			_mm_storeu_si128((__m128i*)(score_matrix + offset_y + offset_x), diag_score_xmm);
			
			if(j>=vector_len){
				v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_xmm, 0b0000);
			}
			
			actualizar_posicion_maxima(best_global,best_x,best_y,vector_len,i,j,diag_score_xmm);
		}	
	}

	// printScoreMatrix2(score_matrix,&alignment,vector_len);
	if(debug){
		alignment.matrix = new_alignment_matrix(vector_len, seq1_len, seq2_len);
		alignment.matrix->matrix = score_matrix;
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
		//ejecuto la implementación en c sin debug
		SW_C_LIN(*alignment, false);
	
	}else if(implementation.compare("LIN") == 0){
		
		//ejecuto el algoritmo en asm lineal
		SW_ASM_LIN(alignment);
	}
	else{
		throw "No existe la implementación ingresada.";
	}


	//devuelvo la estructura modificada
	return alignment;
}
}