#include <iostream>
#include "NW_C.hpp"

#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;


void NW_C_LIN(Alignment& alignment){
	
	unsigned int seq1_len = alignment.sequence_1->length;
	unsigned int seq2_len = alignment.sequence_2->length;
	char* seq1 = alignment.sequence_1->sequence;
	char* seq2 = alignment.sequence_2->sequence;
	
	short** scores = (short**)malloc((seq2_len)*sizeof(short*));
	
	for(unsigned int y = 0;y < seq2_len;y++)
		scores[y] = (short*)malloc((seq1_len)*sizeof(short));
	
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
	
	cerr << "Score Matrix" << endl;
	for(unsigned int y = 0;y < seq2_len;y++){
		for(unsigned int x = 0;x < seq1_len;x++){
			cerr << (int)scores[y][x] << " ";
		}
		cerr << endl << endl;
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

	
	for(unsigned int y = 0;y < seq2_len;y++){
		free(scores[y]);
	}
	
	free(scores);

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

void NW_C_withLogicSSE (Alignment& alignment){
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
		alignment.sequence_1->length-1,alignment.sequence_2->length-1,
		false,
		(score_fun_t)get_score_SSE,
		false
	);

	free(score_matrix);
}

void NW_C_SSE (Alignment& alignment){
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
	__m128i constant_gap_xmm;
	constant_gap_xmm = _mm_insert_epi16(constant_gap_xmm,alignment.parameters->gap,0);
	constant_gap_xmm = _mm_broadcastw_epi16(constant_gap_xmm);
	__m128i constant_missmatch_xmm = _mm_insert_epi16(constant_missmatch_xmm,alignment.parameters->missmatch,0);
	constant_missmatch_xmm = _mm_broadcastw_epi16(constant_missmatch_xmm);
	__m128i constant_match_xmm = _mm_insert_epi16(constant_match_xmm,alignment.parameters->match,0);
	constant_match_xmm = _mm_broadcastw_epi16(constant_match_xmm);

	//arrays auxiliares para el calculo de los scores
	__m128i str_row_xmm;
	__m128i str_col_xmm;
	__m128i left_score_xmm;
	__m128i up_score_xmm;
	__m128i diag_score_xmm;
	
	char identity_shift_mask[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF};

	//|0001|0000|0003|0002|0005|0004|0007|0006|0009|0008|000B|000A|000D|000C|000F|000E|
	char reverse_mask[16] = {0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1};
	__m128i reverse_mask_xmm = _mm_loadu_si128((__m128i*)reverse_mask);
	__m128i identity_shift_mask_xmm =  _mm_loadu_si128((__m128i*)identity_shift_mask);

	for( int i = 0 ; i < height ; i++){
		int offset_y = i * width * vector_len;

		if((i+1)*vector_len >= (int)seq2_len){
			
			int offset_col = (i+1)*vector_len - seq2_len;
			
			//simd : leer de memoria (movdqu)
			str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len - offset_col) );
			__m128i zeroes_xmm = _mm_setzero_si128();
			str_col_xmm = _mm_unpacklo_epi8(str_col_xmm,zeroes_xmm);
			
			//realizamos el shift con shuffle de a bytes, una mascara con que preserva la identidad del registro de byte
			//agarramos el offset, lo broadcasteamos en un registro de 128 enpaquetado de bytes
			//para realizar un shift a derecha se realiza la suma del offset y la mascara y luego con eso realizamos el shuffle b
			//simd : shift right
			
			__m128i offset_str_col_xmm = _mm_insert_epi8(offset_str_col_xmm,2*offset_col,0);
			offset_str_col_xmm = _mm_broadcastb_epi8(offset_str_col_xmm);
			offset_str_col_xmm = _mm_add_epi8(identity_shift_mask_xmm,offset_str_col_xmm);
			str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_col_xmm);
			

		}else{
			//simd : leer de memoria (movdqu)
			str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len) );
			__m128i zeroes_xmm = _mm_setzero_si128();
			str_col_xmm = _mm_unpacklo_epi8(str_col_xmm,zeroes_xmm);
			
		}

		str_col_xmm = _mm_shuffle_epi8(str_col_xmm,reverse_mask_xmm);
		
		for( int j = 2; j < width ; j++){
			int offset_x = j * vector_len;
			//emulamos simd
			if(j-vector_len < 0){ //desborde por izquierda
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = vector_len - j;
				//simd : leer de memoria (movdqu)
				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len + offset_str_row) );
				__m128i zeroes_xmm = _mm_setzero_si128();
				str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);

				//realizamos el shift con shuffle de a bytes, una mascara con que preserva la identidad del registro de byte
				//agarramos el offset, lo broadcasteamos en un registro de 128 enpaquetado de bytes
				//para realizar un shift a izquierda se realiza la resta del offset y la mascara y luego con eso realizamos el shuffle b
				//simd : shift left
				__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,2*offset_str_row,0);
				offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
				offset_str_row_xmm = _mm_sub_epi8(identity_shift_mask_xmm,offset_str_row_xmm);
				str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
				
			}else if(j > width-vector_len){ // desborde por derecha
				//simd : desplazamiento de puntero y levantar datos de memoria
				int offset_str_row = j - (width-vector_len);
				
				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len - offset_str_row) );
				__m128i zeroes_xmm = _mm_setzero_si128();
				str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);

				//realizamos el shift con shuffle de a bytes, una mascara con que preserva la identidad del registro de byte
				//agarramos el offset, lo broadcasteamos en un registro de 128 enpaquetado de bytes
				//para realizar un shift a izquierda se realiza la resta del offset y la mascara y luego con eso realizamos el shuffle b
				//simd : shift right
				__m128i offset_str_row_xmm = _mm_insert_epi8(offset_str_row_xmm,2*offset_str_row,0);
				offset_str_row_xmm = _mm_broadcastb_epi8(offset_str_row_xmm);
				offset_str_row_xmm = _mm_add_epi8(identity_shift_mask_xmm,offset_str_row_xmm);
				str_row_xmm = _mm_shuffle_epi8(str_row_xmm,offset_str_row_xmm);
				
			}else{ //caso feliz
				str_row_xmm = _mm_loadl_epi64((__m128i*)(seq1 + j - vector_len) );
				__m128i zeroes_xmm = _mm_setzero_si128();
				str_row_xmm = _mm_unpacklo_epi8(str_row_xmm,zeroes_xmm);

			}

			//left score
			left_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - vector_len));
			left_score_xmm = _mm_add_epi16(left_score_xmm, constant_gap_xmm);
			
			//up score
			up_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - vector_len));
			up_score_xmm = _mm_srli_si128(up_score_xmm, 2);
			up_score_xmm = _mm_insert_epi8(up_score_xmm,v_aux[j-1],0b1111);
			up_score_xmm = _mm_add_epi16(up_score_xmm, constant_gap_xmm);
			

			//diag score
			diag_score_xmm = _mm_loadu_si128 ((__m128i const*) (score_matrix + offset_y + offset_x - 2*vector_len));
			diag_score_xmm = _mm_srli_si128(diag_score_xmm, 2);
			diag_score_xmm = _mm_insert_epi8(diag_score_xmm,v_aux[j-1],0b1111);

			//compare the 2 strings and put the right penalty (match or missmatch) on each position
			__m128i cmp_match_xmm = str_col_xmm;
			cmp_match_xmm = _mm_cmpeq_epi16(str_col_xmm,str_row_xmm);
			str_row_xmm = _mm_andnot_si128(cmp_match_xmm,constant_missmatch_xmm);
			cmp_match_xmm = _mm_and_si128(cmp_match_xmm,constant_match_xmm);
			
			//get the max score of diag,up,left
			diag_score_xmm = _mm_add_epi16(diag_score_xmm, cmp_match_xmm);
			diag_score_xmm = _mm_add_epi16(diag_score_xmm, str_row_xmm);

			diag_score_xmm = _mm_max_epi16(diag_score_xmm,up_score_xmm);
			diag_score_xmm = _mm_max_epi16(diag_score_xmm,left_score_xmm);

			//save the max score in the right position of score matrix
			_mm_storeu_si128((__m128i*)(score_matrix + offset_y + offset_x), diag_score_xmm);

			v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_xmm, 0b0000);

		}	
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

	free(score_matrix);
}

