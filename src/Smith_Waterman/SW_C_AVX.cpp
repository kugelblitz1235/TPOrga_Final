#include "SW_C_AVX.hpp"

namespace SW{
namespace C{
namespace AVX{
    // Inicializar los valores del vector auxiliar y la matriz de puntajes
	void inicializar_casos_base(
		int width, 
		int height, 
		int vector_len, 
		short* v_aux, 
		short* score_matrix,
		__m256i zeroes_ymm, 
		Alignment& alignment){
		// Llenar vector auxiliar con un valor inicial negativo grande para no afectar los calculos
		for(int i = 0;i < width-1;i++){
			v_aux[i] = SHRT_MIN/2;
		}

		// Inicializar por cada franja las primeras 2 diagonales. 
		// Se pone en cada posicion un cero para no afectar los calculos posteriores			
		for(int i = 0 ; i < height ; i++){
			unsigned int offset_y = i * width * vector_len;
			_mm256_storeu_si256((__m256i*)(score_matrix + offset_y), zeroes_ymm);
			_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + vector_len), zeroes_ymm);
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


		
void actualizar_posicion_maxima(
		int &best_global,
		int &best_x,
		int &best_y,
		int vector_len,
		int i,
		int j,
		__m256i diag_score_ymm
	){
		__m256i nums_ymm =  diag_score_ymm;
		__m256i nums_copy_ymm = nums_ymm;
		__m256i nums_s_ymm;
																			// nums_mm = |WWWWWWWW|WWWWWWWW|
		nums_s_ymm = _mm256_srli_si256 (nums_ymm,1*2);		
		nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);					// nums_mm = | W W W W| W W W W|	
		nums_s_ymm = _mm256_srli_si256  (nums_ymm,2*2);
		nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);					// nums_mm = |   W   W|   W   W|
		nums_s_ymm = _mm256_srli_si256  (nums_ymm,4*2);
		nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);					// nums_mm = |       W|       W|

		nums_s_ymm = _mm256_permute4x64_epi64(nums_ymm, 0b00000010);
		nums_ymm = _mm256_max_epi16 (nums_ymm,nums_s_ymm);					// nums_mm = |        |       W|

		__m128i max_ymm = _mm256_extracti128_si256(nums_ymm, 0);
		nums_ymm = _mm256_broadcastw_epi16(max_ymm);
		
		// Obtener el índice del valor máximo en el registro
		__m256i index_ymm = _mm256_cmpeq_epi16(nums_ymm,nums_copy_ymm);
		int mask = _mm256_movemask_epi8(index_ymm);
		int max_index = __builtin_ffsll(mask)/2;							// Obtener la posicion del word de valor máximo entre todos los de la diagonal
		
		short max_local_score =  _mm256_extract_epi16 (nums_ymm, 0b0000);	// Obtener el valor máximo
		if(best_global < max_local_score){
			
			best_global = max_local_score;
			best_y = vector_len * i + (vector_len-1) - max_index;
			best_x = j - vector_len + max_index;
		}
	}

	void SW(Alignment& alignment, bool debug){
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
		
		inicializar_casos_base(width, height, vector_len, v_aux, score_matrix, zeroes_ymm, alignment);
		/******************************************************************************************************/
		
		int best_global = 0;								// Máximo puntaje obtenido
		int best_x = 0;										// Coordenada x del máximo puntaje
		int best_y = 0;										// Coordenada y del máximo puntaje

		for( int i = 0 ; i < height; i++){
			int offset_y = i * width * vector_len;

			// Levantar strings -------------------------------------------------------------------
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
				
				// Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba, diagonalmente y cero
				diag_score_ymm = _mm256_max_epi16(diag_score_ymm,up_score_ymm);
				diag_score_ymm = _mm256_max_epi16(diag_score_ymm,left_score_ymm);
				diag_score_ymm = _mm256_max_epi16(diag_score_ymm,zeroes_ymm);

				// Almacenamos el puntaje máximo en la posición correcta de la matriz
				_mm256_storeu_si256((__m256i*)(score_matrix + offset_y + offset_x), diag_score_ymm);

				if(j>=vector_len){
					v_aux[j - vector_len] = _mm256_extract_epi16 (diag_score_ymm, 0x0);
				}
				
				// Actualizar las 2 diagonales anteriores para la siguiente iteracion, actualizando diag2 con el valor de la actual
				diag1_ymm = diag2_ymm;
				diag2_ymm = diag_score_ymm;
					
				actualizar_posicion_maxima(best_global,best_x,best_y,vector_len,i,j,diag_score_ymm);
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
			best_x,best_y,
			true,
			(score_fun_t)get_score_SSE,
			false
		);

		if(!debug) free(score_matrix);
	}
}
}
}