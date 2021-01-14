#include "SW_C_SSE.hpp"

namespace SW{
namespace C{
namespace SSE{

    // Equivalentes a registros nombrados:
    __m128i constant_gap_xmm, constant_missmatch_xmm, constant_match_xmm, zeroes_xmm;
    __m128i str_row_xmm, str_col_xmm, left_score_xmm, up_score_xmm, diag_score_xmm;
    __m128i reverse_mask_xmm;
    __m128i diag1_xmm, diag2_xmm;
    
    char* seq1;
    char* seq2;
    unsigned int seq1_len;
    unsigned int seq2_len;

    // Máscara utilizada para invertir el string almacenado en un registro
    char reverse_mask[16] = {0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1};
    
    int vector_len = 8;

    int height;                                         // Cantidad de "franjas" de diagonales
    int width; 						                // Cantidad de diagonales por franja
    int score_matrix_sz; 					            // Cantidad de celdas de la matriz

    // Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores
    short* score_matrix;
    short* v_aux;


	int best_global;									// Máximo puntaje obtenido
	int best_x;											// Coordenada x del máximo puntaje
	int best_y;											// Coordenada y del máximo puntaje

    // Inicializar los valores del vector auxiliar y la matriz de puntajes
	void inicializar_casos_base() {
		// Llenar vector auxiliar con un valor inicial negativo grande para no afectar los calculos
		for(int i = 0;i < width-1;i++){
			v_aux[i] = SHRT_MIN/2;
		}

		// Inicializar por cada franja las primeras 2 diagonales. 
		// Se pone en cada posicion cero para no afectar los calculos posteriores
		for(int i = 0 ; i < height ; i++){
			unsigned int offset_y = i * width * vector_len;
			_mm_storeu_si128((__m128i*)(score_matrix + offset_y), zeroes_xmm);
			_mm_storeu_si128((__m128i*)(score_matrix + offset_y + vector_len), zeroes_xmm);
		}
	}

	// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
	void leer_secuencia_columna(int i){
		if((i+1)*vector_len >= (int)seq2_len){ 	// Caso de desborde por abajo

			// Desborde por abajo
			// Evita levantar de mas en la secuencia vertical
			// lee el tamanio del vector sin pasarse
			// y corrije shifteando
			int offset_col = (i+1)*vector_len - seq2_len;										// Indica cuanto hay que shiftear el string columna luego de levantarlo
			//simd : leer de memoria (movdqu)
			str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len - offset_col) );

			__m128i shift_count_xmm = zeroes_xmm;												    // Se utiliza junto con offset_col para shiftear str_col_xmm
			__m128i shift_mask_xmm = zeroes_xmm;
			// Acomodar los caracteres en str_col_xmm para que queden en las posiciones correctas para la comparación mas adelante
			shift_count_xmm = _mm_insert_epi8(shift_count_xmm, offset_col*8, 0);
			str_col_xmm = _mm_srl_epi64(str_col_xmm, shift_count_xmm);								// str_col_xmm = |0...0|str_col|

			// shift_mask_xmm va a tener 1's donde haya caracteres inválidos
			shift_mask_xmm = _mm_cmpeq_epi8(shift_mask_xmm, shift_mask_xmm);
			shift_count_xmm = _mm_insert_epi8(shift_count_xmm, (char)(8-offset_col)*8, 0);
			shift_mask_xmm = _mm_sll_epi64(shift_mask_xmm, shift_count_xmm);						// shift_mask_xmm = |1...1|0...0|

			// Combinar la mascara de caracteres invalidos con los caracteres validos
			str_col_xmm = _mm_or_si128(str_col_xmm, shift_mask_xmm);								// str_col_xmm = |1...1|str_col|
			
		}else{// Caso sin desborde
			// Levantar directamente de memoria, no hay desborde
			str_col_xmm = _mm_loadl_epi64((__m128i*)(seq2 + i * vector_len) );
			
		}
		// Desempaquetar los caracteres en str_col_xmm para trabajar con words
		str_col_xmm = _mm_unpacklo_epi8(str_col_xmm, zeroes_xmm);
		// Invertir la secuencia de caracteres para compararlos correctamente mas adelante
		str_col_xmm = _mm_shuffle_epi8(str_col_xmm,reverse_mask_xmm);
	}

	// Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
	void leer_secuencia_fila(int j) {
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
	}


	// Calcula los puntajes resultantes de las comparaciones entre caracteres
    void calcular_scores(int j){

        // Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
        left_score_xmm = diag2_xmm;
        left_score_xmm = _mm_add_epi16(left_score_xmm, constant_gap_xmm);
        
        // Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
        up_score_xmm = diag2_xmm;
        up_score_xmm = _mm_srli_si128(up_score_xmm, 2);
        up_score_xmm = _mm_insert_epi16(up_score_xmm,v_aux[j-1],0b111);
        up_score_xmm = _mm_add_epi16(up_score_xmm, constant_gap_xmm);
        
        // Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
        // si coinciden o no los caracteres de la fila y columna correspondientes
        diag_score_xmm = diag1_xmm;
        diag_score_xmm = _mm_srli_si128(diag_score_xmm, 2);
        diag_score_xmm = _mm_insert_epi16(diag_score_xmm,v_aux[j-2],0b111);
        // Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
        __m128i cmp_match_xmm;
        cmp_match_xmm = _mm_cmpeq_epi16(str_col_xmm,str_row_xmm); 									// Mascara con unos en las posiciones donde coinciden los caracteres
        cmp_match_xmm = _mm_blendv_epi8(constant_missmatch_xmm,constant_match_xmm,cmp_match_xmm); 	// Seleccionar para cada posicion el puntaje correcto basado en la mascara previa
        diag_score_xmm = _mm_add_epi16(diag_score_xmm, cmp_match_xmm);
    }

	void actualizar_posicion_maxima(int i,int j){
		// Encontrar el índice del máximo word en el registro xmm
		__m128i nums_xmm =  diag_score_xmm;
		__m128i nums_copy_xmm = nums_xmm;
		__m128i nums_s_xmm;
		
														// nums_mm = |WWWWWWWW|	
		nums_s_xmm = _mm_srli_si128(nums_xmm,1*2);	
		nums_xmm = _mm_max_epi16(nums_xmm,nums_s_xmm);  // nums_mm = | W W W W|	
		nums_s_xmm = _mm_srli_si128 (nums_xmm,2*2);
		nums_xmm = _mm_max_epi16(nums_xmm,nums_s_xmm);	// nums_mm = |   W   W|
		nums_s_xmm = _mm_srli_si128 (nums_xmm,4*2);
		nums_xmm = _mm_max_epi16(nums_xmm,nums_s_xmm);  // nums_mm = |       W|
		
		// Broadcastear el máximo en todo el registro
		nums_xmm = _mm_shufflelo_epi16(nums_xmm,0b0);
		nums_xmm = _mm_shuffle_epi32 (nums_xmm,0b0);

		// Obtener el índice del valor máximo en el registro
		__m128i index_xmm = _mm_cmpeq_epi16(nums_xmm,nums_copy_xmm);
		index_xmm = _mm_packs_epi16(index_xmm,index_xmm);
		int64_t index_mask = _mm_extract_epi64(index_xmm,0);
		
		int max_index = __builtin_ffsll(index_mask)/8;						// Obtener la posicion del word de valor máximo entre todos los de la diagonal
		short max_local_score =  _mm_extract_epi16 (nums_xmm, 0b0000);		// Obtener el valor máximo
		if(best_global < max_local_score){
			
			best_global = max_local_score;
			best_y = vector_len * i + (vector_len-1) - max_index;
			best_x = j - vector_len + max_index;
		}
	}

	void SW(Alignment& alignment, bool debug){
		// Tamaños y punteros a las secuencias
        seq1 = alignment.sequence_1->sequence;	
        seq2 = alignment.sequence_2->sequence;
        seq1_len = alignment.sequence_1->length;
        seq2_len = alignment.sequence_2->length;
        
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
        
        reverse_mask_xmm = _mm_loadu_si128((__m128i*)reverse_mask);
        
        // El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
        vector_len = 8;
        
        height = ((seq2_len + vector_len - 1)/ vector_len); 			// Cantidad de "franjas" de diagonales
        width = (1 + seq1_len + vector_len - 1); 						// Cantidad de diagonales por franja
        score_matrix_sz = height * width * vector_len; 					// Cantidad de celdas de la matriz
        
        // Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores
        score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
        v_aux = (short*)malloc((width-1)*sizeof(short));
            
        inicializar_casos_base();

	/******************************************************************************************************/

		best_global = 0;								// Máximo puntaje obtenido
		best_x = 0;										// Coordenada x del máximo puntaje
		best_y = 0;										// Coordenada y del máximo puntaje
		
		for( int i = 0 ; i < height ; i++){
			int offset_y = i * width * vector_len;

			// Levantar strings -------------------------------------------------------------------
			// String vertical --------------------------------------------------------------------
			leer_secuencia_columna(i);
			
			// Cargar las primeras 2 diagonales de la franja actual necesarios para el primer calculo dentro de la franja
			diag1_xmm = _mm_loadu_si128((__m128i const*) (score_matrix + offset_y));
			diag2_xmm = _mm_loadu_si128((__m128i const*) (score_matrix + offset_y + vector_len));
                
			for( int j = 2; j < width ; j++){
				int offset_x = j * vector_len;
				// String horizontal ------------------------------------------------------------------
				leer_secuencia_fila(j);
				// Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
				calcular_scores(j);
				// Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba, diagonalmente y cero
				diag_score_xmm = _mm_max_epi16(diag_score_xmm,up_score_xmm);
				diag_score_xmm = _mm_max_epi16(diag_score_xmm,left_score_xmm);
				diag_score_xmm = _mm_max_epi16(diag_score_xmm,zeroes_xmm);

				// Almacenamos el puntaje máximo en la posición correcta de la matriz
				_mm_storeu_si128((__m128i*)(score_matrix + offset_y + offset_x), diag_score_xmm);
				
				if(j>=vector_len){
					v_aux[j - vector_len] =  _mm_extract_epi16 (diag_score_xmm, 0b0000);
				}
                // Actualizar las 2 diagonales anteriores para la siguiente iteracion, actualizando diag2 con el valor de la actual
                diag1_xmm = diag2_xmm;
                diag2_xmm = diag_score_xmm;
				
				// Actualizamos el valor máximo y su posicion 
				actualizar_posicion_maxima(i,j);
			}	
		}

		if(debug){// Utilizar para debuggear los valores en la matriz de puntajes
			alignment.matrix = score_matrix;
		}

		// Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion correspondiente al valor máximo obtenido
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