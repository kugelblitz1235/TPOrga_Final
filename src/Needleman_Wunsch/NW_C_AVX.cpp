#include "NW_C_AVX.hpp"

namespace NW{
namespace C{
namespace AVX{
    // Directiva utilizada para emular la representación de un ymm con su parte xmm
    // |    ymm    |
    // | 128 | xmm |
    union SIMDreg{
        __m128i x;
        __m256i y;
    };

    // Tamaños y punteros a las secuencias
    char* seq1;	
    char* seq2;
    unsigned int seq1_len;
    unsigned int seq2_len;
    
    // Equivalentes a registros nombrados:
    SIMDreg constant_gap_mm, constant_missmatch_mm, constant_match_mm, zeroes_mm;
    SIMDreg str_row_mm, str_col_mm, left_score_mm, up_score_mm, diag_score_mm;
    SIMDreg reverse_mask_mm, shift_mask_right_mm, shift_mask_left_mm;
    SIMDreg diag1_mm, diag2_mm;

    // Cada posicion en la máscara tiene 0x70 sumado, para que luego de la suma el bit más significativo esté activado para las posiciones que caen fuera de los límites
    char shift_mask_right[16] = {0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F};
    char shift_mask_left[16] = {0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF};
    char reverse_mask[16] = {0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0};
   
    // El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
    int vector_len = 16;
    
    int height;                                         // Cantidad de "franjas" de diagonales
    int width; 						                // Cantidad de diagonales por franja
    int score_matrix_sz; 					            // Cantidad de celdas de la matriz

    // Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores
    short* score_matrix;
    short* v_aux;
    
    // Inicializar los valores del vector auxiliar y la matriz de puntajes
    void inicializar_casos_base(Alignment& alignment){
        // Llenar vector auxiliar con un valor inicial negativo grande para no afectar los calculos
        for(int i = 0;i < width-1;i++){
            v_aux[i] = SHRT_MIN;
        }

        SIMDreg diag_mm;
        // Inicializar por cada franja las primeras 2 diagonales. 
        // Se pone en cada posicion un valor negativo grande para no afectar los calculos posteriores
        for(int i = 0 ; i < height ; i++){
            unsigned int offset_y = i * width * vector_len;
            // Broadcastear el valor SHRT_MIN/2 a nivel word en el registro diag
            diag_mm.x = _mm_insert_epi16(diag_mm.x,SHRT_MIN,0);
            diag_mm.y = _mm256_broadcastw_epi16(diag_mm.x);
            _mm256_storeu_si256((__m256i*)(score_matrix + offset_y), diag_mm.y);
            // Insertar la penalidad adecuada para la franja en la posicion mas alta de la diagonal
            diag_mm.y = _mm256_insert_epi16(diag_mm.y,i * vector_len * alignment.parameters->gap, 15); // 15 = vector_len - 1
            _mm256_storeu_si256((__m256i*)(score_matrix + offset_y + vector_len), diag_mm.y);
        }
    }

    // Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
    void leer_secuencia_columna(int i){
     
        if((i+1)*vector_len >= (int)seq2_len){// Caso de desborde por abajo
            int offset_str_col = (i+1)*vector_len - seq2_len;										// Indica cuanto hay que shiftear el string columna luego de levantarlo
                                                                                        
            str_col_mm.x = _mm_loadu_si128((__m128i*)(seq2 + seq2_len - vector_len));				// str_col_mm.x = |0|str_col|	

            // Poner el valor offset_str_col en todos los bytes de offset_str_col_mm.x, 
            // el cual indica cuantos bytes hay que shiftear los caracteres para dejarlos en la posición correcta
            SIMDreg offset_str_col_mm;
            offset_str_col_mm.x = _mm_insert_epi8(offset_str_col_mm.x,offset_str_col,0);
            offset_str_col_mm.x = _mm_broadcastb_epi8(offset_str_col_mm.x);						// offset_str_col_mm.x = |offset|...|offset|
            // Sumar con shift_mask_right_mm.x para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
            offset_str_col_mm.x = _mm_add_epi8(shift_mask_right_mm.x,offset_str_col_mm.x);			// offset_str_col_mm.x = |0x7F + offset|0x7E + offset|...|0x70 + offset|
            // Shiftear utilizando la mascara offset_str_col_mm.x 
            str_col_mm.x = _mm_shuffle_epi8(str_col_mm.x,offset_str_col_mm.x);
            // Todos los elementos que sean basura van a convertirse en un valor que arranca con 1, haciendo que nunca matcheen mas adelante
            str_col_mm.x = _mm_blendv_epi8(str_col_mm.x,offset_str_col_mm.x, offset_str_col_mm.x);	// str_col_mm.x = |1...1|str_col|
        }else{ // Caso sin desborde
            str_col_mm.x = _mm_loadu_si128((__m128i*)(seq2 + i * vector_len));
        }

        // Invertir el string almacenado en str_col_mm.x
        str_col_mm.x = _mm_shuffle_epi8 (str_col_mm.x, reverse_mask_mm.x);
        // Desempaquetar los caracteres almacenados en str_col_mm.x para trabajar con words
        SIMDreg str_col_hi_mm;
        str_col_hi_mm.x = _mm_unpackhi_epi8(str_col_mm.x, zeroes_mm.x);
        SIMDreg str_col_lo_mm;
        str_col_lo_mm.x = _mm_unpacklo_epi8(str_col_mm.x, zeroes_mm.x);
        str_col_mm.y = _mm256_set_m128i(str_col_hi_mm.x,str_col_lo_mm.x);
    }

    // Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
    void leer_secuencia_fila(int j) {
        if(j-vector_len < 0){ // Caso de desborde por izquierda
            int offset_str_row = vector_len - j;													// Indica cuanto hay que shiftear a la izquierda el string fila luego de levantarlo
            str_row_mm.x = _mm_loadu_si128((__m128i*)(seq1) );										
            
            // Broadcastear el offset_str_row en offset_str_row_mm.x 
            SIMDreg offset_str_row_mm;
            offset_str_row_mm.x = _mm_insert_epi8(offset_str_row_mm.x,offset_str_row,0);
            offset_str_row_mm.x = _mm_broadcastb_epi8(offset_str_row_mm.x);							// offset_str_col_mm.x = |offset|...|offset|
            // Sumar con shift_mask_left_mm.x para indicar cuanto hay que shiftear a la izquierda los caracteres para que queden en las posiciones correctas
            offset_str_row_mm.x = _mm_sub_epi8(shift_mask_left_mm.x,offset_str_row_mm.x);				// offset_str_col_mm.x = |0xF - offset|0xE - offset|...|0x0 - offset|
            
            // Acomodar los caracteres en str_row_mm.x para que queden en las posiciones correctas para la comparación mas adelante
            str_row_mm.x = _mm_shuffle_epi8(str_row_mm.x,offset_str_row_mm.x);							// str_row_mm.x = |str_row|0...0|
                
        }else if(j > width-vector_len){ // Caso de desborde por derecha
            //Desplazamiento de puntero a derecha y levantar datos de memoria
            int offset_str_row = j - (width-vector_len);											// Indica cuanto hay que shiftear a la derecha el string fila luego de levantarlo
            str_row_mm.x = _mm_loadu_si128((__m128i*)(seq1 + j - vector_len - offset_str_row) );
            
            // Broadcastear el offset_str_row en offset_str_row_mm.x 
            SIMDreg offset_str_row_mm;
            offset_str_row_mm.x = _mm_insert_epi8(offset_str_row_mm.x,offset_str_row,0);
            offset_str_row_mm.x = _mm_broadcastb_epi8(offset_str_row_mm.x);							// offset_str_col_mm.x = |offset|...|offset|
            // Sumar con shift_mask_right_mm.x para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
            offset_str_row_mm.x = _mm_add_epi8(shift_mask_right_mm.x,offset_str_row_mm.x);				// offset_str_col_mm.x = |0x7F + offset|0x7E + offset|...|0x70 + offset|
            // Acomodar los caracteres en str_row_mm.x para que queden en las posiciones correctas para la comparación mas adelante
            str_row_mm.x = _mm_shuffle_epi8(str_row_mm.x,offset_str_row_mm.x);							// str_row_mm.x = |0...0|str_row|
        
        }else{ // Caso sin desborde
            str_row_mm.x = _mm_loadu_si128((__m128i*)(seq1 + j - vector_len));
        }
        
        // Desempaquetar los caracteres en str_row_mm.x para trabajar con words
        SIMDreg str_row_hi_mm;
        str_row_hi_mm.x = _mm_unpackhi_epi8(str_row_mm.x, zeroes_mm.x);
        SIMDreg str_row_lo_mm;
        str_row_lo_mm.x = _mm_unpacklo_epi8(str_row_mm.x, zeroes_mm.x);
        str_row_mm.y = _mm256_set_m128i(str_row_hi_mm.x,str_row_lo_mm.x);
    }

    //Calcula los puntajes resultantes de las comparaciones entre caracteres
    void calcular_scores(int j){

        // Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
        left_score_mm.y = diag2_mm.y;
        left_score_mm.y = _mm256_adds_epi16 (left_score_mm.y, constant_gap_mm.y);
        
        // Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
        up_score_mm.y = diag2_mm.y;
        // El shift de a bytes se hace en cada linea de 128 bits,
        // por lo que hay que mover manualmente el word en la posicion mas baja en la linea de 128 bits
        // mas alta, e insertarlo en la posicion mas alta de la linea de 128 bits mas baja
        short medium_score = _mm256_extract_epi16( up_score_mm.y, 0b1000);					// |0xF...0x8|0x7...0x0| -> 0x8
        up_score_mm.y = _mm256_bsrli_epi128(up_score_mm.y, 2); 								// |0xF...0x8|0x7...0x0| -> |0x0...0x9|0x0...0x1|
        up_score_mm.y = _mm256_insert_epi16(up_score_mm.y,medium_score,0b0111); 				// |0x0...0x9|0x0...0x1| -> |0x0...0x9|0x8...0x1|
        up_score_mm.y = _mm256_insert_epi16(up_score_mm.y,v_aux[j-1],15); 					// 15 = vector_len - 1
        up_score_mm.y = _mm256_adds_epi16(up_score_mm.y, constant_gap_mm.y);
        
        // Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
        // si coinciden o no los caracteres de la fila y columna correspondientes
        diag_score_mm.y = diag1_mm.y;                                             
        medium_score = _mm256_extract_epi16( diag_score_mm.y, 0b1000);                       // |0xF...0x8|0x7...0x0| -> 0x8
        diag_score_mm.y = _mm256_bsrli_epi128(diag_score_mm.y, 2);                            // |0xF...0x8|0x7...0x0| -> |0x0...0x9|0x0...0x1|
        diag_score_mm.y = _mm256_insert_epi16(diag_score_mm.y,medium_score,0b0111);           // |0x0...0x9|0x0...0x1| -> |0x0...0x9|0x8...0x1|
        
        diag_score_mm.y = _mm256_insert_epi16(diag_score_mm.y,v_aux[j-2],15); 				// 15 = vector_len - 1
        
        // Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
        SIMDreg cmp_match_mm;
        cmp_match_mm.y = _mm256_cmpeq_epi16(str_col_mm.y,str_row_mm.y);
        cmp_match_mm.y = _mm256_blendv_epi8(constant_missmatch_mm.y,constant_match_mm.y,cmp_match_mm.y); 
        diag_score_mm.y = _mm256_adds_epi16(diag_score_mm.y, cmp_match_mm.y);

    }

    void NW (Alignment& alignment, bool debug){
        // Tamaños y punteros a las secuencias
        seq1 = alignment.sequence_1->sequence;	
        seq2 = alignment.sequence_2->sequence;
        seq1_len = alignment.sequence_1->length;
        seq2_len = alignment.sequence_2->length;
       
        // Broadcastear el valor de gap, a nivel word, en el registro
        constant_gap_mm.x = _mm_insert_epi16(constant_gap_mm.x,alignment.parameters->gap,0);
        constant_gap_mm.y = _mm256_broadcastw_epi16(constant_gap_mm.x);
        // Broadcastear el valor de missmatch, a nivel word, en el registro
        constant_missmatch_mm.x = _mm_insert_epi16(constant_missmatch_mm.x,alignment.parameters->missmatch,0);
        constant_missmatch_mm.y = _mm256_broadcastw_epi16(constant_missmatch_mm.x);
        // Broadcastear el valor de match, a nivel word, en el registro
        constant_match_mm.x = _mm_insert_epi16(constant_match_mm.x,alignment.parameters->match,0);
        constant_match_mm.y = _mm256_broadcastw_epi16(constant_match_mm.x);
        // Máscara de ceros
        zeroes_mm.y = _mm256_setzero_si256();
        
        // Máscara utilizada para invertir el string almacenado en un registro
        reverse_mask_mm.x = _mm_loadu_si128 ((__m128i*)reverse_mask);
        // Máscara utilizada para shiftear a derecha los caracteres 
        shift_mask_right_mm.x =  _mm_loadu_si128((__m128i*)shift_mask_right);
        // Máscara utilizada para shiftear a izquierda los caracteres
        shift_mask_left_mm.x =  _mm_loadu_si128((__m128i*)shift_mask_left);
     
        height = ((seq2_len + vector_len - 1)/ vector_len); // Cantidad de "franjas" de diagonales
        width = (1 + seq1_len + vector_len - 1); 			// Cantidad de diagonales por franja
        score_matrix_sz = height * width * vector_len; 		// Cantidad de celdas de la matriz

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
            diag1_mm.y = _mm256_loadu_si256((__m256i const*) (score_matrix + offset_y));
            diag2_mm.y = _mm256_loadu_si256((__m256i const*) (score_matrix + offset_y + vector_len));

            for( int j = 2; j < width ; j++){
                int offset_x = j * vector_len;
                // String horizontal ------------------------------------------------------------------
                leer_secuencia_fila(j);
                // Calculo scores de izquierda, arriba y diagonal --------------------------------------------------------------------
                calcular_scores(j);
                // Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
                diag_score_mm.y = _mm256_max_epi16(diag_score_mm.y,up_score_mm.y);
                diag_score_mm.y = _mm256_max_epi16(diag_score_mm.y,left_score_mm.y);

                // Almacenamos el puntaje máximo en la posición correcta de la matriz
                _mm256_storeu_si256((__m256i*)(score_matrix + offset_y + offset_x), diag_score_mm.y);

                if(j>=vector_len){
                    v_aux[j - vector_len] =  _mm256_extract_epi16 (diag_score_mm.y, 0x0);
                }

                // Actualizar las 2 diagonales anteriores para la siguiente iteracion, actualizando diag2 con el valor de la actual
                diag1_mm.y = diag2_mm.y;
                diag2_mm.y = diag_score_mm.y;
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
}
}