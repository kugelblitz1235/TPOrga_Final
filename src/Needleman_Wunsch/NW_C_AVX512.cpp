#include "NW_C_AVX512.hpp"

namespace NW{
namespace C{
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
    SIMDreg str_reverse_mask_mm;
    SIMDreg str_512_unpacklo_epi8_mask_mm;
    SIMDreg score_512_rot_right_word_mask_mm;
    SIMDreg diag1_mm, diag2_mm;
    
    //tamaños y punteros a las secuencias
    char *seq1, *seq2;
    unsigned int seq1_len, seq2_len;
    // Máscara utilizada para invertir el string almacenado en un registro
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
            v_aux[i] = SHRT_MIN;
        }

        // Inicializar por cada franja las primeras 2 diagonales. 
        // Se pone en cada posicion un valor negativo grande para no afectar los calculos posteriores
        for(int i = 0 ; i < height ; i++){
            unsigned int offset_y = i * width * vector_len;
            SIMDreg temp_mm;
            SIMDreg diag_mm;
            // Broadcastear el valor SHRT_MIN/2 a nivel word en el registro diag_zmm
            diag_mm.z = _mm512_maskz_set1_epi16(0xFFFFFFFF, SHRT_MIN);
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
            // A su vez poner en los lugares de posiciones invalidas todos 1s, para evitar que coincida con algun caracter de la secuencia columna
            str_col_mm.y = _mm256_mask_loadu_epi8(ones_mm.y, shift_right_mask, (__m256i*)(seq2 + seq2_len - vector_len + offset_str_col));  	// str_row_mm = |1...1|str_col|
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
            // Seleccionar los caracteres validos a izquierda con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
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
        
        // Desempaquetamos los caracteres en str_row_ymm para trabajar a nivel word
        str_row_mm.z = _mm512_permutexvar_epi64(str_512_unpacklo_epi8_mask_mm.z, str_row_mm.z);
        str_row_mm.z = _mm512_unpacklo_epi8(str_row_mm.z, zeroes_mm.z); 
    }

    //Calcula los puntajes resultantes de las comparaciones entre caracteres
    void calcular_scores(int j){
        // Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
        left_score_mm.z = diag2_mm.z;
        left_score_mm.z = _mm512_adds_epi16 (left_score_mm.z, constant_gap_mm.z);
        
        // Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
        up_score_mm.z = diag2_mm.z;
        up_score_mm.x = _mm_insert_epi16(up_score_mm.x,v_aux[j-1],0b0);
        up_score_mm.z = _mm512_permutexvar_epi16(score_512_rot_right_word_mask_mm.z, up_score_mm.z);
        up_score_mm.z = _mm512_adds_epi16(up_score_mm.z, constant_gap_mm.z);
        
        // Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
        // si coinciden o no los caracteres de la fila y columna correspondientes
        diag_score_mm.z = diag1_mm.z;
        diag_score_mm.x = _mm_insert_epi16(diag_score_mm.x,v_aux[j-2],0b0);
        diag_score_mm.z = _mm512_permutexvar_epi16 (score_512_rot_right_word_mask_mm.z, diag_score_mm.z);
        
        // Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
        SIMDreg cmp_match_mm;
        __mmask32 cmp_mask = _mm512_cmpeq_epi16_mask(str_col_mm.z,str_row_mm.z);
        cmp_match_mm.z = _mm512_mask_blend_epi16(cmp_mask, constant_missmatch_mm.z, constant_match_mm.z); 
        diag_score_mm.z = _mm512_adds_epi16(diag_score_mm.z, cmp_match_mm.z);
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
}