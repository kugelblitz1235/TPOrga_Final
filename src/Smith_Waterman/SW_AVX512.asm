global SW_ASM_AVX512

extern malloc
extern free
extern printf

extern backtracking_C
extern new_alignment_matrix
extern get_score_SSE

section .rodata
malloc_error_str : db `No se pudo reservar memoria suficiente.\nMalloc: %d\nIntente reservar: %d bytes\n`,0

; Máscara utilizada para invertir el string almacenado en un registro
str_reverse_mask: Dw 0x1F,0x1E,0x1D,0x1C,0x1B,0x1A,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0
; Máscara utilizada para desempaquetar los caracteres a nivel word 
str_512_unpacklo_epi8_mask: DQ 0x0,0xFF,0x1,0xFF,0x2,0xFF,0x3,0xFF
; Mascara para rotar a derecha a nivel word un zmm
score_512_rot_right_word_mask: DW 0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x0

; Registros globales utilizados
%define constant_gap_zmm zmm0
%define constant_gap_ymm ymm0
%define constant_gap_xmm xmm0
%define constant_missmatch_zmm zmm1
%define constant_missmatch_ymm ymm1
%define constant_missmatch_xmm xmm1
%define constant_match_zmm zmm2
%define constant_match_ymm ymm2
%define constant_match_xmm xmm2
%define zeroes_zmm zmm3
%define zeroes_ymm ymm3
%define zeroes_xmm xmm3
%define str_row_zmm zmm4
%define str_row_ymm ymm4
%define str_row_xmm xmm4
%define str_col_zmm zmm5
%define str_col_ymm ymm5
%define str_col_xmm xmm5
%define left_score_zmm zmm6
%define left_score_ymm ymm6
%define left_score_xmm xmm6
%define up_score_zmm zmm7
%define up_score_ymm ymm7
%define up_score_xmm xmm7
%define diag_score_zmm zmm8
%define diag_score_ymm ymm8
%define diag_score_xmm xmm8
%define diag1_zmm zmm9
%define diag1_ymm ymm9
%define diag1_xmm xmm9
%define diag2_zmm zmm10
%define diag2_ymm ymm10
%define diag2_xmm xmm10
%define str_reverse_mask_zmm zmm11
%define str_reverse_mask_ymm ymm11
%define str_reverse_mask_xmm xmm11
%define str_512_unpacklo_epi8_mask_zmm zmm12
%define str_512_unpacklo_epi8_mask_ymm ymm12
%define str_512_unpacklo_epi8_mask_xmm xmm12
%define score_512_rot_right_word_mask_zmm zmm13
%define score_512_rot_right_word_mask_ymm ymm13
%define score_512_rot_right_word_mask_xmm xmm13


; ------- best position and best global score ----------
; | x | y | g | 0 |
%define best_x_y_global xmm14
%define best_x_xmm_pos 0b00
%define best_y_xmm_pos 0b01
%define best_global_xmm_pos 0b10
; -----------------------------------------------------

%define height r8 
%define width r9
%define score_matrix r10
%define v_aux r11
%define seq1 r12
%define seq2 r13
%define seq1_len r14
%define seq2_len r15

; Alignment offsets
; struct Alignment{
;   Sequence* sequence_1;
;   Sequence* sequence_2;
;   Parameters* parameters;
;   Result* result;
;   short* matrix;
; };
%define alignment_offset_sequence_1 0
%define alignment_offset_sequence_2 8
%define alignment_offset_parameters 16
%define alignment_offset_result 24
%define alignment_offset_matrix 32

; Sequence offsets
; struct Sequence{
;   unsigned int length;
;   char* sequence; //data
; };
%define sequence_offset_length 0
%define sequence_offset_sequence 8

; Parameters offsets
; struct Parameters{
;   char* algorithm;
;   short match;
;   short missmatch;
;   short gap;
; };
%define parameters_offset_algorithm 0
%define parameters_offset_match 8
%define parameters_offset_missmatch 10
%define parameters_offset_gap 12

; Result offsets
;struct Result{
;  Sequence* sequence_1;
;  Sequence* sequence_2;
;  short score;
;};
%define result_offset_sequence_1 0
%define result_offset_sequence_2 8
%define result_offset_score 16

; Este valor se usa para calcular el tamanio de la matriz
; y poder navegarla.
; El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
%define vector_len 32
%define vector_len_log 5

section .text

; Funciones auxiliares
inicializar_casos_base:
%define offset_y rbx
%define diag_zmm zmm15
%define diag_ymm ymm15
%define diag_xmm xmm15
%define temp_zmm zmm16
%define temp_ymm ymm16
%define temp_xmm xmm16

mov rdi, [rdi + alignment_offset_parameters]
mov di, [rdi + parameters_offset_gap]

; Inicializamos los valores del vector auxiliar y la matriz de puntajes	
mov rsi, 0
mov rax, width
dec rax
.loop:
    mov word [v_aux + rsi*2], -32768 ; SHRT_MIN
    inc rsi
    cmp rax, rsi
    jne .loop

; Inicializar por cada franja las primeras 2 diagonales. 
; Se pone en cada posicion un cero para no afectar los calculos posteriores			
mov rsi, 0
.loop1:
    mov rax, rsi
    mul width
    shl rax, vector_len_log
    mov offset_y, rax                                           ;offset_y = i * width * vector_len;
   
    vmovdqu16 [score_matrix + 2*offset_y], zeroes_zmm
    vmovdqu16 [score_matrix + 2*offset_y + 2*vector_len], zeroes_zmm
    
    inc rsi
    mov rax, height
    cmp rsi, rax
    jne .loop1
ret

; Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
leer_secuencia_columna:
; rdi = i
%define str_col_temp_zmm zmm15
%define str_col_temp_ymm ymm15
%define str_col_temp_xmm xmm15
;es el mismo registro pero lo usamos una vez que dejamos de usar el otro

%define shift_right_mask k1
%define i_index rdi

    mov rcx, i_index
    inc rcx
    shl rcx, vector_len_log ; rdx = (i+1) * vector_len
    cmp rcx, seq2_len
    jl .else ; Caso de desborde por abajo   
        sub rcx, seq2_len ; rdx = offset_str_col
        ; Shiftear a derecha shift_right_mask una cantidad de posiciones equivalente a caracteres invalidos,
        ; para evitar levantar memoria invalida
        mov edx, 0xFFFFFFFF
        shr edx, cl
        kmovd shift_right_mask, edx
        add rcx, seq2_len
        vpcmpeqw str_col_ymm, str_col_ymm, str_col_ymm
        ; Seleccionar los caracteres validos a derecha con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
        ; A su vez poner en los lugares de posiciones invalidas todos 1s, para evitar que coincida con algun caracter de la secuencia columna
        vmovdqu8 str_col_ymm{shift_right_mask}, [seq2 + rcx - vector_len]

        jmp .end
    .else:; Caso sin desborde
        shl i_index, vector_len_log
        vmovdqu str_col_ymm, [seq2 + i_index]
        ; Desempaquetar los caracteres almacenados en str_col_xmm para trabajar con words

        jmp .end
    .end:
    vpermq str_col_zmm, str_512_unpacklo_epi8_mask_zmm, str_col_zmm
    vpunpcklbw str_col_zmm, str_col_zmm, zeroes_zmm

    ; Invertir el string almacenado en str_col_xmm
    vpermw str_col_zmm, str_reverse_mask_zmm, str_col_zmm
    ret


leer_secuencia_fila:
; rdi = j
%define shift_left_mask k1
%define shift_right_mask k2
%define offset_str_row_zmm zmm15
%define offset_str_row_ymm ymm15
%define offset_str_row_xmm xmm15
%define str_row_hi_zmm zmm16
%define str_row_hi_ymm ymm16
%define str_row_hi_xmm xmm16
%define str_row_lo_zmm zmm17
%define str_row_lo_ymm ymm17
%define str_row_lo_xmm xmm17
%define j_index rdi

    mov rdx, j_index 
    sub rdx, vector_len ; rdx = j - vector_len
    cmp rdx, 0
    jge .elseif ; j-vector_len < 0
        mov rcx, vector_len 
        sub rcx, j_index ; rcx = offset_str_row
        ; Shiftear a izquierda shift_left_mask una cantidad de posiciones equivalente a caracteres invalidos,
        ; para evitar levantar memoria invalida
        mov edx, 0xFFFFFFFF
        shl edx, cl
        kmovd shift_left_mask, edx
        mov rdx, seq1
        sub rdx, rcx
        ; Seleccionar los caracteres validos a izquierda con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
        ; A su vez poner en los lugares de posiciones invalidas todos 0s, para evitar que coincida con algun caracter de la secuencia columna
        vmovdqu8 str_row_ymm{shift_left_mask}{z}, [rdx]
        jmp .end
    .elseif:
    mov rdx, width
    sub rdx, vector_len
    cmp j_index, rdx ; j > width-vector_len
    jle .else
        mov rcx, j_index
        sub rcx, rdx ; rcx = offset_str_row
        ; Shiftear a derecha shift_right_mask una cantidad de posiciones equivalente a caracteres invalidos,
        ; para evitar levantar memoria invalida
        mov edx, 0xFFFFFFFF
        shr edx, cl
        kmovd shift_right_mask, edx
        ; Seleccionar los caracteres validos a derecha con el offset adecuado para que queden cargados correctamente para la comparación mas adelante
        ; A su vez poner en los lugares de posiciones invalidas todos 0s, para evitar que coincida con algun caracter de la secuencia columna     
        vmovdqu8 str_row_ymm{shift_right_mask}{z}, [seq1 + j_index - vector_len]
        jmp .end
    .else:
        vmovdqu str_row_ymm, [seq1 + j_index - vector_len]
        jmp .end
    .end:
    ; Desempaquetamr los caracteres en str_row_ymm para trabajar a nivel word
    vpermq str_row_zmm, str_512_unpacklo_epi8_mask_zmm, str_row_zmm
    vpunpcklbw str_row_zmm, str_row_zmm, zeroes_zmm
    ret
; Calcula los puntajes resultantes de las comparaciones entre caracteres
calcular_scores:
; rdi = j
    %define cmp_match_zmm zmm15
    %define cmp_match_ymm ymm15
    %define cmp_match_xmm xmm15
    %define cmp_mask k1

    ; Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
    vmovdqu16 left_score_zmm, diag2_zmm 
    vpaddsw left_score_zmm, left_score_zmm, constant_gap_zmm
    
    ; Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
    vmovdqu16 up_score_zmm, diag2_zmm
    mov bx, word [v_aux + 2*rdi - 2*1] 
    pinsrw up_score_xmm, ebx, 0b0
    vpermw up_score_zmm, score_512_rot_right_word_mask_zmm, up_score_zmm
    vpaddsw up_score_zmm, up_score_zmm, constant_gap_zmm
    
    ; Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
    ; si coinciden o no los caracteres de la fila y columna correspondientes
    vmovdqu16 diag_score_zmm, diag1_zmm
    mov bx, word [v_aux + 2*rdi - 2*2] 
    pinsrw diag_score_xmm, ebx, 0b0
    vpermw diag_score_zmm, score_512_rot_right_word_mask_zmm, diag_score_zmm
    ; Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
    vpcmpw cmp_mask, str_col_zmm, str_row_zmm, 0
    vpblendmw cmp_match_zmm{cmp_mask}, constant_missmatch_zmm, constant_match_zmm
    vpaddsw diag_score_zmm, diag_score_zmm, cmp_match_zmm
    
    ret

actualizar_posicion_maxima:
; rdi : i
; rsi : j
; best_x_y_global = | x | y | g | - | 
%define nums_zmm zmm15
%define nums_ymm ymm15
%define nums_xmm xmm15
%define nums_s_zmm zmm16
%define nums_s_ymm ymm16
%define nums_s_xmm xmm16
%define index_zmm zmm17
%define index_ymm ymm17
%define index_xmm xmm17
%define index_mask k1
%define max_index eax
; Encontrar el índice del máximo word en el registro zmm
vmovdqu16 nums_zmm, diag_score_zmm                  ; nums_mm = |WWWWWWWW|WWWWWWWW|WWWWWWWW|WWWWWWWW|

vpsrldq nums_s_zmm, nums_zmm, 1*2
vpmaxsw nums_zmm, nums_zmm, nums_s_zmm              ; nums_mm = | W W W W| W W W W| W W W W| W W W W|

vpsrldq nums_s_zmm, nums_zmm, 2*2
vpmaxsw nums_zmm, nums_zmm, nums_s_zmm              ; nums_mm = |   W   W|   W   W|   W   W|   W   W|

vpsrldq nums_s_zmm, nums_zmm, 4*2
vpmaxsw nums_zmm, nums_zmm, nums_s_zmm              ; nums_mm = |       W|       W|       W|       W|

mov bx, 0x1111
kmovw index_mask, ebx
vpcompressd nums_zmm {index_mask}, nums_zmm         ; nums_mm = |        |        |        | W W W W|

vpsrldq nums_s_zmm, nums_zmm, 2*2
vpmaxsw nums_zmm, nums_zmm, nums_s_zmm              ; nums_mm = |        |        |        |   W   W|

vpsrldq nums_s_zmm, nums_zmm, 4*2
vpmaxsw nums_zmm, nums_zmm, nums_s_zmm              ; nums_mm = |        |        |        |       W|

vpbroadcastw nums_zmm, nums_xmm

; Obtener el índice del valor máximo en el registro
vpcmpw index_mask, nums_zmm, diag_score_zmm, 0
kmovd max_index, index_mask
bsf max_index, max_index                            ; Obtener la posicion del word de valor máximo entre todos los de la diagonal

vpextrw edx, nums_xmm, 0b0000
vpextrd ecx, best_x_y_global, best_global_xmm_pos
cmp ecx, edx
jge .menor_a_best

vpinsrd best_x_y_global, edx, best_global_xmm_pos   ; best_global = max_local_score
mov rdx, rdi
shl rdx, vector_len_log
add rdx, vector_len
dec rdx
sub rdx, rax                                        ; rdx = vector_len * i + (vector_len-1) - max_index
vpinsrd best_x_y_global, edx, best_y_xmm_pos
mov rdx, rsi
sub rdx, vector_len
add rdx, rax                                        ; rdx = j - vector_len + max_index
vpinsrd best_x_y_global, edx, best_x_xmm_pos

.menor_a_best:
ret

; Funcion principal (global)
SW_ASM_AVX512:
; struct Alignment{
;   Sequence* sequence_1;
;   Sequence* sequence_2;
;   Parameters* parameters;
;   Result* result;
;   AlignmentMatrix* matrix;
; };

; rdi = *alignment, rsi = debug

; prologo ----------------------------------------------------------
push rbp
mov rbp, rsp

push rbx      ;save current rbx
push r12      ;save current r12
push r13      ;save current r13
push r14      ;save current r14
push r15      ;save current r15

; preservo debug --------------------------------------------------
push rsi

; acceso a las subestructuras de alignment ------------------------
mov rax, [rdi + alignment_offset_sequence_1]
mov seq1, [rax + sequence_offset_sequence]
xor seq1_len, seq1_len
mov r14d, [rax + sequence_offset_length]
mov rax, [rdi + alignment_offset_sequence_2]
mov seq2, [rax + sequence_offset_sequence]
xor seq2_len, seq2_len
mov r15d, [rax + sequence_offset_length]

;------------------------------------------------------------------

; Calculo height, width y score_matrix. Malloc matrix y v_aux -----
mov rax, seq2_len
add rax, vector_len
dec rax
shr rax, vector_len_log
mov height,rax

mov rax, seq1_len
add rax, vector_len
mov width, rax

mov rax, height
mul width
shl rax, vector_len_log


; -----------------------------------------------------------------
; Reservar memoria para la matriz de puntajes y el vector auxiliar, luego inicializamos sus valores
push rdi ; conserva *alignment
push r8
push r9

mov rdi, rax
shl rdi, 1 ; score_matrix_sz*sizeof(short)
sub rsp, 8
call malloc
add rsp, 8
mov rsi, 0
cmp rax, 0
je .malloc_error
pop r9
pop r8

mov score_matrix, rax
push r8
push r9
push r10
mov rdi, width
dec rdi
shl rdi,1
call malloc
mov rsi, 1
cmp rax, 0
je .malloc_error
pop r10
pop r9
pop r8
pop rdi
mov v_aux, rax
;------------------------------------------------------------------
; asignacion de datos en los registros xmm nombrados --------------

; Broadcastear el valor de gap, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_gap]
vpbroadcastw constant_gap_zmm, eax

; Broadcastear el valor de missmatch, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_missmatch]
vpbroadcastw constant_missmatch_zmm, eax

; Broadcastear el valor de match, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_match]
vpbroadcastw constant_match_zmm, eax

; Máscara de ceros
vpxorq zeroes_zmm, zeroes_zmm, zeroes_zmm
;------------------------------------------------------------------

; Carga de las mascaras -------------------------------------------
; Máscara utilizada para invertir el string almacenado en un registro
vmovdqu16 str_reverse_mask_zmm, [str_reverse_mask]
vmovdqu8 str_512_unpacklo_epi8_mask_zmm, [str_512_unpacklo_epi8_mask]
vmovdqu16 score_512_rot_right_word_mask_zmm, [score_512_rot_right_word_mask]
;------------------------------------------------------------------

; Casos base ------------------------------------------------------
push rdi
call inicializar_casos_base

; Loop principal --------------------------------------------------
vpxorq best_x_y_global, best_x_y_global
mov rbx, 0 ; i
.loop_i:
    ; Calcular offset_y
    mov rax, rbx
    mul width
    shl rax, vector_len_log
    mov rsi, rax ; rsi = offset_y
    
    mov rdi, rbx ; rdi = i
    call leer_secuencia_columna
    vmovdqu16 diag1_zmm, [score_matrix + 2*rsi]
    vmovdqu16 diag2_zmm, [score_matrix + 2*rsi + 2*vector_len] 
    
    mov rcx, 2 ; j
    .loop_j:
        push rbx
        push rcx
        mov rdi, rcx
        call leer_secuencia_fila 
        pop rcx
        mov rdx, rcx 
        shl rdx, vector_len_log ; rdx = offset_x
        push rsi
        push rcx
        mov rdi, rcx ; rdi = j
        call calcular_scores 
        pop rcx
        pop rsi
        pop rbx
        
        ; Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba ,diagonalmente y cero
        vpmaxsw diag_score_zmm, diag_score_zmm, up_score_zmm 
        vpmaxsw diag_score_zmm, diag_score_zmm, left_score_zmm
        vpmaxsw diag_score_zmm, diag_score_zmm, zeroes_zmm
        ; Almacenamos el puntaje máximo en la posición correcta de la matriz
        mov rax, rsi
        add rax, rdx
        vmovdqu16 [score_matrix + 2*rax], diag_score_zmm

        cmp rcx, vector_len
        jl .menor
        pextrw eax, diag_score_xmm, 0b0000
        mov [v_aux + 2*rcx - 2*vector_len], ax
        .menor:
        
        push rdx
        push rsi
        push rbx
        push rcx
        mov rdi, rbx ; rdi = i  
        mov rsi, rcx ; rsi = j
        call actualizar_posicion_maxima 
        pop rcx
        pop rbx
        pop rsi
        pop rdx

        vmovdqu16 diag1_zmm, diag2_zmm
        vmovdqu16 diag2_zmm, diag_score_zmm
        inc rcx
        cmp rcx, width
        jne .loop_j    

    inc rbx
    cmp rbx, height
    jne .loop_i

; Restaurar *alignment luego de que el algoritmo termina
pop rdi

.debug:
pop rsi
cmp rsi, 0
je .no_debug
mov [rdi + alignment_offset_matrix], score_matrix


.no_debug:
push rsi
push score_matrix
mov rsi, rdi
mov rdi, score_matrix
mov rdx, vector_len

vpextrd ecx, best_x_y_global, best_x_xmm_pos
vpextrd r8d, best_x_y_global, best_y_xmm_pos
mov r9, 1 ; true
push 1 ; false
push get_score_SSE
call backtracking_C
add rsp, 0x10
pop score_matrix
pop rsi
cmp rsi, 0
jne .epilogo
mov rdi, score_matrix
call free
;------------------------------------------------------------------
; epilogo
.epilogo:
pop r15
pop r14
pop r13
pop r12
pop rbx
pop rbp
ret

.malloc_error:
mov rdi, malloc_error_str
mov rax, 0
call printf
jmp .epilogo



