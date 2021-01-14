global NW_ASM_AVX512

extern malloc
extern free
extern printf

extern backtracking_C
extern new_alignment_matrix
extern get_score_SSE

section .rodata
malloc_error_str : db `No se pudo reservar memoria suficiente.\nMalloc: %d\nIntente reservar: %d bytes\n`,0

; Máscara utilizada para shiftear a derecha los caracteres 
str_shift_right_mask: DW 0x7FE0,0x7FE1,0x7FE2,0x7FE3,0x7FE4,0x7FE5,0x7FE6,0x7FE7,0x7FE8,0x7FE9,0x7FEA,0x7FEB,0x7FEC,0x7FED,0x7FEE,0x7FEF,0x7FF0,0x7FF1,0x7FF2,0x7FF3,0x7FF4,0x7FF5,0x7FF6,0x7FF7,0x7FF8,0x7FF9,0x7FFA,0x7FFB,0x7FFC,0x7FFD,0x7FFE,0x7FFF
; Máscara utilizada para shiftear a izquierda los caracteres
str_shift_left_mask: DW 0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F
; Máscara utilizada para invertir el string almacenado en un registro
str_reverse_mask: Dw 0x1F,0x1E,0x1D,0x1C,0x1B,0x1A,0x19,0x18,0x17,0x16,0x15,0x14,0x13,0x12,0x11,0x10,0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0
str_512_unpacklo_epi8_mask: DQ 0x0,0xFF,0x1,0xFF,0x2,0xFF,0x3,0xFF
; Mascara para rotar a derecha a nivel word un zmm
score_512_rot_right_word_mask: DW 0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF,0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,0x18,0x19,0x1A,0x1B,0x1C,0x1D,0x1E,0x1F,0x0

; Variables globales
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
%define str_shift_right_mask_zmm zmm12
%define str_shift_right_mask_ymm ymm12
%define str_shift_right_mask_xmm xmm12
%define str_shift_left_mask_zmm zmm13
%define str_shift_left_mask_ymm ymm13
%define str_shift_left_mask_xmm xmm13
%define str_512_unpacklo_epi8_mask_zmm zmm14
%define str_512_unpacklo_epi8_mask_ymm ymm14
%define str_512_unpacklo_epi8_mask_xmm xmm14
%define score_512_rot_right_word_mask_zmm zmm15
%define score_512_rot_right_word_mask_ymm ymm15
%define score_512_rot_right_word_mask_xmm xmm15

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
; y poder navegarla. Es necesario actualizarlo si cambia.
%define vector_len 32
%define vector_len_log 5

section .text

; Funciones auxiliares
; Inicializar los valores del vector auxiliar y la matriz de puntajes
inicializar_casos_base:
%define offset_y rbx
%define diag_zmm zmm16
%define diag_ymm ymm16
%define diag_xmm xmm16
%define temp_zmm zmm17
%define temp_ymm ymm17
%define temp_xmm xmm17

    mov rdi, [rdi + alignment_offset_parameters]
    mov di, [rdi + parameters_offset_gap]

    ; Llenamos el vector auxiliar
    mov rsi, 0
    mov rax, width
    dec rax
    .loop:
        mov word [v_aux + rsi*2], -16384 ; SHRT_MIN/2
        inc rsi
        cmp rax, rsi
        jne .loop

    ; Inicializar casos base en matriz
    mov rsi, 0
    .loop1:
        ;offset_y = i * width * vector_len;
        mov rax, rsi
        mul width
        shl rax, vector_len_log
        mov offset_y, rax                                                   ; offset_y = i * width * vector_len
    
        mov ax, -16384 ; SHRT_MIN/2
        vpbroadcastw diag_zmm, eax                                          ; diag_xmm = | -16384 | -16384 | ... | -16384 | -16384 |
        vmovdqu16 [score_matrix + 2*offset_y], diag_zmm
        mov rax, rdi
        mul rsi
        shl rax, vector_len_log
        vmovdqu16 temp_xmm, diag_xmm
        vpinsrw temp_xmm, eax, 7
        vinserti64x2 diag_zmm, diag_zmm, temp_xmm, 3                        ; diag_xmm = | gap * i * vector_len | -16384 | ... | -16384 | -16384 |
        vmovdqu16 [score_matrix + 2*offset_y + 2*vector_len], diag_zmm
        inc rsi
        mov rax, height
        cmp rsi, rax
        jne .loop1
    ret

; Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
leer_secuencia_columna:
; rdi = i
%define str_col_temp_zmm zmm16
%define str_col_temp_ymm ymm16
%define str_col_temp_xmm xmm16

%define shift_right_mask k1
%define i_index rdi

    mov rdx, i_index
    inc rdx
    shl rdx, vector_len_log ; rdx = (i+1) * vector_len
    cmp rdx, seq2_len
    jl .else ; Caso de desborde por abajo   
        sub rdx, seq2_len ; rdx = offset_str_col
        vmovdqu8 str_col_temp_ymm, [seq2 + seq2_len - vector_len]

        ; Unpack 256 -> 512
        vpermq str_col_temp_zmm, str_512_unpacklo_epi8_mask_zmm, str_col_temp_zmm
        vpunpcklbw str_col_temp_zmm, str_col_temp_zmm, zeroes_zmm

        vpbroadcastw str_col_zmm, edx
        ; Sumar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
        vpaddw str_col_zmm, str_col_zmm, str_shift_right_mask_zmm
        ; La mascara tiene un bit en 1 si el caracter asociado a esa posicion es valido
        vpcmpw shift_right_mask, zeroes_zmm, str_col_zmm, 2
        ; Shiftear utilizando la mascara offset_str_col_xmm 
        ; Todos los elementos que sean basura van a tener el valor de la mascara, haciendo que nunca matcheen mas adelante
        vpermw str_col_zmm{shift_right_mask}, str_col_zmm, str_col_temp_zmm
        jmp .end
    .else:; Caso sin desborde
        shl i_index, vector_len_log
        vmovdqu str_col_ymm, [seq2 + i_index]
        ; Desempaquetar los caracteres almacenados en str_col_xmm para trabajar con words
        vpermq str_col_zmm, str_512_unpacklo_epi8_mask_zmm, str_col_zmm
        vpunpcklbw str_col_zmm, str_col_zmm, zeroes_zmm

        jmp .end
    .end:
    ; Invertir el string almacenado en str_col_xmm
    vpermw str_col_zmm, str_reverse_mask_zmm, str_col_zmm
    ret

leer_secuencia_fila:
; rdi = j
%define shift_left_mask k1
%define shift_right_mask k2
%define offset_str_row_zmm zmm16
%define offset_str_row_ymm ymm16
%define offset_str_row_xmm xmm16
%define str_row_hi_zmm zmm17
%define str_row_hi_ymm ymm17
%define str_row_hi_xmm xmm17
%define str_row_lo_zmm zmm18
%define str_row_lo_ymm ymm18
%define str_row_lo_xmm xmm18
    
    mov rdx, rdi 
    sub rdx, vector_len ; rdx = j - vector_len
    cmp rdx, 0
    jge .elseif ; j-vector_len < 0

        mov rcx, vector_len 
        sub rcx, rdi ; rcx = offset_str_row
        vmovdqu str_row_ymm, [seq1]

        ; Desempaquetamr los caracteres en str_row_ymm para trabajar a nivel word
        vpermq str_row_zmm, str_512_unpacklo_epi8_mask_zmm, str_row_zmm
        vpunpcklbw str_row_zmm, str_row_zmm, zeroes_zmm
        
        vpbroadcastw offset_str_row_zmm, ecx
        ; Restar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la izquierda los caracteres para que queden en las posiciones correctas
        vpsubw offset_str_row_zmm, str_shift_left_mask_zmm, offset_str_row_zmm
        ; La mascara tiene un bit en 1 si el caracter asociado a esa posicion es valido
        vpcmpw shift_left_mask, zeroes_zmm, offset_str_row_zmm, 2
        ; Shiftear utilizando la mascara offset_str_row_zmm 
        ; Todos los elementos que sean basura van a tener el valor 0, haciendo que nunca matcheen mas adelante
        vpermw str_row_zmm{shift_left_mask}{z}, offset_str_row_zmm, str_row_zmm
        jmp .end

    .elseif:
    mov rdx, width
    sub rdx, vector_len
    cmp rdi, rdx ; j > width-vector_len
    jle .else

        mov rcx, rdi
        sub rcx, rdx ; rcx = offset_str_row

        mov rdx, rdi
        sub rdx, rcx
        vmovdqu str_row_ymm, [seq1 + rdx - vector_len]
        ; Desempaquetamr los caracteres en str_row_ymm para trabajar a nivel word
        vpermq str_row_zmm, str_512_unpacklo_epi8_mask_zmm, str_row_zmm
        vpunpcklbw str_row_zmm, str_row_zmm, zeroes_zmm
        
        vpbroadcastw offset_str_row_zmm, ecx
        ; Sumar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
        vpaddw offset_str_row_zmm, str_shift_right_mask_zmm, offset_str_row_zmm
        ; La mascara tiene un bit en 1 si el caracter asociado a esa posicion es valido
        vpcmpw shift_right_mask, zeroes_zmm, offset_str_row_zmm, 2
        ; Shiftear utilizando la mascara offset_str_row_zmm 
        ; Todos los elementos que sean basura van a tener el valor 0, haciendo que nunca matcheen mas adelante
        vpermw str_row_zmm{shift_right_mask}{z}, offset_str_row_zmm, str_row_zmm
        jmp .end

    .else:
        vmovdqu str_row_ymm, [seq1 + rdi - vector_len]
        ; Desempaquetamr los caracteres en str_row_ymm para trabajar a nivel word
        vpermq str_row_zmm, str_512_unpacklo_epi8_mask_zmm, str_row_zmm
        vpunpcklbw str_row_zmm, str_row_zmm, zeroes_zmm
        jmp .end

    .end:
    ret
; Calcula los puntajes resultantes de las comparaciones entre caracteres
calcular_scores:
; rdi = j
; rsi = offset_y
; rdx = offset_x
    %define cmp_match_zmm zmm16
    %define cmp_match_ymm ymm16
    %define cmp_match_xmm xmm16
    %define cmp_mask k1

    ; Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
    vmovdqu16 left_score_zmm, diag2_zmm 
    vpaddw left_score_zmm, left_score_zmm, constant_gap_zmm
    
    ; Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
    vmovdqu16 up_score_zmm, diag2_zmm
    mov bx, word [v_aux + 2*rdi - 2*1] 
    pinsrw up_score_xmm, ebx, 0b0
    vpermw up_score_zmm, score_512_rot_right_word_mask_zmm, up_score_zmm
    vpaddw up_score_zmm, up_score_zmm, constant_gap_zmm
    
    ; Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
    ; si coinciden o no los caracteres de la fila y columna correspondientes
    vmovdqu16 diag_score_zmm, diag1_zmm
    mov bx, word [v_aux + 2*rdi - 2*2] 
    pinsrw diag_score_xmm, ebx, 0b0
    vpermw diag_score_zmm, score_512_rot_right_word_mask_zmm, diag_score_zmm
    ; Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
    vpcmpw cmp_mask, str_col_zmm, str_row_zmm, 0
    vpblendmw cmp_match_zmm{cmp_mask}, constant_missmatch_zmm, constant_match_zmm
    vpaddw diag_score_zmm, diag_score_zmm, cmp_match_zmm
    
    ret

; Funcion principal (global)
NW_ASM_AVX512:
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
vmovdqu16 str_reverse_mask_zmm, [str_reverse_mask]
vmovdqu16 str_shift_right_mask_zmm, [str_shift_right_mask]
vmovdqu16 str_shift_left_mask_zmm, [str_shift_left_mask]
vmovdqu8 str_512_unpacklo_epi8_mask_zmm, [str_512_unpacklo_epi8_mask]
vmovdqu16 score_512_rot_right_word_mask_zmm, [score_512_rot_right_word_mask]
;------------------------------------------------------------------

; Casos base ------------------------------------------------------
push rdi
call inicializar_casos_base

; Loop principal --------------------------------------------------
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
    push rbx
    .loop_j:
        push rcx
        mov rdi, rcx
        call leer_secuencia_fila 
        pop rcx
        mov rdx, rcx 
        shl rdx, vector_len_log ; rdx = offset_x
        push rsi
        push rcx
        sub rsp, 8
        mov rdi, rcx ; rdi = j
        call calcular_scores 
        add rsp, 8
        pop rcx
        pop rsi
        ; Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
        vpmaxsw diag_score_zmm, diag_score_zmm, up_score_zmm    
        vpmaxsw diag_score_zmm, diag_score_zmm, left_score_zmm  

        ; Almacenamos el puntaje máximo en la posición correcta de la matriz
        mov rax, rsi
        add rax, rdx
        vmovdqu16 [score_matrix + 2*rax], diag_score_zmm
        

        cmp rcx, vector_len
        jl .menor
        pextrw eax, diag_score_xmm, 0b0000
        mov [v_aux + 2*rcx - 2*vector_len], ax
        .menor:
        vmovdqu16 diag1_zmm, diag2_zmm
        vmovdqu16 diag2_zmm, diag_score_zmm
        inc rcx
        cmp rcx, width
        jne .loop_j    
    pop rbx
    inc rbx
    cmp rbx, height
    jne .loop_i

; Restaurar *alignment luego de que el algoritmo termina
pop rdi

.debug:;  Utilizar para debuggear los valores en la matriz de puntajes
pop rsi
cmp rsi, 0
je .no_debug
mov [rdi + alignment_offset_matrix], score_matrix


.no_debug:
; Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion inferior derecha
push rsi
push score_matrix
mov rsi, rdi
mov rdi, score_matrix
mov rdx, vector_len

mov rcx, seq1_len
dec rcx
mov r8, seq2_len
dec r8
mov r9, 0 ; false
push 0 ; false
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



