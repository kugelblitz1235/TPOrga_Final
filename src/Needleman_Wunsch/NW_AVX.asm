global NW_ASM_AVX

extern malloc
extern free
extern printf

extern backtracking_C
extern new_alignment_matrix
extern get_score_SSE

section .rodata
malloc_error_str : db `No se pudo reservar memoria suficiente.\nMalloc: %d\nIntente reservar: %d bytes\n`,0

; Máscara utilizada para shiftear a derecha los caracteres 
shift_mask_right : DB 0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F
; Máscara utilizada para shiftear a izquierda los caracteres
shift_mask_left : DB 0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF
; Máscara utilizada para invertir el string almacenado en un registro
reverse_mask : DB 0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0

; Variables globales
%define constant_gap_ymm ymm0
%define constant_gap_xmm xmm0
%define constant_missmatch_ymm ymm1
%define constant_missmatch_xmm xmm1
%define constant_match_ymm ymm2
%define constant_match_xmm xmm2
%define zeroes_ymm ymm3
%define zeroes_xmm xmm3
%define str_row_ymm ymm4
%define str_row_xmm xmm4
%define str_col_ymm ymm5
%define str_col_xmm xmm5
%define left_score_ymm ymm6
%define left_score_xmm xmm6
%define up_score_ymm ymm7
%define up_score_xmm xmm7
%define diag_score_ymm ymm8
%define diag_score_xmm xmm8
%define reverse_mask_xmm xmm9
%define shift_mask_right_xmm xmm10
%define shift_mask_left_xmm xmm11
%define diag1_ymm ymm12
%define diag2_ymm ymm13

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
%define vector_len 16
%define vector_len_log 4

section .text

; Funciones auxiliares
; Inicializar los valores del vector auxiliar y la matriz de puntajes
inicializar_casos_base:
%define offset_y rbx
%define diag_xmm xmm14
%define diag_ymm ymm14
%define temp_xmm xmm15

    mov rdi, [rdi + alignment_offset_parameters]
    mov di, [rdi + parameters_offset_gap]

    ; Llenamos el vector auxiliar
    mov rsi, 0
    mov rax, width
    dec rax
    .loop:
        mov word [v_aux + rsi*2], -32768 ; SHRT_MIN/2
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
        mov offset_y, rax                                               ; offset_y = i * width * vector_len
    
        mov ax, -32768 ; SHRT_MIN/2
        pinsrw diag_xmm, eax, 0
        vpbroadcastw diag_ymm, diag_xmm                                 ; diag_xmm = | -16384 | -16384 | ... | -16384 | -16384 |
        vmovdqu [score_matrix + 2*offset_y], diag_ymm
        mov rax, rdi
        mul rsi
        shl rax, vector_len_log
        movdqu temp_xmm, diag_xmm
        pinsrw temp_xmm, eax, 7                                         ; diag_xmm = | gap * i * vector_len | -16384 | ... | -16384 | -16384 |
        vinserti128 diag_ymm, diag_ymm, temp_xmm, 1 
        vmovdqu [score_matrix + 2*offset_y + 2*vector_len], diag_ymm
        
        inc rsi
        mov rax, height
        cmp rsi, rax
        jne .loop1
    ret

; Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
leer_secuencia_columna:
; rdi = i
%define offset_str_col_xmm xmm15
; Es el mismo registro pero lo usamos una vez que dejamos de usar el otro
%define str_col_hi_xmm xmm15

%define i_index rdi

    mov rdx, i_index
    inc rdx
    shl rdx, vector_len_log ; rdx = (i+1) * vector_len
    cmp rdx, seq2_len
    jl .else ; Caso de desborde por abajo                                           ; (i+1)*vector_len < seq2_len ?
    sub rdx, seq2_len                                                               ; rdx = offset_col = (i+1) * vector_len
    movdqu str_col_xmm, [seq2 + seq2_len - vector_len]                              ; str_col_xmm = |0|str_col|

    pinsrb offset_str_col_xmm, edx, 0
    vpbroadcastb offset_str_col_xmm, offset_str_col_xmm                             ; offset_str_col_xmm = |offset|...|offset|
    ; Sumar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
    paddb offset_str_col_xmm, shift_mask_right_xmm                                  ; offset_str_col_xmm = |0x7F + offset|0x7E + offset|...|0x70 + offset|
    ; Shiftear utilizando la mascara offset_str_col_xmm 
    pshufb str_col_xmm, offset_str_col_xmm
    ; Todos los elementos que sean basura van a convertirse en el valor 0xFFFF, haciendo que nunca matcheen mas adelante
    vpblendvb str_col_xmm, str_col_xmm, offset_str_col_xmm, offset_str_col_xmm      ; str_col_xmm = |1...1|str_col|

    jmp .end

    .else:
    ; está accediendo fuera de memoria acá
    shl i_index, vector_len_log
    movdqu str_col_xmm, [seq2 + i_index]
    jmp .end

    .end:
    ; Invertir el string almacenado en str_col_xmm
    pshufb str_col_xmm, reverse_mask_xmm
    movdqu str_col_hi_xmm, str_col_xmm
    ; Desempaquetar los caracteres almacenados en str_col_xmm para trabajar con words
    punpckhbw str_col_hi_xmm, zeroes_xmm
    punpcklbw str_col_xmm, zeroes_xmm
    vinserti128 str_col_ymm, str_col_ymm, str_col_hi_xmm, 1
    ret

; Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
leer_secuencia_fila:
; rdi = j
%define shift_mask_left_copy_xmm xmm14
%define offset_str_row_xmm xmm15
%define str_row_hi_xmm xmm15

%define j_index rdi

    mov rdx, j_index 
    sub rdx, vector_len ; rdx = j - vector_len
    cmp rdx, 0
    jge .elseif ; Caso de desborde por izquierda                    ; j-vector_len < 0 ?

        mov rcx, vector_len 
        sub rcx, j_index ; rcx = offset_str_row
        movdqu str_row_xmm, [seq1]
        ; Broadcastear el offset_str_row en offset_str_row_xmm 
        pinsrb offset_str_row_xmm, ecx, 0
        vpbroadcastb offset_str_row_xmm, offset_str_row_xmm         ; offset_str_col_xmm = |offset|...|offset|
        movdqu shift_mask_left_copy_xmm, shift_mask_left_xmm
        ; Sumar con shift_mask_left_xmm para indicar cuanto hay que shiftear a la izquierda los caracteres para que queden en las posiciones correctas
        psubb shift_mask_left_copy_xmm, offset_str_row_xmm          ; offset_str_col_xmm = |0xF - offset|0xE - offset|...|0x0 - offset|
        ; Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
        pshufb str_row_xmm, shift_mask_left_copy_xmm                ; str_row_xmm = |str_row|0...0|
    jmp .end

    .elseif: ; Caso de desborde por derecha
    mov rdx, width
    sub rdx, vector_len
    cmp j_index, rdx ; j > width-vector_len
    jle .else
        ; Desplazamiento de puntero a derecha y levantar datos de memoria
        mov rcx, j_index
        sub rcx, rdx ; rcx = offset_str_row                         ; Indica cuanto hay que shiftear a la derecha el string fila luego de levantarlo

        mov rdx, j_index
        sub rdx, rcx
        movdqu str_row_xmm, [seq1 + rdx - vector_len]
        
        pinsrb offset_str_row_xmm, ecx, 0
        vpbroadcastb offset_str_row_xmm, offset_str_row_xmm         ; offset_str_col_xmm = |offset|...|offset|
        ; Sumar con shift_mask_right_xmm para indicar cuanto hay que shiftear a la derecha los caracteres para que queden en las posiciones correctas
        paddb offset_str_row_xmm, shift_mask_right_xmm              ; offset_str_col_xmm = |0xF + offset|0xE + offset|...|0x0 + offset|
        ; Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
        pshufb str_row_xmm, offset_str_row_xmm                      ; str_row_xmm = |0...0|str_row|

    jmp .end

    .else:  ; Caso sin desborde
        movdqu str_row_xmm, [seq1 + j_index - vector_len]

    jmp .end

    .end:
    ; Desempaquetar los caracteres en str_row_xmm para trabajar con words
    movdqu str_row_hi_xmm, str_row_xmm
    punpckhbw str_row_hi_xmm, zeroes_xmm
    punpcklbw str_row_xmm, zeroes_xmm
    vinserti128 str_row_ymm, str_row_ymm, str_row_hi_xmm, 1
    ret
; Calcula los puntajes resultantes de las comparaciones entre caracteres
calcular_scores:
; rdi = j
    %define cmp_match_ymm ymm14
    %define temp_xmm xmm15
    
    ; Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
    vmovdqu left_score_ymm, diag2_ymm 
    vpaddsw left_score_ymm, left_score_ymm, constant_gap_ymm
    
    ; Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
    vmovdqu up_score_ymm, diag2_ymm
    ; El shift de a bytes se hace en cada linea de 128 bits,
    ; por lo que hay que mover manualmente el word en la posicion mas baja en la linea de 128 bits
    ; mas alta, e insertarlo en la posicion mas alta de la linea de 128 bits mas baja
    vextracti128 temp_xmm, up_score_ymm, 1                      ; |0xF...0x8|0x7...0x0| -> 0x8
    pextrw esi, temp_xmm, 0b0                                  
    vpsrldq  up_score_ymm, up_score_ymm, 2                      ; |0xF...0x8|0x7...0x0| -> |0x0...0x9|0x0...0x1|
    pinsrw up_score_xmm, esi, 0b0111                           
   
    mov bx, word [v_aux + 2*rdi - 2*1] 
    vextracti128 temp_xmm, up_score_ymm, 1
    pinsrw temp_xmm, ebx, 0b111                                 ; |0x0...0x9|0x0...0x1| -> |0x0...0x9|0x8...0x1|
    vinserti128 up_score_ymm, up_score_ymm, temp_xmm, 1
    vpaddsw up_score_ymm, up_score_ymm, constant_gap_ymm

    ; Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
    ; si coinciden o no los caracteres de la fila y columna correspondientes
    vmovdqu diag_score_ymm, diag1_ymm
    
    vextracti128 temp_xmm, diag_score_ymm, 1                    ; |0xF...0x8|0x7...0x0| -> 0x8
    pextrw esi, temp_xmm, 0b0
    vpsrldq  diag_score_ymm, diag_score_ymm, 2                  ; |0xF...0x8|0x7...0x0| -> |0x0...0x9|0x0...0x1|
    pinsrw diag_score_xmm, esi, 0b0111                          ; |0x0...0x9|0x0...0x1| -> |0x0...0x9|0x8...0x1|
    ; Insert v_aux
    mov cx, word [v_aux + 2*rdi - 2*2] 
    vextracti128 temp_xmm, diag_score_ymm, 1
    pinsrw temp_xmm, ecx, 0b111
    vinserti128 diag_score_ymm, diag_score_ymm, temp_xmm, 1
    
    ; Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
    vmovdqu cmp_match_ymm, str_col_ymm
    vpcmpeqw cmp_match_ymm, cmp_match_ymm, str_row_ymm                                      ; Mascara con unos en las posiciones donde coinciden los caracteres
    vpblendvb cmp_match_ymm, constant_missmatch_ymm, constant_match_ymm, cmp_match_ymm      ; Seleccionar para cada posicion el puntaje correcto basado en la mascara previa

    vpaddsw diag_score_ymm, diag_score_ymm, cmp_match_ymm
    
    ret

; Funcion principal (global)
NW_ASM_AVX:
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
pinsrw constant_gap_xmm, eax, 0
vpbroadcastw constant_gap_ymm, constant_gap_xmm

; Broadcastear el valor de missmatch, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_missmatch]
pinsrw constant_missmatch_xmm, eax, 0
vpbroadcastw constant_missmatch_ymm, constant_missmatch_xmm

; Broadcastear el valor de match, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_match]
pinsrw constant_match_xmm, eax, 0
vpbroadcastw constant_match_ymm, constant_match_xmm

; Máscara de ceros
vpxor zeroes_ymm, zeroes_ymm, zeroes_ymm
;------------------------------------------------------------------

; Carga de las mascaras -------------------------------------------
; Máscara utilizada para invertir el string almacenado en un registro
movdqu reverse_mask_xmm, [reverse_mask]
movdqu shift_mask_right_xmm, [shift_mask_right]
movdqu shift_mask_left_xmm, [shift_mask_left]
;------------------------------------------------------------------

; Casos base ------------------------------------------------------
; Preservar *alignment durante todo el algoritmo
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
    vmovdqu diag1_ymm, [score_matrix + 2*rsi]
    vmovdqu diag2_ymm, [score_matrix + 2*rsi + 2*vector_len] 
    
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
        mov rdi, rcx ; rdi = j
        call calcular_scores 
        pop rcx
        pop rsi
        ; Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
        vpmaxsw diag_score_ymm, diag_score_ymm, up_score_ymm ; debugging
        vpmaxsw diag_score_ymm, diag_score_ymm, left_score_ymm ; debugging
        
        ; Almacenamos el puntaje máximo en la posición correcta de la matriz
        mov rax, rsi
        add rax, rdx
        vmovdqu [score_matrix + 2*rax], diag_score_ymm
        

        cmp rcx, vector_len
        jl .menor
        pextrw eax, diag_score_xmm, 0b0000
        mov [v_aux + 2*rcx - 2*vector_len], ax
        .menor:
        vmovdqu diag1_ymm, diag2_ymm
        vmovdqu diag2_ymm, diag_score_ymm
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
; Traigo debug
pop rsi
cmp rsi, 0
je .no_debug
mov [rdi + alignment_offset_matrix], score_matrix


.no_debug:
; Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion mas inferior derecha
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



