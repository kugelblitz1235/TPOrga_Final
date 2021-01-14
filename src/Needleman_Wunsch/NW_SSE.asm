global NW_ASM_SSE

extern malloc
extern free
extern printf

extern backtracking_C
extern new_alignment_matrix
extern get_score_SSE

section .rodata
malloc_error_str : db `No se pudo reservar memoria suficiente.\nMalloc: %d\nIntente reservar: %d bytes\n`,0
; Máscara utilizada para invertir el string almacenado en un registro
reverse_mask : DB 0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1

; Variables globales
%define constant_missmatch_xmm xmm1
%define constant_match_xmm xmm2
%define constant_gap_xmm xmm3
%define str_row_xmm xmm4
%define str_col_xmm xmm5
%define left_score_xmm xmm6
%define up_score_xmm xmm7
%define diag_score_xmm xmm8
%define reverse_mask_xmm xmm9
%define zeroes_xmm xmm10
%define diag1_xmm xmm11
%define diag2_xmm xmm12

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

; Este valor se usa para calcular el tamaño de la matriz
; y poder navegarla. Es necesario actualizarlo si cambia.
; El tamaño del vector auxiliar se corresponde con la cantidad de caracteres que vamos a procesar simultáneamente
%define vector_len 8
%define vector_len_log 3

section .text

; Funciones auxiliares
; Inicializar los valores del vector auxiliar y la matriz de puntajes
inicializar_casos_base:

%define diag_xmm xmm13

%define offset_y rbx
%define i_index rsi
    mov rdi, [rdi + alignment_offset_parameters]
    mov di, [rdi + parameters_offset_gap]

    ; Llenar el vector auxiliar con el valor SHRT_MIN / 2
    mov i_index, 0
    mov rax, width
    dec rax
    .loop:
        mov word [v_aux + i_index*2], -16384 ; SHRT_MIN/2
        inc i_index
        cmp rax, i_index
        jne .loop

    ; Inicializar casos base en matriz
    mov i_index, 0
    .loop1:
        mov rax, i_index
        mul width
        shl rax, vector_len_log
        mov offset_y, rax                                                   ; offset_y = i * width * vector_len
    
        mov ax, -16384 ; SHRT_MIN/2
        pinsrw diag_xmm, eax, 0                       
        pshuflw diag_xmm, diag_xmm, 0b0                   
        pshufd diag_xmm, diag_xmm, 0b0                                      ; diag_xmm = | -16384 | -16384 | ... | -16384 | -16384 |
        movdqu [score_matrix + 2*offset_y], diag_xmm
        mov rax, rdi
        mul i_index
        shl rax, vector_len_log
        pinsrw diag_xmm, eax, vector_len-1                                  ; diag_xmm = | gap * i * vector_len | -16384 | ... | -16384 | -16384 |
        movdqu [score_matrix + 2*offset_y + 2*vector_len], diag_xmm
            
        inc i_index
        mov rax, height
        cmp i_index, rax
        jne .loop1
    ret

; Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia columna a utilizar en la comparación
leer_secuencia_columna:
; rdi = i
%define shift_count xmm13
%define shift_mask xmm14
%define i_index rdi

    mov rdx, i_index
    inc rdx
    shl rdx, vector_len_log
    cmp rdx, seq2_len
    jl .else; Caso de desborde por abajo                                ; (i+1)*vector_len < seq2_len ?
        sub rdx, seq2_len                                               ; rdx = offset_col = (i+1) * vector_len
        movq str_col_xmm, [seq2 + seq2_len - vector_len]                

        pxor shift_count, shift_count
        mov rcx, rdx
        shl rcx, 3
        pinsrb shift_count, ecx, 0
        psrlq str_col_xmm, shift_count                                  ; str_col_xmm = | 0...0 | str_col |

        mov ecx, 0xFF
        pcmpeqb shift_mask, shift_mask                                  ; shift_mask = | 1...1 |
        mov rcx, 8
        sub rcx, rdx
        shl rcx, 3
        pinsrb shift_count, ecx, 0
        psllq shift_mask, shift_count                                   ; shift_mask = | 1...1 | 0...0 |

        por str_col_xmm, shift_mask                                     ; str_col_xmm = | 1...1 | str_col |
        jmp .end

    .else:; Caso sin desborde
        movq str_col_xmm, [seq2 + i_index * vector_len]
        jmp .end

    .end:
    ; Desempaquetar los caracteres en str_col_xmm para trabajar con words
    punpcklbw str_col_xmm, zeroes_xmm
    ; Invertir la secuencia de caracteres para compararlos correctamente mas adelante
    pshufb str_col_xmm, reverse_mask_xmm
    ret

; Lee de memoria y almacena correctamente en los registros los caracteres de la secuencia fila a utilizar en la comparación
leer_secuencia_fila:
; rdi = j
%define shift_count xmm13
%define shift_mask xmm14

%define j_index rdi

    mov rdx, j_index 
    sub rdx, vector_len ; rdx = j - vector_len
    cmp rdx, 0; Caso de desborde por izquierda
    jge .elseif                                             ; j - vector_len >= 0 ?

        mov rcx, vector_len 
        sub rcx, j_index                                    ; rcx = offset_str_row = vector_len - j
        movq str_row_xmm, [seq1]

        pxor shift_count, shift_count
        shl rcx, 3
        ; Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
        pinsrb shift_count, ecx, 0
        psllq str_row_xmm, shift_count                      ; str_row_xmm = | str_row | 0...0 |
        jmp .end

    .elseif:; Caso de desborde por derecha                  ; j > width - vector_len
        ; Desplazamiento de puntero a derecha y levantar datos de memoria
        mov rdx, width
        sub rdx, vector_len
        cmp j_index, rdx ; j > width-vector_len
        jle .else

        mov rcx, j_index
        sub rcx, rdx                                        ; rcx = offset_str_row = j - (width - vector_len)

        mov rdx, j_index
        sub rdx, rcx
        movq str_row_xmm, [seq1 + rdx - vector_len]
        pxor shift_count, shift_count
        shl rcx, 3
        ; Acomodar los caracteres en str_row_xmm para que queden en las posiciones correctas para la comparación mas adelante
        pinsrb shift_count, ecx, 0
        psrlq str_row_xmm, shift_count                      ; str_row_xmm = | 0...0 | str_row |

        jmp .end

    .else:; Caso sin desborde
        movq str_row_xmm, [seq1 + j_index - vector_len]

        jmp .end

    .end:

    ; Desempaquetar los caracteres en str_row_xmm para trabajar con words
    punpcklbw str_row_xmm, zeroes_xmm
    ret

; Calcula los puntajes resultantes de las comparaciones entre caracteres
calcular_scores:
    ; rdi = j
    ; rsi = offset_y
    ; rdx = offset_x
    %define cmp_match_xmm xmm0
    %define offset_y rsi
    %define offset_x rdx

    mov rcx, rsi
    add rcx, rdx 
    ; Calcular los scores viniendo por izquierda, sumandole a cada posicion la penalidad del gap
    movdqu left_score_xmm, diag2_xmm
    paddw left_score_xmm, constant_gap_xmm                          
    ; Calcular los scores viniendo por arriba, sumandole a cada posicion la penalidad del gap
    movdqu up_score_xmm, diag2_xmm
    psrldq  up_score_xmm, 2                                         ; up_score_xmm = | 0 | up_score |
    mov bx, word [v_aux + 2*rdi - 2*1]
    pinsrw up_score_xmm, ebx, 0b111                                 ; up_score_xmm = | v_aux[j-1] | up_score |
    paddw up_score_xmm, constant_gap_xmm
    ; Calcular los scores viniendo diagonalmente, sumando en cada caso el puntaje de match o missmatch 
    movdqu diag_score_xmm, diag1_xmm
    psrldq  diag_score_xmm, 2                                       ; up_score_xmm = | 0 | diag_score |
    mov cx, word [v_aux + 2*rdi - 2*2]
    pinsrw diag_score_xmm, ecx, 0b111                               ; diag_score_xmm = | v_aux[j-2] | up_score |

    ; Comparar los dos strings y colocar según corresponda el puntaje correcto (match o missmatch) en cada posición
    movdqu cmp_match_xmm, str_col_xmm
    pcmpeqw cmp_match_xmm, str_row_xmm                              ; Mascara con unos en las posiciones donde coinciden los caracteres
    movdqu str_row_xmm, constant_missmatch_xmm
    pblendvb  str_row_xmm, constant_match_xmm                       ; Seleccionar para cada posicion el puntaje correcto basado en la mascara previa

    ; Obtener el máximo puntaje entre venir por la diagonal, por izquierda y por arriba
    paddw diag_score_xmm, str_row_xmm
    ret

; Funcion principal (global)
NW_ASM_SSE:
; void NW_C_SSE (Alignment& alignment, bool debug)

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
; Punteros y tamaños de las secuencias
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
pshuflw constant_gap_xmm, constant_gap_xmm, 0b0
pshufd constant_gap_xmm, constant_gap_xmm, 0b0

; Broadcastear el valor de missmatch, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_missmatch]
pinsrw constant_missmatch_xmm, eax, 0
pshuflw constant_missmatch_xmm, constant_missmatch_xmm, 0b0
pshufd constant_missmatch_xmm, constant_missmatch_xmm, 0b0

; Broadcastear el valor de match, a nivel word, en el registro
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_match]
pinsrw constant_match_xmm, eax, 0
pshuflw constant_match_xmm, constant_match_xmm, 0b0
pshufd constant_match_xmm, constant_match_xmm, 0b0

; Máscara de ceros
pxor zeroes_xmm,zeroes_xmm
;------------------------------------------------------------------

; Carga de las mascaras -------------------------------------------
; Máscara utilizada para invertir el string almacenado en un registro
movdqu reverse_mask_xmm, [reverse_mask]
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
     
    vmovdqu diag1_xmm, [score_matrix + 2*rsi]
    vmovdqu diag2_xmm, [score_matrix + 2*rsi + 2*vector_len] 

    mov rcx, 2 ; j
    push rbx
    .loop_j:
        push rcx
        mov rdi, rcx
        call leer_secuencia_fila
        pop rcx
        
        mov rdx, rcx 
        shl rdx, vector_len_log ; rdx = offset_x
        push rdx
        push rcx
        mov rdi, rcx ; rdi = j
        call calcular_scores 
        pop rcx
        pop rdx
        ; Guardar en cada posicion de la diagonal el maximo entre los puntajes de venir por izquierda, arriba y diagonalmente
        pmaxsw diag_score_xmm, up_score_xmm
        pmaxsw diag_score_xmm, left_score_xmm

        ; Almacenamos el puntaje máximo en la posición correcta de la matriz
        mov rax, rsi
        add rax, rdx
        movdqu [score_matrix + 2*rax], diag_score_xmm
        

        cmp rcx, vector_len
        jl .menor
            pextrw eax, diag_score_xmm, 0b0000
            mov [v_aux + 2*rcx - 2*vector_len], ax
        .menor:

        vmovdqu diag1_xmm, diag2_xmm
        vmovdqu diag2_xmm, diag_score_xmm
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
mov r9, 0           ; false
push 0              ; false
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

; Se utiliza en el caso de que haya un error en el malloc
.malloc_error:
mov rdi, malloc_error_str
mov rax, 0
call printf
jmp .epilogo