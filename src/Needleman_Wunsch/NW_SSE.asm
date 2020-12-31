global NW_ASM_SSE
extern malloc
extern printf
extern backtracking_c

section .rodata
reverse_mask : DB 0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1
malloc_error_str : db `No se pudo reservar memoria suficiente.\n`,0

; Variables globales
%define constant_gap_xmm xmm0
%define constant_missmatch_xmm xmm1
%define constant_match_xmm xmm2
%define zeroes_xmm xmm3
%define str_row_xmm xmm4
%define str_col_xmm xmm5
%define left_score_xmm xmm6
%define up_score_xmm xmm7
%define diag_score_xmm xmm8
%define reverse_mask_xmm xmm9

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
;   AlignmentMatrix* matrix;
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
%define sequence_offset_sequence 4

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


; AlignmentMatrix offsets
; struct AlignmentMatrix{
;   short* matrix;
; };
%define alignmentmatrix_offset_matrix 0

; Este valor se usa para calcular el tamanio de la matriz
; y poder navegarla. Es necesario actualizarlo si cambia.
%define vector_len 8
%define vector_len_log 3

section .text

; Funciones auxiliares
inicializar_casos_base:
%define offset_y rbx
mov rdi, [rdi + alignment_offset_parameters]
mov di, [rdi + parameters_offset_gap]

; llenamos el vector auxiliar
mov rsi, 0
.loop:
    mov word [v_aux + rsi], -16384 ; SHRT_MIN/2
    inc rsi
    mov rax, width
    dec rax
    cmp rax, rsi
    jne .loop

; inicializar casos base en matriz
mov rsi, 0
.loop1:
    ;offset_y = i * width * vector_len;
    mov rax, rsi
    mul width
    shl rax, vector_len_log
    mov offset_y, rax
   
    mov ax, -16384 ; SHRT_MIN/2
    pinsrw xmm10, eax, 0
    vpbroadcastw xmm10, xmm10
    movdqu [score_matrix + offset_y], xmm10
    pinsrw xmm10, edi, vector_len-1
    movdqu [score_matrix + offset_y + vector_len], xmm10
    
    inc rsi
    mov rax, height
    cmp rsi, rax
    jne .loop1
ret


leer_secuencia_columna:
; rdi = i
%define shift_count xmm10
%define shift_mask xmm11
mov rdx, rdi
inc rdx
shl rdx, vector_len_log ; rdx = (i+1) * vector_len
cmp rdx, seq2_len
jl .else 
sub rdx, seq2_len ; rdx = offset_col
mov rcx, rdi
shl rcx, vector_len_log
sub rcx, rdx
movq str_col_xmm, [seq2 + rcx]

pxor shift_count, shift_count
mov rcx, rdx
shl rcx, 3
pinsrb shift_count, ecx, 0
psrlq str_col_xmm, shift_count

mov ecx, 0xFF
pinsrb shift_mask, ecx, 0
vpbroadcastb shift_mask, shift_mask
mov rcx, 8
sub rcx, rdx
shl rcx, 3
pinsrb shift_count, ecx, 0
psllq shift_mask, shift_count

por str_col_xmm, shift_mask
jmp .end

.else:
movq str_col_xmm, [seq2 + rdi * vector_len]
jmp .end

.end:
punpcklbw str_col_xmm, zeroes_xmm
pshufb str_col_xmm, reverse_mask_xmm
ret


leer_secuencia_fila:
; rdi = j
%define shift_count xmm10
%define shift_mask xmm11

mov rdx, rdi 
sub rdx, vector_len ; rdx = j - vector_len
cmp rdx, 0
jge .elseif ; j-vector_len < 0

mov rcx, vector_len 
sub rcx, rdi ; rcx = offset_str_row
add rdx, rcx
movq str_row_xmm, [seq1 + rcx]

pxor shift_count, shift_count
shl rcx, 3
pinsrb shift_count, ecx, 0
psllq str_row_xmm, shift_count

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
movq str_row_xmm, [seq1 + rdx - vector_len]
pxor shift_count, shift_count
shl rcx, 3
pinsrb shift_count, ecx, 0
psrlq str_row_xmm, shift_count

jmp .end

.else:
movq str_row_xmm, [seq1 + rdi - vector_len]

jmp .end

.end:
punpcklbw str_row_xmm, zeroes_xmm
ret

calcular_scores:
; rdi = j
; rsi = offset_y
; rdx = offset_x


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
;------------------------------------------------------------------

; acceso a las subestructuras de alignment ------------------------
mov rax, [rdi + alignment_offset_sequence_1]
mov seq1, [rax + sequence_offset_sequence]
xor seq1_len, seq1_len
mov r14, [rax + sequence_offset_length]
mov rax, [rdi + alignment_offset_sequence_2]
mov seq2, [rax + sequence_offset_sequence]
xor seq2_len, seq2_len
mov r15, [rax + sequence_offset_length]
;------------------------------------------------------------------

; asignacion de datos en los registros xmm nombrados --------------

; _mm_insert_epi16(constant_gap_xmm,alignment.parameters->gap,0);
; _mm_broadcastw_epi16(constant_gap_xmm);
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_gap]
pinsrw constant_gap_xmm, eax, 0
vpbroadcastw constant_gap_xmm, constant_gap_xmm

; mm_insert_epi16(constant_missmatch_xmm,alignment.parameters->missmatch,0);
; _mm_broadcastw_epi16(constant_missmatch_xmm);
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_missmatch]
pinsrw constant_missmatch_xmm, eax, 0
vpbroadcastw constant_missmatch_xmm, constant_missmatch_xmm

; _mm_insert_epi16(constant_match_xmm,alignment.parameters->match,0);
; _mm_broadcastw_epi16(constant_match_xmm);
mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_match]
pinsrw constant_match_xmm, eax, 0
vpbroadcastw constant_match_xmm, constant_match_xmm

pxor zeroes_xmm,zeroes_xmm
;------------------------------------------------------------------

; Carga de las mascaras -------------------------------------------
movdqu reverse_mask_xmm, [reverse_mask]
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
shr rax, vector_len_log
; Ignoramos rdx porque no podriamos manejar un nro tan grande
; Estaria bueno chequear que no pase
push rdi ; conserva *matrix y deja alineado antes de llamar
mov rdi, rax
shl rdi, 1 ; score_matrix_sz*sizeof(short)
call malloc
pop rdi
cmp rax, 0
je .malloc_error
mov score_matrix,rax

push rdi ; conserva *matrix y deja alineado antes de llamar
mov rdi, width
dec rdi
shl rdi,1
call malloc
pop rdi
cmp rax, 0
je .malloc_error
mov v_aux, rax
;------------------------------------------------------------------

; Casos base ------------------------------------------------------
push rdi
call inicializar_casos_base
pop rdi

; Loop principal --------------------------------------------------
mov rdx, 0 ; i
.loop_i:
    ; Calcular offset_y
    mov rax, rdx
    mul width
    shr rax, vector_len_log
    mov rsi, rax ; rsi = offset_y
    
    push rdi
    push rsi
    push rdx
    push rcx
    mov rdi, rdx ; rdi = i
    call leer_secuencia_columna
    pop rcx
    pop rdx
    pop rsi
    pop rdi
    
    mov rcx, 2 ; j
    .loop_j:
        push rdi
        push rsi
        push rdx
        push rcx
        mov rdi, rcx
        call leer_secuencia_fila
        pop rcx
        pop rdx
        pop rsi
        pop rdi

        push rdi
        push rsi
        push rdx
        push rcx
        mov rdi, rcx ; rdi = j
        mov rdx, rcx 
        shl rdx, vector_len_log ; rdx = offset_x
        call calcular_scores
        pop rcx
        pop rdx
        pop rsi
        pop rdi

        inc rcx
        cmp rcx, width
        jne .loop_j
    inc rdx
    cmp rdx, height
    jne .loop_i
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
mov  rax, 0
sub rsp, 8 ; alineo antes de llamar
call printf
jmp .epilogo