global NW_ASM_SSE
extern malloc
extern printf

section .rodata
shift_mask_col : DB 0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F
shift_mask_row : DB 0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF
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
%define shift_mask_col_xmm xmm10
%define shift_mask_row_xmm xmm11


%define height r8 
%define width r9
%define score_matrix r10
%define v_aux r11
%define seq1 r12
%define seq2 r13
%define seq1_len r14d
%define seq2_len r15d

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
%define offset_y rbx
inicializar_casos_base:
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
    pinsrw xmm12, eax, 0
    vpbroadcastw xmm12, xmm12
    movdqu [score_matrix + offset_y], xmm12
    pinsrw xmm12, edi, vector_len-1
    movdqu [score_matrix + offset_y + vector_len], xmm12
    
    inc rsi
    mov rax, height
    cmp rsi, rax
    jne .loop1
ret

leer_secuencia_columna:

leer_secuencia_fila:

calcular_scores:


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
mov seq1_len, [rax + sequence_offset_length]
mov rax, [rdi + alignment_offset_sequence_2]
mov seq2, [rax + sequence_offset_sequence]
mov seq2_len, [rax + sequence_offset_length]
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
movdqu shift_mask_col_xmm, [shift_mask_col]
movdqu shift_mask_row_xmm, [shift_mask_row]
;------------------------------------------------------------------

; Calculo height, width y score_matrix. Malloc matrix y v_aux -----
xor rax, rax
mov eax, seq2_len
add rax, vector_len
dec rax
shr rax, vector_len_log
mov height,rax

xor rax, rax
mov eax, seq1_len
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
push rdi
call inicializar_casos_base
pop rdi
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