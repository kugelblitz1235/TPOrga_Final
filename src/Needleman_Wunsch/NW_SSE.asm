global NW_ASM_SSE
extern malloc
extern free
extern printf
extern backtracking_C
extern new_alignment_matrix
extern get_score_SSE
extern print_registers
extern print_xmm

section .rodata
reverse_mask : DB 0xE,0xF,0xC,0xD,0xA,0xB,0x8,0x9,0x6,0x7,0x4,0x5,0x2,0x3,0x0,0x1
malloc_error_str : db `No se pudo reservar memoria suficiente.\nMalloc: %d\nIntente reservar: %d bytes\n`,0
dbg_str : db `rsp: 0x%llx-> 0x%llx\n`
%macro dbg_rsp 0
push rdi ; conserva *matrix
push rsi
push rdx
push rcx
push r8
push r9
push r10
push r11

mov rdi, dbg_str
lea rsi, [rsp+0x40]
mov rdx, [rsp+0x40]
mov rax, 0
call printf

pop r11
pop r10
pop r9
pop r8
pop rcx
pop rdx
pop rsi
pop rdi
%endmacro

%macro dbg_print 0
push rax
push rbx
push rdi ; conserva *matrix
push rsi
push rdx
push rcx
push r8
push r9
push r10
push r11
push r12
push r13
push r14
push r15

call print_registers

pop r15
pop r14
pop r13
pop r12
pop r11
pop r10
pop r9
pop r8
pop rcx
pop rdx
pop rsi
pop rdi
pop rbx
pop rax
%endmacro

%macro dbg_print_xmm 0
push rax
push rbx
push rdi ; conserva *matrix
push rsi
push rdx
push rcx
push r8
push r9
push r10
push r11
push r12
push r13
push r14
push r15

movdqu [rsp-16], xmm15
sub rsp, 16
movdqu [rsp-16], xmm14
sub rsp, 16
movdqu [rsp-16], xmm13
sub rsp, 16
movdqu [rsp-16], xmm12
sub rsp, 16
movdqu [rsp-16], xmm11
sub rsp, 16
movdqu [rsp-16], xmm10
sub rsp, 16
movdqu [rsp-16], xmm9
sub rsp, 16
movdqu [rsp-16], xmm8
sub rsp, 16
movdqu [rsp-16], xmm7
sub rsp, 16
movdqu [rsp-16], xmm6
sub rsp, 16
movdqu [rsp-16], xmm5
sub rsp, 16
movdqu [rsp-16], xmm4
sub rsp, 16
movdqu [rsp-16], xmm3
sub rsp, 16
movdqu [rsp-16], xmm2
sub rsp, 16
movdqu [rsp-16], xmm1
sub rsp, 16
movdqu [rsp-16], xmm0
sub rsp, 16
mov rdi, rsp
mov rsi, 16
call print_xmm
add rsp, 16*16

pop r15
pop r14
pop r13
pop r12
pop r11
pop r10
pop r9
pop r8
pop rcx
pop rdx
pop rsi
pop rdi
pop rbx
pop rax
%endmacro

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
mov rax, width
dec rax
.loop:
    mov word [v_aux + rsi*2], -16384 ; SHRT_MIN/2
    inc rsi
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
    movdqu [score_matrix + 2*offset_y], xmm10
    mov rax, rdi
    mul rsi
    shl rax, vector_len_log
    pinsrw xmm10, eax, vector_len-1
    movdqu [score_matrix + 2*offset_y + 2*vector_len], xmm10
        
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
; está accediendo fuera de memoria acá
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
    %define cmp_match_xmm xmm11
    mov rcx, rsi
    add rcx, rdx 
    ; left score
    movdqu left_score_xmm, [score_matrix + 2*rcx - 2*vector_len] 
    paddw left_score_xmm, constant_gap_xmm
    ; up score 
    movdqu up_score_xmm, [score_matrix + 2*rcx - 2*vector_len]
    psrldq  up_score_xmm, 2
    mov bx, word [v_aux + 2*rdi - 2*1] ; mov ebx, dword [v_aux + 2*rdi - 2*1]
    pinsrw up_score_xmm, ebx, 0b111
    paddw up_score_xmm, constant_gap_xmm
    ;diag score
    movdqu diag_score_xmm, [score_matrix + 2*rcx - 2*vector_len]
    psrldq  diag_score_xmm, 2
    mov cx, word [v_aux + 2*rdi - 2*2] ; mov ecx, dword [v_aux + 2*rdi - 2*2]
    pinsrw diag_score_xmm, ecx, 0b111

; %define constant_gap_xmm xmm0
; %define constant_missmatch_xmm xmm1
; %define constant_match_xmm xmm2
; %define zeroes_xmm xmm3
; %define str_row_xmm xmm4
; %define str_col_xmm xmm5
; %define left_score_xmm xmm6
; %define up_score_xmm xmm7
; %define diag_score_xmm xmm8
; %define reverse_mask_xmm xmm9

    dbg_print_xmm
    
    ;compare the 2 strings and put the right penalty (match or missmatch) on each position
    movdqu cmp_match_xmm, str_col_xmm
    pcmpeqw cmp_match_xmm, str_row_xmm
    movdqu str_row_xmm, cmp_match_xmm
    pandn str_row_xmm, constant_missmatch_xmm
    pand cmp_match_xmm, constant_match_xmm

    ;get the max score of diag, up, left
    paddw diag_score_xmm, cmp_match_xmm
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
mov rax, [rdi + alignment_offset_sequence_1]
mov seq1, [rax + sequence_offset_sequence]
xor seq1_len, seq1_len
mov r14d, [rax + sequence_offset_length]
mov rax, [rdi + alignment_offset_sequence_2]
mov seq2, [rax + sequence_offset_sequence]
xor seq2_len, seq2_len
mov r15d, [rax + sequence_offset_length]
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
shl rax, vector_len_log


; -----------------------------------------------------------------
; Ignoramos rdx porque no podriamos manejar un nro tan grande
; Estaria bueno chequear que no pase
push rdi ; conserva *matrix
push rsi
push rdx
push rcx
push r8
push r9
push r10
push r11

mov rdi, rax
shl rdi, 1 ; score_matrix_sz*sizeof(short)
push rdi
sub rsp, 8
call malloc
add rsp, 8
pop rdx
mov rsi, 0
cmp rax, 0
je .malloc_error
pop r11
pop r10
pop r9
pop r8
pop rcx
pop rdx
pop rsi
pop rdi
mov score_matrix,rax
push rdi ; conserva *matrix
push rsi
push rdx
push rcx
push r8
push r9
push r10
push r11
mov rdi, width
dec rdi
shl rdi,1
push rdi
sub rsp, 8
call malloc
add rsp, 8
pop rdx
mov rsi, 1
cmp rax, 0
je .malloc_error
pop r11
pop r10
pop r9
pop r8
pop rcx
pop rdx
pop rsi
pop rdi
mov v_aux, rax
;------------------------------------------------------------------

; Casos base ------------------------------------------------------
push rdi
call inicializar_casos_base
pop rdi

; Loop principal --------------------------------------------------
mov rbx, 0 ; i
.loop_i:
    ; Calcular offset_y
    mov rax, rbx
    mul width
    shr rax, vector_len_log
    mov rsi, rax ; rsi = offset_y
    
    push rdi
    push rsi
    push rbx
    push rcx
    mov rdi, rbx ; rdi = i
    call leer_secuencia_columna
    pop rcx
    pop rbx
    pop rsi
    pop rdi

    mov rcx, 2 ; j
    .loop_j:
        push rdi
        push rsi
        push rbx
        push rcx
        mov rdi, rcx
        call leer_secuencia_fila
        pop rcx
        pop rbx
        pop rsi
        pop rdi
        mov rdx, rcx 
        shl rdx, vector_len_log ; rbx = offset_x
        push rdx
        push rdi
        push rsi
        push rbx
        push rcx
        sub rsp, 8
        mov rdi, rcx ; rdi = j
        call calcular_scores 
        add rsp, 8
        pop rcx
        pop rbx
        pop rsi
        pop rdi
        pop rdx
        pmaxsw diag_score_xmm, up_score_xmm
        pmaxsw diag_score_xmm, left_score_xmm

        ;save the max score in the right position of score matrix
        mov rax, rsi
        add rax, rdx
        movdqu [score_matrix + 2*rax], diag_score_xmm
        

        cmp rcx, vector_len
        jl .menor
        pextrw eax, diag_score_xmm, 0b0000
        mov [v_aux + 2*rcx - 2*vector_len], ax
        .menor:
        inc rcx
        cmp rcx, width

        jmp .debug
        jne .loop_j    
    inc rbx
    cmp rbx, height
    jne .loop_i

; Traigo debug
.debug
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



; old rbp 0x00007fffffffd8f0 rsp 0x7fffffffd608: 0x004055d8