global NW_ASM_AVX
extern malloc
extern free
extern printf
extern backtracking_C
extern new_alignment_matrix
extern get_score_SSE
extern print_registers
extern print_xmm

section .rodata
shift_mask_right : DB 0x70,0x71,0x72,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x7B,0x7C,0x7D,0x7E,0x7F
shift_mask_left : DB 0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF
reverse_mask : DB 0xF,0xE,0xD,0xC,0xB,0xA,0x9,0x8,0x7,0x6,0x5,0x4,0x3,0x2,0x1,0x0
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
call print_registers


movdqu xmm0, [rsp]
add rsp, 16
movdqu xmm1, [rsp]
add rsp, 16
movdqu xmm2, [rsp]
add rsp, 16
movdqu xmm3, [rsp]
add rsp, 16
movdqu xmm4, [rsp]
add rsp, 16
movdqu xmm5, [rsp]
add rsp, 16
movdqu xmm6, [rsp]
add rsp, 16
movdqu xmm7, [rsp]
add rsp, 16
movdqu xmm8, [rsp]
add rsp, 16
movdqu xmm9, [rsp]
add rsp, 16
movdqu xmm10, [rsp]
add rsp, 16
movdqu xmm11, [rsp]
add rsp, 16
movdqu xmm12, [rsp]
add rsp, 16
movdqu xmm13, [rsp]
add rsp, 16
movdqu xmm14, [rsp]
add rsp, 16
movdqu xmm15, [rsp]
add rsp, 16

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

movdqu xmm0, [rsp]
add rsp, 16
movdqu xmm1, [rsp]
add rsp, 16
movdqu xmm2, [rsp]
add rsp, 16
movdqu xmm3, [rsp]
add rsp, 16
movdqu xmm4, [rsp]
add rsp, 16
movdqu xmm5, [rsp]
add rsp, 16
movdqu xmm6, [rsp]
add rsp, 16
movdqu xmm7, [rsp]
add rsp, 16
movdqu xmm8, [rsp]
add rsp, 16
movdqu xmm9, [rsp]
add rsp, 16
movdqu xmm10, [rsp]
add rsp, 16
movdqu xmm11, [rsp]
add rsp, 16
movdqu xmm12, [rsp]
add rsp, 16
movdqu xmm13, [rsp]
add rsp, 16
movdqu xmm14, [rsp]
add rsp, 16
movdqu xmm15, [rsp]
add rsp, 16

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
inicializar_casos_base:
%define offset_y rbx
%define diag_xmm xmm14
%define diag_ymm ymm14
%define temp_xmm xmm15

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
    pinsrw diag_xmm, eax, 0
    vpbroadcastw diag_ymm, diag_xmm
    vmovdqu [score_matrix + 2*offset_y], diag_ymm
    mov rax, rdi
    mul rsi
    shl rax, vector_len_log
    movdqu temp_xmm, diag_xmm
    pinsrw temp_xmm, eax, 7
    vinserti128 diag_ymm, diag_ymm, temp_xmm, 1 
    vmovdqu [score_matrix + 2*offset_y + 2*vector_len], diag_ymm
    
    inc rsi
    mov rax, height
    cmp rsi, rax
    jne .loop1
ret


leer_secuencia_columna:
; rdi = i
%define offset_str_col_xmm xmm15
;es el mismo registro pero lo usamos una vez que dejamos de usar el otro
%define str_col_hi_xmm xmm15

mov rdx, rdi
inc rdx
shl rdx, vector_len_log ; rdx = (i+1) * vector_len
cmp rdx, seq2_len
jl .else 
sub rdx, seq2_len ; rdx = offset_str_col
movdqu str_col_xmm, [seq2 + seq2_len - vector_len]

pinsrb offset_str_col_xmm, edx, 0
vpbroadcastb offset_str_col_xmm, offset_str_col_xmm
paddb offset_str_col_xmm, shift_mask_right_xmm

pshufb str_col_xmm, offset_str_col_xmm
vpblendvb str_col_xmm, str_col_xmm, offset_str_col_xmm, offset_str_col_xmm

jmp .end

.else:
; está accediendo fuera de memoria acá
shl rdi, vector_len_log
movdqu str_col_xmm, [seq2 + rdi]
jmp .end

.end:
pshufb str_col_xmm, reverse_mask_xmm
movdqu str_col_hi_xmm, str_col_xmm
punpckhbw str_col_hi_xmm, zeroes_xmm
punpcklbw str_col_xmm, zeroes_xmm
vinserti128 str_col_ymm, str_col_ymm, str_col_hi_xmm, 1
ret

leer_secuencia_fila:
; rdi = j
%define shift_mask_left_copy_xmm xmm14
%define offset_str_row_xmm xmm15
%define str_row_hi_xmm xmm15

mov rdx, rdi 
sub rdx, vector_len ; rdx = j - vector_len
cmp rdx, 0
jge .elseif ; j-vector_len < 0

    mov rcx, vector_len 
    sub rcx, rdi ; rcx = offset_str_row
    movdqu str_row_xmm, [seq1]
    pinsrb offset_str_row_xmm, ecx, 0
    vpbroadcastb offset_str_row_xmm, offset_str_row_xmm
    movdqu shift_mask_left_copy_xmm, shift_mask_left_xmm
    psubb shift_mask_left_copy_xmm, offset_str_row_xmm
    pshufb str_row_xmm, shift_mask_left_copy_xmm
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
    movdqu str_row_xmm, [seq1 + rdx - vector_len]
    
    pinsrb offset_str_row_xmm, ecx, 0
    vpbroadcastb offset_str_row_xmm, offset_str_row_xmm
    paddb offset_str_row_xmm, shift_mask_right_xmm
    pshufb str_row_xmm, offset_str_row_xmm

jmp .end

.else:
    movdqu str_row_xmm, [seq1 + rdi - vector_len]

jmp .end

.end:

movdqu str_row_hi_xmm, str_row_xmm
punpckhbw str_row_hi_xmm, zeroes_xmm
punpcklbw str_row_xmm, zeroes_xmm
vinserti128 str_row_ymm, str_row_ymm, str_row_hi_xmm, 1
ret

calcular_scores:
; rdi = j
; rsi = offset_y
; rdx = offset_x
    %define cmp_match_ymm ymm14
    %define temp_xmm xmm15
    
    ; left score
    vmovdqu left_score_ymm, diag2_ymm 
    vpaddw left_score_ymm, left_score_ymm, constant_gap_ymm
    
    ; up score 
    vmovdqu up_score_ymm, diag2_ymm
    ; insert shift
    vextracti128 temp_xmm, up_score_ymm, 1
    pextrw esi, temp_xmm, 0b0 
    vpsrldq  up_score_ymm, up_score_ymm, 2
    pinsrw up_score_xmm, esi, 0b0111
    ; insert v_aux
    mov bx, word [v_aux + 2*rdi - 2*1] 
    vextracti128 temp_xmm, up_score_ymm, 1
    pinsrw temp_xmm, ebx, 0b111
    vinserti128 up_score_ymm, up_score_ymm, temp_xmm, 1
    
    vpaddw up_score_ymm, up_score_ymm, constant_gap_ymm
    ;diag score
    vmovdqu diag_score_ymm, diag1_ymm
    ;insert shit
    vextracti128 temp_xmm, diag_score_ymm, 1
    pextrw esi, temp_xmm, 0b0
    vpsrldq  diag_score_ymm, diag_score_ymm, 2
    pinsrw diag_score_xmm, esi, 0b0111
    ; insert v_aux
    mov cx, word [v_aux + 2*rdi - 2*2] 
    vextracti128 temp_xmm, diag_score_ymm, 1
    pinsrw temp_xmm, ecx, 0b111
    vinserti128 diag_score_ymm, diag_score_ymm, temp_xmm, 1
    
    ;compare the 2 strings and put the right penalty (match or missmatch) on each position
    vmovdqu cmp_match_ymm, str_col_ymm
    vpcmpeqw cmp_match_ymm, cmp_match_ymm, str_row_ymm
    vpblendvb cmp_match_ymm, constant_missmatch_ymm, constant_match_ymm, cmp_match_ymm

    ;get the max score of diag, up, left
    vpaddw diag_score_ymm, diag_score_ymm, cmp_match_ymm
    
    
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
; asignacion de datos en los registros xmm nombrados --------------

mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_gap]
pinsrw constant_gap_xmm, eax, 0
vpbroadcastw constant_gap_ymm, constant_gap_xmm

mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_missmatch]
pinsrw constant_missmatch_xmm, eax, 0
vpbroadcastw constant_missmatch_ymm, constant_missmatch_xmm

mov rax, [rdi + alignment_offset_parameters]
mov ax, [rax + parameters_offset_match]
pinsrw constant_match_xmm, eax, 0
vpbroadcastw constant_match_ymm, constant_match_xmm

vpxor zeroes_ymm, zeroes_ymm, zeroes_ymm
;------------------------------------------------------------------

; Carga de las mascaras -------------------------------------------
movdqu reverse_mask_xmm, [reverse_mask]
movdqu shift_mask_right_xmm, [shift_mask_right]
movdqu shift_mask_left_xmm, [shift_mask_left]
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
    shl rax, vector_len_log
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
    vmovdqu diag1_ymm, [score_matrix + 2*rsi]
    vmovdqu diag2_ymm, [score_matrix + 2*rsi + 2*vector_len] 
    
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
        shl rdx, vector_len_log ; rdx = offset_x
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
        vpmaxsw diag_score_ymm, diag_score_ymm, up_score_ymm ; debugging
        vpmaxsw diag_score_ymm, diag_score_ymm, left_score_ymm ; debugging

        ;save the max score in the right position of score matrix
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

    inc rbx
    cmp rbx, height
    jne .loop_i

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



