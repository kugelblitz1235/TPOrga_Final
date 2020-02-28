section .data:

section .rodata:

section .text:
;Enumero los pasos del algoritmo
;1_Determinar el esquema de penalizacion de gaps y la matriz de sustitucion
;2_Inicializar la matriz de scores
;3_Darle puntaje a cada elemento recorriendo de izq-> der, top-> bottom
;4_Traceback para obtener el mejor alineamiento local

global SW_SSE

;podriamos colocar aca unos defines 
%define seq1 RBX
%define seq2 R12
%define rows R13
%define columns R14
%define alignment_ptr R15
%define match_score RDX
%define missmatch_pen RCX
%define gap_pen R8
%define sizeOf_short 2
%define best_y 0
%define best_x 8
;Sequence
%define lenght 8
%define sequence 0
;Parameters
%define algorithm 0
%define match 8
%define missmatch 10
%define gap 12
;Result 
%define sequence_1 0
%define sequence_2 8
%define score 16
;Metrics
;Alignment
%define sequence_1 0
%define sequence_2 8
%define parameters 16
%define result 24

;struct Sequence{
;  unsigned int length;-->4
;  char* sequence;-->8
;};
;
;struct Parameters{
;  char* algorithm;
;  short match;
;  short missmatch;
;  short gap;
;};
;
;struct Result{
;  Sequence* sequence_1;
;  Sequence* sequence_2;
;  short score;
;};
;
;struct Metrics{
;  // long long cycles;
;  
;};
;
;struct Alignment{
;  Sequence* sequence_1;
;  Sequence* sequence_2;
;  Parameters* parameters;
;  Result* result;
;};
SW_SSE:

;parametros 
;char* seq1 RDI
;char* seq2 RSI
;int16_t match DX
;int16_t missmatch CX
;int16_t gap R8W
;uint16_t len(seq1) R9W
;uint16_t len(seq2) [RBP+8]
;void* alignment [RBP+16]
    push rbp 
    mov rbp, rsp
    push seq1
    push seq2
    push rows
    push columns
    push alignment_ptr
    sub rsp, 16

    ;preservlos parametros
    mov seq1, rdi
    mov seq2, rsi
    xor rows, rows 
    mov rows, r9
    xor columns, columns
    mov columns, word [rbp+8]
    lea alignment_ptr, qword [rbp+16]
    shl rdx, 48
    shr rdx, 48             ;limpio la parte alta del registro rdx
    shl rcx, 48
    shr rcx, 48             ;limpio la parte alta del registro rcx
    mov [rsp],0             ;best_y=0
    mov [rsp+8],0           ;best_x=0
    
    ;para mas comodidad de los calculos
    inc columns             ;rows=rows+1
    inc rows                ;columns=columns+1

    ;inicializar la matriz de puntajes
    mov rsi, rdx
    mov rax, columns
    mul rows                 ;rax = (rows+1)*(columns+1)
    shl rax, 2               ;rax = (rows+1)*(columns+1)*2<--tamano en bytes de los shorts
    mov rdx, rsi             ;restauro rdx
    
    mov rdi, rax
    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score          ;rax=scores[rows+1][columns+1]
    
    mov word [rax],  0
    ;primera fila y primera columna deben ser seteadas en cero
    mov r10, rax              ;r8 = scores**
    mov r9, rdx
    mov rax, columns
    mul rows
    mov rdx, r9
    sub rax, columns         ;rax = (columns+1)*(rows+1)-columns-1
    xor rsi, rsi
    inc rsi
    .setFirstRow:
    mov r9w, word [r10+rsi*2-1]
    add r9w, gap_pen
    mov word [r10+rsi*2], r9w
    inc rsi
    cmp rsp, columns
    jl .setFirstRow

    
    xor rsi, rsi
    add rsi, columns
    .setFirstColumn:
    mov r9w, word [r10+rsi*2-1]
    add r9w, gap_pen
    mov word [r10+rsi*2], r9w
    add rsi, columns
    cmp rsi, rax
    jl .setFirstColumn

    ;doble ciclo
    xor rsi, rsi
    .rows:
    mov rdi, 1;rdi=1
    inc rsi ;rsi=indice actual 
    cmp rsi, rows
    jge .EndCicle
        .cols:
        cmp rdi, columns
        jge .rows
        
        ;score_diag
        dec rsi                 ;rsi=(rsi-1)*2=y-1
        shl rsi,2          
        mov r11, rdx
        mov rax, columns        
        mul rsi                 
        mov rdx,r11             ;rax=columns*indice_anterior
        dec rdi                 ;rdi=(rdi-1)*2=x-1            
        shl rdi,2          
        mov r9w, word [r10+rax+rdi]
        mov r11, [seq1+rsi]
        cmp r11, [seq2+rdi]
        je .matchScore
        add r9w, missmatch_pen
        jmp .Continue
        .matchScore:
        add r9w, match_score
        
        .Continue:
        cmp r9w, 0
        jg .score_left
        mov r9w, 0          ;r9w= best_score

        .score_left:
        add rsi,2
        add rax, columns        ;rax=columns*indice_actual
        
        mov r11w, word [r10+rax+rdi]
        add r11w, gap_pen
        cmp r11w, r9w
        jle .score_up
        mov r9w, r11w
        .score_up:
        sub rsi, 2
        add rdi,2 
        sub rax, columns        ;rax=columns*indice_anterior
        
        mov r11w, word [r10+rax+rdi]
        add r11w, gap_pen
        cmp r11w, r9w
        jle .UpdateScoreAndPos
        mov r9w, r11w
        .UpdateScoreAndPos:
        shr rsi,2
        inc rsi
        add rax, columns        ;rax=columns*indice_actual
        mov word [r10+rax*2+rdi], r9w    ;score[y][x]=best_score
        
        ;update best position
        mov r11, rdx
        mov r9, [rsp]
        mul columns
        mov rdx, r11
        mov r11, [rsp+8]
        shl r11,2                        ;r11=best_x*2
        mov r11,[r10+r9+r11]
        cmp r11, [r10+rax+rdi]
        shr rdi,2
        inc rdi
        jge .cols
        mov [rsp], rsi
        mov [rsp+8], rdi           ;update best index
        jmp .cols

    .EndCicle: 
        
    ;generate best_sequences
    mov rax, columns
    add rax, rows 
    mov rdi, rax             ;rdi=columns+1+rows+1
    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score          ;rax=best_sequence1
    lea r9, [alignment_ptr+result]
    mov [r9+sequence_1], rax              ;(alignment->result)->sequence_1= best_sequence1

    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score          
    lea r9, [alignment_ptr+result]
    mov [r9+sequence_2], rax              ;(alignment->result)->sequence_2= best_sequence2

    mov rsi, [rsp]           ;rsi=best_y
    mov rdi, [rsp+8]         ;rdi=best_x
    mov [rsp],0              ;[rsp]=lenght

    mov r11, rdx
    mov rax, columns
    mul rsi
    mov rdx,r11          ;rax=(columns+1)*best_y

    ;ambos positivos (distintos de cero)
    .tracingBack:
        mov r9, rdi
        add r9, rsi
        cmp r9, 0
        je .reverseSequences
        
        cmp rsi,0
        je .bestXPositive

        cmp rdi, 0
        je .bestYPositive

        ;ambos son positivos
        ;score_diag
        dec rsi                 ;rsi=(rsi-1)*2=y-1
        shl rsi,2          
        mov r11, rdx
        mov rax, columns        
        mul rsi                 
        mov rdx,r11             ;rax=columns*indice_anterior
        dec rdi                 ;rdi=(rdi-1)*2=x-1            
        shl rdi,2          
        mov r9w, word [r10+rax+rdi]
        mov r11, [seq1+rsi]
        cmp r11, [seq2+rdi]
        je .matchScoreTracing
        add r9w, missmatch_pen
        jmp .ContinueTracing
        .matchScoreTracing:
        add r9w, match_score
        
        .ContinueTracing:
        add rsi, 2
        add rdi, 2
        add rax, columns        ;rax=columns*indice_actual
        cmp r9w, [r10+rax+rdi]
        jne .score_left_tracing
            sub rsi, 2
            sub rdi, 2
            mov r9, [alignment_ptr+result]
            mov r9, [r9+sequence_1]
            mov r11, [seq1+rsi]
            mov byte [r9+rsi], r11
            
            mov r11, [alignment_ptr+result]
            mov r11, [r11+sequence_2]
            mov r9, [seq2+rdi]
            mov byte [r11+rax], r9
            
            shr rdi, 2
            shr rsi, 2
            dec rdi             ;x--
            dec rsi             ;y--
            add [rsp],1         ;lenght++
            jmp .tracingBack
        .score_left_tracing:
        mov r11w, word [r10+rax+rdi]
        add r11w, gap_pen
        add rdi,2
        cmp r11w, word [r10+rax+rdi]
        jne .score_up_tracing
        jmp .bestXPositive
        .score_up_tracing:
        sub rsi, 2
        add rdi,2 
        sub rax, columns        ;rax=columns*indice_anterior
        
        mov r11w, word [r10+rax+rdi]
        add r11w, gap_pen
        add rax, columns
        cmp r11w, [r10+rax+rdi]
        jne .reverseSequences       ;pues el unico caso que queda es que score[y][x]==0
        jmp .bestYPositive
        .bestXPositive:
            mov rax, [rsp]
            mov r9, [alignment_ptr+result]
            mov r9, [r9+sequence_1]
            mov byte [r9+rax*2], '_'
            
            mov r11, [alignment_ptr+result]
            mov r11, [r11+sequence_2]
            mov r9, [seq2+rdi-2]
            mov byte [r11+rax*2], r9
            shr rdi, 2
            shr rsi, 2
            dec rdi                     ;x--
            add [rsp],1                 ;lenght++
            jmp .tracingBack

        .bestYPositive:
            mov rax, [rsp]
            mov r9, [alignment_ptr+result]
            mov r9, [r9+sequence_2]
            mov byte [r9+rax*2], '_'
            
            mov r11, [alignment_ptr+result]
            mov r11, [r11+sequence_1]
            mov r9, [seq1+rdi-2]
            mov byte [r11+rax*2], r9
            dec rsi                     ;y--
            add [rsp],1                 ;lenght++
            jmp .tracingBack
                

    .reverseSequences:
    mov r9, [rsp]
    shr r9,2                    ;r9=lenght/2
    xor rsi, rsi                ;rsi=i=0
    .reversingCicle:
    cmp rsi, r9
    jge .end
    mov rdi, [alignment_ptr+result]
    mov rdi, [rdi+sequence_1]
    mov al, byte[rdi+rsi*2]
    mov byte[rsp+8], al
    mov rax, [rsp]
    dec rax
    sub rax, rsi
    mov  r11b, byte [rdi+rax*2]
    mov  byte[rdi+rsi*2], r11b 
    mov r11b , [rsp+8]
    mov byte [rdi+rax*2], r11b

    mov rdi, [alignment_ptr+result]
    mov rdi, [rdi+sequence_2]
    mov al, byte[rdi+rsi*2]
    mov byte[rsp+8], al
    mov rax, [rsp]
    dec rax
    sub rax, rsi
    mov  r11b, byte [rdi+rax*2]
    mov  byte[rdi+rsi*2], r11b 
    mov r11b , [rsp+8]
    mov byte [rdi+rax*2], r11b
    
    inc rsi
    jmp .reverseCicle

    .end:
    ;FALTA HACER FREE SOBRE TODA LA MATRIZ DE SCORES
    add rsp, 16
    pop alignment_ptr
    pop columns
    pop rows
    pop seq2
    pop seq1
    pop rbp 
    ret
