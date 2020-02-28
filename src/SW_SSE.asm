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
    mov alignment_ptr, qword [rbp+16]
    shl rdx, 48
    shr rdx, 48             ;limpio la parte alta del registro rdx
    shl rcx, 48
    shr rcx, 48             ;limpio la parte alta del registro rcx
    mov [rsp+best_y],0             ;best_y=0
    mov [rsp+best_x],0           ;best_x=0
    
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
    mov rdi, 1
    inc rsi
    mov r11, rdx
    mov rax, columns        ;rdi=1
    mul rsi                 ;rsi=indice actual 
    mov rdx,r11             ;rax=columns*indice_actual
    cmp rsi, rows
    jge .EndCicle
        .cols:
        cmp rdi, columns
        jge .rows
        
        ;score_diag
        dec rax            ;rax=columns
        shl rdi,2          
        dec rdi            ;rdi=rdi*2-1            
        shl rsi,2          
        dec rsi            ;rsi=rsi*2-1
        mov r9w, word [r10+rax*2+rdi*2]
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
        inc rax
        mov r11w, word [r10+rax*2+rdi]
        add r11w, gap_pen
        cmp r11w, r9w
        jle .score_up
        mov r9w, r11w
        .score_up:
        dec rax
        inc rdi 
        mov r11w, word [r10+rax*2+rdi]
        add r11w, gap_pen
        cmp r11w, r9w
        jle .UpdateScoreAndPos
        mov r9w, r11w
        .UpdateScoreAndPos:
        ;restauro los indices y los avanzo
        inc rsi
        inc rax 
        mov word [r10+rax*2+rdi], r9w    ;score[y][x]=best_score
        
        ;update best position
        mov r9, rax
        mov r11, rdx
        mov rax, [rsp+best_y]
        mul columns
        mov rdx, r11
        mov r11, [rsp+best_x]
        shl r11,2                        ;r11=best_x*2
        mov r11,[r10+rax*2+r11]
        mov rax, r9
        cmp r11, [r10+rax*2+rdi]
        shr rdi,2
        shr rsi,2
        inc rdi
        jge .cols
        mov [rsp+best_y], rsi
        mov [rsp+best_x], rdi           ;update best index
        jmp .cols

    .EndCicle: 
        
    ;generate best_sequences 

    add rsp, 16
    pop alignment_ptr
    pop columns
    pop rows
    pop seq2
    pop seq1
    pop rbp 
    ret
