section .data:

section .rodata:
debuggingW: DB '%i',10,0
indexView: DB '%hi ,%hi',10,0

section .text:
;Enumero los pasos del algoritmo
;1_Determinar el esquema de penalizacion de gaps y la matriz de sustitucion
;2_Inicializar la matriz de scores
;3_Darle puntaje a cada elemento recorriendo de izq-> der, top-> bottom
;4_Traceback para obtener el mejor alineamiento local
extern malloc 
extern free 
extern new_string
extern new_string
extern new_result
extern destroy_result
extern new_alignment
extern destroy_alignment
global SWLIN

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

extern printf

SWLIN:

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
    sub rsp, 24

    ;preservo los parametros
    mov seq1, rdi
    mov seq2, rsi
    xor rows, rows 
    mov rows, r9
    xor columns, columns
    mov columns, [rbp+16]
    lea alignment_ptr, [rbp+24]
    mov qword[rsp],0             ;[rsp]=best_y=0
    mov qword[rsp+8],0           ;[rsp+8]=best_x=0
    ;para mas comodidad de los calculos
    inc columns             ;columns=columns+1
    inc rows                ;rows=rows+1

    ;inicializar la matriz de puntajes
    mov rsi, rdx
    mov rax, columns
    mul rows                 ;rax = (rows+1)*(columns+1)
    shl rax, 2               ;rax = (rows+1)*(columns+1)*2<--tamano en bytes de los shorts
    mov rdx, rsi             ;restauro rdx
    
    
    mov rdi, rax
    sub rsp, 8
    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score          ;rax=scores[rows+1][columns+1]
    add rsp, 8
    
    mov word [rax],  0
    mov r10, rax              ;r10 = scores**
    mov r9, rdx
    mov rax, columns
    mul rows
    mov rdx, r9                 ;rax = (columns+1)*(rows+1)
    xor rsi, rsi
    xor r9, r9
    inc rsi                  ;rsi=1

    xor r9,r9
    .setFirstRow:
    ;add r9w, r8w
    mov word [r10+rsi*2], r9w
    inc rsi 
    cmp rsi, columns
    jl .setFirstRow

    xor rsi, rsi
    add rsi, columns        ;rsi=columns+1
    xor r9, r9
    .setFirstColumn:
    ;add r9w, r8w
    mov word [r10+rsi*2], r9w
    add rsi, columns
    cmp rsi, rax
    jl .setFirstColumn  
    
    ;doble ciclo
    mov rax, columns                
    mov r11, rdx
    mul rows
    shl rax, 1
    mov qword[rsp+16], rax
    mov rdx, r11                  ;[rsp+16] = columns*rows*sizeOf(short)=final
    shl columns, 1                ;columns=columns*2 (para una mejor comparacíon)
    xor rax, rax
    xor rsi, rsi                  ;rsi=0
    .rows:
    xor rdi, rdi                  ;rdi=0 reseteado  
    add rsi, columns              ;rsi=indice_actual_fila*sizeOf(short)
    cmp rsi, [rsp+16]
    jge .EndCicle
    .cols:
        add rdi,2
                         ;rdi=indice_actual_columna*sizeOf(short)
        cmp rdi, columns
        je .rows
        ;score_diag
        sub rsi, columns        ;rsi=indice_anterior_fila*columns*2
        sub rdi, 2              ;rdi=indice_anterior_columna            
        add rsi, rdi
        mov r9w, [r10+rsi]
        sub rsi, rdi

        shr rdi, 1              ; rdi=rdi/2=indice_columna_anterior
        cmp rsi, 0
        je .esCero     
        mov r11, rdx
        xor rdx, rdx
        mov rax, rsi
        div columns
        mov rdx, r11            ;rax=rsi/columns=indice_anterior_fila
        
        .esCero:
        mov r11b, [seq1+rax]     ;seq1[indice_fila_anterior*2]
        cmp r11b, [seq2+rdi]     ;seq2[indice_columna_anterior*2]
        je .matchScore
        add r9w, cx             ;+=missmatch_pen
        jmp .Continue
        .matchScore:
        add r9w, dx             ;+=match_score
        
        .Continue:
        xor r11, r11
        shl rdi, 1
        cmp r9w, 0
        jg .score_left
        mov r9w, 0               ;r9w= best_score
        .score_left:
        add rsi, columns
        add rsi, rdi
        mov r11w, [r10+rsi]     ;scores[columns*indice_fila_actual*2+indice_columna_anterior*2]
        sub rsi, rdi
        add r11w, r8w           ;+=gap_pen                
        cmp r11w, r9w
        jle .score_up
        mov r9w, r11w
        
        .score_up:
        sub rsi, columns
        add rdi, 2
        add rsi, rdi            ;scores[rsi] (lo mismo que lo de abajo)
        mov r11w, [r10+rsi]     ;scores[columns*indice_fila_anterior*2+indice_columna_actual*2]
        add r11w, r8w           ;+=gap_pen
        
        cmp r11w, r9w
        jle .UpdateScoreAndPos
        mov r9w, r11w
        
        
        .UpdateScoreAndPos:
        add rsi, columns
        mov word [r10+rsi], r9w    ;score[y][x]=best_score

        ;update best position
        mov rax, [rsp]
        mov r11, rdx
        mul columns
        mov rdx, r11                ;rax=best_y*2*columns
        mov r11, [rsp+8]            ;r11=best_x*2
        shl r11, 1
        add rax, r11
        mov r9w, [r10+rax]
        cmp r9w, [r10+rsi]
        sub rsi, rdi                    ;vuelvo a poner el valor debido en rsi
        jge .cols       

        mov rax, rsi                    ;rax=rsi/columns=indice_actual*2
        mov r11, rdx        
        div columns     
        mov rdx, r11        
        mov [rsp],   rax                ;best_y=y
        mov [rsp+8], rdi                ;best_x=x
      
        jmp .cols

    .EndCicle:  
    ;generate best_sequences
    shr columns, 1
    mov rdi, columns
    add rdi, rows            
    shl rdi, 1                          ;rdi=(columns+1+rows+1)*2

    
    sub rsp, 8
    push rdi
    push r10
    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score
    pop r10                             ;rax=best_sequence1
    pop rdi
    add rsp, 8
    

    mov r9, [alignment_ptr+result]
    ;CUANDO INTENTA ACCEDER A LA SECUENCIA EN RESULT SEG_FAULT
    mov [r9+sequence_1], rax            ;(alignment->result)->sequence_1= best_sequence1
    ;debbug
    mov rax, r10
    add rsp, 24
    pop alignment_ptr
    pop columns
    pop rows
    pop seq2
    pop seq1
    pop rbp 
    ret
    ;debbug 

    sub rsp, 8
    push rdi
    push r10
    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score
    pop r10                             ;rax=best_sequence1
    pop rdi
    add rsp, 8
    
    mov r9, [alignment_ptr+result]
    mov [r9+sequence_2], rax            ;(alignment->result)->sequence_2= best_sequence2

    ;TENGO QUE REVISAR QUE EFECTIVAMENTE SE HAYAN GUARDADO LAS BEST_Y Y BEST_X
    mov rsi, [rsp]                      ;rsi=best_y
    mov rdi, [rsp+8]                    ;rdi=best_x
    mov word [rsp],0                    ;[rsp]=lenght=0
    
    mov r11, rdx
    mov rax, columns
    mul rsi
    mov rdx,r11                         ;rax=(columns+1)*best_y

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
        shl rsi,1                
        mov r11, rdx
        mov rax, columns        
        mul rsi                 
        mov rdx,r11             ;rax=columns*indice_anterior
        dec rdi                 ;rdi=(rdi-1)*2=x-1            
        shl rdi,1    
        add rax, rdi      
        mov r9w, [r10+rax]
        mov r11, [seq1+rsi]
        cmp r11, [seq2+rdi]
        je .matchScoreTracing
        add r9, missmatch_pen
        jmp .ContinueTracing
        .matchScoreTracing:
        add r9, match_score
        
        .ContinueTracing:
        add rsi, 2
        add rax, 2
        add rax, columns        ;rax=columns*indice_actual
        cmp r9w, [r10+rax]
        jne .score_left_tracing
            sub rsi, 2
            sub rax, 2
            mov r9, [alignment_ptr+result]
            mov r9, [r9+sequence_1]
            mov r11, [seq1+rsi]
            mov [r9+rsi], r11b
            
            mov r11, [alignment_ptr+result]
            mov r11, [r11+sequence_2]
            mov r9, [seq2+rdi]
            mov [r11+rax], r9b
            
            shr rdi, 1
            shr rsi, 1
            dec rdi             ;x--
            dec rsi             ;y--
            add word [rsp],1         ;lenght++
            jmp .tracingBack
        .score_left_tracing:
        mov r11w, word [r10+rax]
        add r11w, r8w
        add rax,2
        cmp r11w, word [r10+rax]
        jne .score_up_tracing
        jmp .bestXPositive
        .score_up_tracing:
        sub rsi, 2
        add rax,2 
        sub rax, columns        ;rax=columns*indice_anterior
        
        mov r11w, word [r10+rax]
        add r11w, r8w
        add rax, columns
        cmp r11w, [r10+rax]
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
            mov [r11+rax*2], r9b
            shr rdi, 1
            shr rsi, 1
            dec rdi                     ;x--
            add word [rsp],1                 ;lenght++
            jmp .tracingBack

        .bestYPositive:
            mov rax, [rsp]
            mov r9, [alignment_ptr+result]
            mov r9, [r9+sequence_2]
            mov byte [r9+rax*2], '_'
            
            mov r11, [alignment_ptr+result]
            mov r11, [r11+sequence_1]
            mov r9, [seq1+rdi-2]
            mov [r11+rax*2], r9b
            dec rsi                     ;y--
            add word [rsp],1                 ;lenght++
            jmp .tracingBack
                

    .reverseSequences:
    mov r9, [rsp]
    shr r9,1                    ;r9=lenght/2
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
    jmp .reversingCicle

    .end:
    xor rsi, rsi
    mov rax, columns
    mul rsi
    ;hago free sobre cada fila
    .freeScoreMatrix:
        mov rdi, [r10+rax*2]
        push rsi
        push rax
        call free
        pop rax
        pop rsi

        inc rsi
        add rax, columns
        cmp rsi, rows
        jl .freeScoreMatrix
    ;hago free sobre el puntero a la matriz
    mov rdi, r10
    call free

    add rsp, 24
    pop alignment_ptr
    pop columns
    pop rows
    pop seq2
    pop seq1
    pop rbp 
    ret
    ;push r11
    ;push r9
    ;sub rsp, 8
    ;push r10
    ;push rdx
    ;push rsi
    ;push rdi
    ;push match_score
    ;push missmatch_pen
    ;push gap_pen
    ;xor rdx, rdx
    ;xor rsi, rsi
    ;mov dx, r9w
    ;mov si, r11w
    ;mov rdi, indexView
    ;call printf 
    ;pop gap_pen
    ;pop missmatch_pen
    ;pop match_score          
    ;pop rdi
    ;pop rsi
    ;pop rdx
    ;pop r10
    ;add rsp, 8
    ;pop r9
    ;pop r11
    