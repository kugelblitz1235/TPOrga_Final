
section .rodata:
debuggingW: DB '%c',10,0
indexView: DB '%hi ,%hi',10,0

section .text:
;Enumero los pasos del algoritmo
;1_Determinar el esquema de penalizacion de gaps y la matriz de sustitucion
;2_Inicializar la matriz de scores
;3_Darle puntaje a cada elemento recorriendo de izq-> der, top-> bottom
;4_Traceback para obtener el mejor alineamiento local
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
;String
%define lenght 0
%define sequence 8
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
extern malloc 
extern free 
extern new_string
extern new_result
extern destroy_result
extern new_alignment
extern destroy_alignment

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
    sub rsp, 32
    push seq1
    push seq2
    push rows
    push columns
    push alignment_ptr
    
    ;preservo los parametros
    mov seq1, rdi
    mov seq2, rsi
    xor rows, rows 
    mov rows, r9
    xor columns, columns
    mov columns, [rbp+16]
    lea alignment_ptr, [rbp-24]
    mov qword[rbp-8],0             ;[rsp]=best_y=0
    mov qword[rbp-16],0           ;[rsp+8]=best_x=0
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
    push match_score
    push missmatch_pen
    push gap_pen
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score          ;rax=scores[rows+1][columns+1]
    
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
    mov qword[rbp-24], rax
    mov rdx, r11                  ;[rsp+16] = columns*rows*sizeOf(short)=final
    shl columns, 1                ;columns=columns*2 (para una mejor comparacíon)
    xor rax, rax
    xor rsi, rsi                  ;rsi=0
    .rows:
    xor rdi, rdi                  ;rdi=0 reseteado  
    add rsi, columns              ;rsi=indice_actual_fila*sizeOf(short)
    cmp rsi, [rbp-24]
    je .EndCicle
    .cols:
        add rdi,2               ;rdi=indice_actual_columna*sizeOf(short)
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
        mov rax, [rbp-8]              ;rax=best_y
        mov r11, rdx                ;r11=match_score
        mul columns                 ;columns=> (columns+1)*2
        mov rdx, r11                ;rax=best_y*2*(columns+1)
        mov r11, [rbp-16]            ;r11=best_x*2
        shl r11, 1
        add rax, r11
        mov r9w, [r10+rsi] 
        sub rsi, rdi                ;
        cmp r9w, [r10+rax]; actual - best
                           
        jle .cols      
        mov rax, rsi                    ;rax=rsi/columns=indice_actual*2
        mov r11, rdx
        xor rdx, rdx        
        div columns     
        mov rdx, r11
        shr rdi, 1
        mov [rbp-8],   rax               ;best_y=y
        mov [rbp-16], rdi                ;best_x=x

        shl rdi, 1
        jmp .cols
    .EndCicle:  
    ;generate best_sequences
    shr columns, 1
    mov r11, columns
    add r11, rows                         ;r11=(columns+1+rows+1)*2
    mov rsi, [rbp-8]                      ;rsi=best_y
    mov rdi, [rbp-16]                     ;rdi=best_x

;preservo r11 para conservar el valor y a su vez rdi, rsi que tienen el valor de la mejor posicion
    push rsi
    push r11
    push rdi
    push r10
    push match_score
    push missmatch_pen
    push gap_pen
    mov rdi, r11
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score
    pop r10                             
    pop rdi
    pop r11
    pop rsi

    mov [rbp-8], rax                  ;[rbp-8]=best_sequence1
     
    sub rsp, 8
    push rsi
    push rdi
    push r10
    push match_score
    push missmatch_pen
    push gap_pen
    mov rdi, r11
    call malloc 
    pop gap_pen
    pop missmatch_pen
    pop match_score
    pop r10                             
    pop rdi
    pop rsi
    add rsp, 8
    
    mov [rbp-16], rax                       ;[rbp-16]=best_sequence2
    mov word [rbp-24],0                    ;[rbp-24]=lenght=0
    
    mov r11, rdx
    mov rax, columns
    mul rsi
    mov rdx,r11                         ;rax=(columns+1)*best_y

    ;doble ciclo
    mov rax, columns                
    mov r11, rdx
    mul rsi
    shl rax, 1
    mov rsi, rax
    mov rdx, r11                  ;rsi = (columns+1)*sizeOf(short)*best_y
    shl columns, 1                ;columns=(columns+1)*2 (para una mejor comparacíon)
    shl rdi, 1                    ;rdi= best_x*2

    .tracingBack:
        mov r9, rdi
        or  r9, rsi               ;debería chequear que ambos valores son cero
        jz .reverseSequences
        
        cmp rsi,0
        je .bestXPositive

        cmp rdi, 0
        je .bestYPositive

        ;ambos son positivos
        sub rsi, columns            ;rsi=(y-1)*2*(columns+1)
        sub rdi, 2                  ;rdi=(x-1)*2
        add rsi, rdi  
        xor r9, r9    
        mov r9w, [r10+rsi]
        sub rsi, rdi
        ;score_diag 
        mov r11, rdx            
        xor rdx, rdx
        mov rax, rsi        
        div columns 
        mov rdx,r11             ;rax=y-1
        shr rdi, 1              ;rdi=x-1 
        
        mov r11b, [seq1+rax]
        cmp r11b, byte [seq2+rdi]

        je .matchScoreTracing
        add r9w, cx
        jmp .ContinueTracing
        .matchScoreTracing:
        add r9w, dx
        
        .ContinueTracing:
        shl rdi, 1
        add rsi, columns        ;rsi=y*2*(columns+1)
        add rdi, 2              ;rdi=x*2
        mov rax, rsi
        add rax, rdi            ;utilizo rax, porque si utilizara rsi
        cmp r9w, [r10+rax]      ;no podría restaurar su valor
        jne .score_left_tracing
            sub rdi, 2
            sub rsi, columns
            mov r11, rdx        ;entra en este caso cuando el score en la diagonal es el mismo que en la posicion actual de la score matrix            
            xor rdx, rdx        
            mov rax, rsi        
            div columns 
            mov rdx,r11             ;rax=y-1
            mov r9, [rbp-8]         ;r9=best_sequence1*
            mov r11b, [seq1+rax]    ;r11b=seq1[y-1]
            xor rax, rax            ;limpio rax
            mov ax, [rbp-24]        ;ax=length (word)
            mov [r9+rax], r11b      ;best_sequence1[length]=r11b
            shr rdi, 1              ;rdi=rdi/2=x-1 
            mov r9, [rbp-16]        ;r9=best_sequence2*       
            mov r11b, [seq2+rdi]    ;r11b=seq2[x-1]
            mov [r9+rax], r11b      ;best_sequence2[length]=r11b (rax=length sigue igual)
            
            shl rdi, 1
            add word [rbp-24],1     ;lenght++
            jmp .tracingBack
        .score_left_tracing:
        sub rax, 2                  ;rax=rax-2=y*2*(columns+1)+x*2-2=y*2*(columns+1)+(x-1)*2
        xor r9, r9
        mov r9w, [r10+rax]          
        add r9w, r8w                ;r9w=score[y][x-1]+gap_penalty
        add rax, 2                  ;rax=rax+2=y*2*(columns+1)+(x-1)*2+2=y*2*(columns+1)+x*2
        cmp r9w, [r10+rax]          ;(r9w=>score_left)==scores[x][y]?
        jne .score_up_tracing
        jmp .bestXPositive          ;dado que cuando ocurre el caso en el que x>0 solamente, se ejecuta el mismo procedimiento utilizo el mismo código
        .score_up_tracing:
        sub rax, columns            ;rax=rax-(columns+1)*2=y*2*(columns+1)+x*2-(columns+1)*2=(y-1)*2*(columns+1)+x*2
        xor r9, r9
        mov r9w, [r10+rax]
        add r9w, r8w
        add rax, columns
        cmp r9w, [r10+rax]
        jne .reverseSequences       ;pues el unico caso que queda es que score[y][x]==0
        jmp .bestYPositive
        .bestXPositive:             ;entra en este caso cuando x>0
            xor rax, rax
            mov ax, [rbp-24]               ;rax=length
            mov r9, [rbp-8]                 ;r9=best_sequence1*
            mov byte [r9+rax], '_'          ;best_sequence1[length]='_'

            mov r11, [rbp-16]               ;r11=best_sequence2*
            sub rdi, 2
            shr rdi, 1
            xor r9, r9
            mov r9b, [seq2+rdi]             ;r9b=seq2[x-1]
            mov byte [r11+rax], r9b         ;best_sequence2[length]=seq2[x-1]
            shl rdi, 1                      ;rdi=(x-1)*2=> x--
            add word [rbp-24],1             ;lenght++
            jmp .tracingBack

        .bestYPositive:             ;entra en este caso cuando y>0
            sub rsi, columns
            xor rax, rax
            mov ax,  [rbp-24]       ;rax=[0000|0000|0000|length]
            mov r9,  [rbp-16]       ;r9=best_sequence2*
            mov byte [r9+rax], '_'
            
            mov rax, rsi
            mov r9, rdx
            xor rdx, rdx
            div columns
            mov rdx, r9             ;rax=rsi/columns=y*2*(columns+1)/((columns+1)*2)=y
            xor r9, r9 
            mov r9b, [seq1+rax]     ;r9b=seq1[y-1]    
            mov r11, [rbp-8]        ;r11=best_sequence1*
            xor rax, rax
            mov ax,  [rbp-24]       ;rax=[0000|0000|0000|length]
            mov [r11+rax], r9b      ;best_sequence1[length]=seq1[y-1]
            add word [rbp-24],1     ;lenght++
            jmp .tracingBack
                
    .reverseSequences:
    shr columns, 1              ;(columns+1)*2/2=columns+1
    xor r9, r9
    mov r9w, [rbp-24]
    shr r9w,1                    ;r9=lenght/2
    xor rsi, rsi                ;rsi=i=0
    .reversingCicle:
    xor rax, rax
    cmp rsi, r9
    jge .end
    mov rdi, [rbp-8]            ;rdi=best_sequence1*
    mov al,  [rdi+rsi]
    mov byte[rbp-8], al       ;piso el puntero a best_sequence1, porque no tengo mas registros
    mov rax, [rbp-24]       
    dec rax
    sub rax, rsi                ;rax=length-1-i
    mov r11b, [rdi+rax]
    mov byte[rdi+rsi], r11b     ;best_sequence1[i]=best_sequence1[length-1-i] 
    mov r11b , byte[rbp-8]      ;r11=[0000|0000|0000|best_sequence1[i](el valor antes de ser pisado)]
    mov byte [rdi+rax], r11b    ;best_sequence1[length-1-i]=(valor anterior)best_sequence1[i]
    mov qword[rbp-8], rdi       ;vuelvo a poner el best_sequence1* en su lugar

    mov rdi, [rbp-16]           ;rdi=best_sequence2*
    mov al, [rdi+rsi]
    mov byte[rbp-16], al
    mov rax, [rbp-24]
    dec rax
    sub rax, rsi                ;rax=length-1-i
    mov r11b, byte [rdi+rax]   
    mov byte[rdi+rsi], r11b 
    mov r11b , byte[rbp-16]
    mov byte [rdi+rax], r11b
    mov qword [rbp-16], rdi
    
    inc rsi
    jmp .reversingCicle

    .end: 
    ;tengo que poner el caracter terminador a los best_sequence* 
    mov rdi, [rbp-8]                ;rdi=best_sequence1*
    mov rsi, [rbp-24]               
    inc rsi                         ;rsi=length+1
    mov byte [rdi+rsi], 0           ;best_sequence1[length+1]=0
    mov rdi, [rbp-16]               ;rdi=best_sequence2*
    mov byte [rdi+rsi], 0           ;best_sequence2[length+1]=0

    xor rsi, rsi
    mov rax, rows
    mul columns
    ;hago free sobre cada fila
    .freeScoreMatrix:
        cmp rsi, rax
        jge .free
        
        mov rdi, [r10+rsi]
        push rax
        push rsi
        push r10
        call free
        pop r10
        pop rsi
        pop rax
        add rsi, columns
        jmp .freeScoreMatrix

    ;hago free sobre el puntero a la matriz
    .free:
    mov rdi, r10
    call free

    ;tengo que asignar a la estructura alignment el length y las secuencias resultantes
    ;mov rdi, [rbp-8]
    ;call new_string     ;rax=String*(best_sequence1)
    ;mov r9, [alignment_ptr+result]
    ;mov [r9+sequence_1], rax

    ;mov rdi, [rbp-16]
    ;call new_string     ;rax=String*(best_sequence2)         
    ;mov r9, [alignment_ptr+result]
    ;mov [r9+sequence_2], rax    
    .skip:
    pop alignment_ptr
    pop columns
    pop rows
    pop seq2
    pop seq1
    add rsp, 32
    pop rbp 
    ret             ;<--------------TERMINA ACA
;========================DEBUG=========================================

        ;debug
        push r11
        push rax
        push r10
        push rdx
        push rsi
        push rdi
        push match_score
        push missmatch_pen
        push gap_pen
        mov rdx, rax
        mov rdi, indexView
        call printf 
        pop gap_pen
        pop missmatch_pen
        pop match_score          
        pop rdi
        pop rsi
        pop rdx
        pop r10
        pop rax
        pop r11
        ;debug



    
    