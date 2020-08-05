
section .rodata:
debuggingW: DB '%c',10,0
indexView: DB '%hi ,%hi',10,0

section .text:
global NWLIN
;podriamos colocar aca unos defines 
%define seq1 RBX
%define seq2 R12
%define rows R13
%define columns R14
%define alignment_ptr R15
%define match_score DX   
%define missmatch_pen CX 
%define gap_pen R8w        
%define sizeOf_short 2
%define best_y 0
%define best_x 8
;Sequence
%define length 0
%define sequence 8
;Parameters
%define algorithm 0
%define match 8
%define missmatch 10
%define gap 12
;Result 
%define result_sequence_1 0
%define result_sequence_2 8
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
extern new_Sequence_from_string
extern new_result
extern destroy_result
extern new_alignment
extern destroy_alignment

NWLIN:
;parametros 
;void* alignment RDI
    push rbp 
    mov rbp, rsp
    sub rsp, 32
    push seq1
    push seq2
    push rows
    push columns
    push alignment_ptr
    
    ;preservo los parametros
    mov alignment_ptr, rdi        
    mov r11, [alignment_ptr+sequence_1]         ;r11=  alignment_ptr*
    mov seq1, [r11+sequence]                    ;seq1= (alignment_ptr->sequence_1)->sequence
    xor rows, rows 
    mov rows, [r11+length]                      ;rows= (alignment_ptr->sequence_1)->length
    mov r11, [alignment_ptr+sequence_2]
    mov seq2, [r11+sequence]                    ;seq2= (alignment_ptr->sequence_2)->sequence
    xor columns, columns
    mov columns, [r11+length]                   ;columns= (alignment_ptr->sequence_2)->length
    mov r11, [alignment_ptr+parameters]
    xor rdx, rdx
    xor rcx, rcx
    xor r8, r8
    mov dx,  word [r11+match]
    mov cx,  word [r11+missmatch]
    mov r8w, word [r11+gap]

    mov qword[rbp-8],0                          ;[rbp-8]=best_y=0
    mov qword[rbp-16],0                         ;[rbp-16]=best_x=0
    
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
    
    mov word [rax],  0        ;scores[0][0]=0
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
    add r9w, gap_pen
    mov word [r10+rsi*2], r9w
    inc rsi 
    cmp rsi, columns
    jl .setFirstRow

    xor rsi, rsi
    add rsi, columns        ;rsi=columns+1
    xor r9, r9
    .setFirstColumn:
    add r9w, gap_pen
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
        inc rax
        inc rdi
        mov r11b, [seq1+rax]     ;seq1[indice_fila_anterior]
        cmp r11b, [seq2+rdi]     ;seq2[indice_columna_anterior]
        je .matchScore
        add r9w, cx             ;+=missmatch_pen
        jmp .score_left
        .matchScore:
        add r9w, dx             ;+=match_score
        
        .score_left:
        dec rdi                 ;rdi-- dado que se incrementó para acceder a la secuencia
        xor r11, r11
        shl rdi, 1
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
        jle .UpdateScore
        mov r9w, r11w
        
        .UpdateScore:
        add rsi, columns
        mov word [r10+rsi], r9w    ;score[y][x]=best_score
        sub rsi, rdi               ;despues de acceder le resto el valor del otro indice que le había sumado
        mov r11, [alignment_ptr+result]
        mov [r11+score], r9w
        
        jmp .cols
    .EndCicle:  
    ;generate best_sequences
    shr columns, 1
    mov r11, columns
    add r11, rows                         ;r11=(columns+1+rows+1)*2
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
    
    mov [rbp-16], rax                      ;[rbp-16]=best_sequence2
    mov word [rbp-24],0                    ;[rbp-24]=lenght=0
    
    mov rsi, rows                          ;rsi=y=sequence1->length=rows 
    dec rsi                                ;porque actualmente rows=rows+1 
    mov rdi, columns                       ;rdi=x=sequence2->length=columns
    dec rdi                                ;porque actualmente columns=columns+1
    
    ;doble ciclo
    mov rax, columns                
    mov r11, rdx
    mul rsi
    shl rax, 1
    mov rsi, rax
    mov rdx, r11                  ;rsi = (columns+1)*sizeOf(short)*y
    shl columns, 1                ;columns=(columns+1)*2 (para una mejor comparacíon)
    shl rdi, 1                    ;rdi= x*2

;asignamos el score a la estructura alignment

    .tracingBack:
        mov r9, rdi
        or  r9, rsi               ;debería chequear que ambos valores son cero
        jz .reverseSequences

        cmp rsi, 0
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
        mov rdx,r11                 ;rax=y-1
        shr rdi, 1                  ;rdi=x-1 
        
        inc rax
        inc rdi
        mov r11b, [seq1+rax]
        cmp r11b, [seq2+rdi]

        je .matchScoreTracing
        add r9w, missmatch_pen
        jmp .ContinueTracing
        .matchScoreTracing:
        add r9w, match_score
        
        .ContinueTracing:
        dec rdi
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
            inc rax
            mov r9, [rbp-8]         ;r9=best_sequence1*
            mov r11b, [seq1+rax]    ;r11b=seq1[y-1]

            xor rax, rax            ;limpio rax
            mov ax, [rbp-24]        ;ax=length (word)
            mov [r9+rax], r11b      ;best_sequence1[length]=r11b
            shr rdi, 1              ;rdi=rdi/2=x-1 
            mov r9, [rbp-16]        ;r9=best_sequence2*
            inc rdi       
            mov r11b, [seq2+rdi]    ;r11b=seq2[x-1]
            dec rdi
            mov [r9+rax], r11b      ;best_sequence2[length]=r11b (rax=length sigue igual)
            shl rdi, 1
            add word [rbp-24],1     ;lenght++
            jmp .tracingBack
        .score_left_tracing:
        sub rax, 2                  ;rax=rax-2=y*2*(columns+1)+x*2-2=y*2*(columns+1)+(x-1)*2
        xor r9, r9
        mov r9w, [r10+rax]          
        add r9w, gap_pen                ;r9w=score[y][x-1]+gap_penalty
        add rax, 2                  ;rax=rax+2=y*2*(columns+1)+(x-1)*2+2=y*2*(columns+1)+x*2
        cmp r9w, [r10+rax]          ;(r9w=>score_left)==scores[x][y]?
        jne .score_up_tracing
        jmp .bestXPositive          ;dado que cuando ocurre el caso en el que x>0 solamente, se ejecuta el mismo procedimiento utilizo el mismo código
        .score_up_tracing:
        sub rax, columns            ;rax=rax-(columns+1)*2=y*2*(columns+1)+x*2-(columns+1)*2=(y-1)*2*(columns+1)+x*2
        xor r9, r9
        mov r9w, [r10+rax]
        add r9w, gap_pen
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
            shr rdi, 1
            xor r9, r9
            mov r9b, [seq2+rdi]             ;r9b=seq2[x-1]
            mov byte [r11+rax], r9b         ;best_sequence2[length]=seq2[x-1]
            shl rdi, 1                      ;rdi=(x-1)*2=> x--
            sub rdi, 2
            add word [rbp-24],1             ;lenght++
            jmp .tracingBack

        .bestYPositive:             ;entra en este caso cuando y>0
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


            sub rsi, columns
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

    mov rdi, r10
    call free
    ;tengo que asignar a la estructura alignment el length y las secuencias resultantes
    mov rdi, [rbp-8]    ;[rbp-8]=best_sequence1
    call new_Sequence_from_string     ;rax=Sequence*(best_sequence1)
    mov r9, [alignment_ptr+result]
    mov rdi, debuggingW
    mov qword [r9+result_sequence_1], rax        ;(result->sequence_1)=new_Sequence_from_string(best_sequence1)
    

    mov rdi, [rbp-16]   ;[rbp-16]=best_sequence2
    call new_Sequence_from_string     ;rax=Sequence*(best_sequence1)
    mov r9, [alignment_ptr+result]
    mov [r9+result_sequence_2], rax      ;(result->sequence_1)=new_Sequence_from_string(best_sequence2)


.skip:
    pop alignment_ptr
    pop columns
    pop rows
    pop seq2
    pop seq1
    add rsp, 32
    pop rbp 
    ret    
