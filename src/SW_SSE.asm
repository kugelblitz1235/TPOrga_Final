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
%define match DX
%define missmatch CX
%define match R8W



SW_SSE:

;parametros 
;char* seq1 RDI
;char* seq2 RSI
;int16_t match DX
;int16_t missmatch CX
;int16_t gap R8W
;uint16_t len(seq1) R9W
;uint16_t len(seq2) [RBP+8]
    push rbp 
    mov rbp, rsp
    push seq1
    push seq2
    push rows
    push columns

    ;preservlos parametros
    mov seq1, rdi
    mov seq2, rsi
    xor rows, rows 
    mov rows, r9w
    xor columns, columns
    mov columns, [rbp+8]

    ;inicializar la matriz de puntajes

    

    ;preservo los parametros

    ret
