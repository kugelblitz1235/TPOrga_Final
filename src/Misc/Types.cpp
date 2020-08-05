#include "Types.hpp"
#include <iostream>
#include <string.h>

using namespace std;

extern "C" RNA_Sequence* new_RNA_Sequence(char *s) {
    RNA_Sequence* str = (RNA_Sequence*) malloc(sizeof(RNA_Sequence));
    char* new_s = (char*) malloc(strlen(s)+2);
    //pongo en primer lugar el guion
    new_s[0]='-';
    //avanzo el puntero
    new_s++;
    strcpy(new_s, s);
    //retrocedo el puntero para que el RNA_Sequence completo tenga el guion
    new_s--;
    str->sequence = new_s;
    str->length = strlen(s);
    //por el guion dado que no lo cuenta el strlen
    str->length++;

    return str;
}

extern "C" RNA_Sequence* new_RNA_Sequence_length(unsigned int length){
    RNA_Sequence* str = (RNA_Sequence*) malloc(sizeof(RNA_Sequence));
    char* new_s = (char*) malloc(sizeof(char)*length);
    //pongo en primer lugar el guion
    new_s[0]='-';
    //avanzo el puntero
    new_s++;
    //retrocedo el puntero para que el RNA_Sequence completo tenga el guion
    new_s--;
    //la secuencia va a apuntar al espacio reservado por el malloc, es un RNA_Sequence vacio
    str->sequence = new_s;
    str->length = length;
    //por el guion dado que no lo cuenta el strlen
    str->length++;

    return str;
}

void destroy_RNA_Sequence(RNA_Sequence* str) {
    if (str != NULL) {
        if (str->sequence != NULL){
            free(str->sequence);
        };
        free(str);
    };
}

void reverse_RNA_Sequence(RNA_Sequence* seq){
	for(unsigned int i = 0;i < seq->length/2;i++){
		char swap = seq->sequence[i];
		seq->sequence[i] = seq->sequence[seq->length-1-i];
		seq->sequence[seq->length-1-i] = swap;
	}
}

Result* new_result(){
    Result* r = (Result*)malloc(sizeof(Result));
    r->sequence_1 = NULL;
    r->sequence_2 = NULL;
    r->score = 0;

    return r;
}

void destroy_result(Result* result){
    if(result->sequence_1 != NULL)
        destroy_RNA_Sequence(result->sequence_1);
    if(result->sequence_2 != NULL)
        destroy_RNA_Sequence(result->sequence_2);
    free(result);
}

Alignment* new_alignment() {
    Parameters* p = (Parameters*) malloc(sizeof(Parameters));
    p->match = 2;
    p->gap = -1;
    p->missmatch = -1;
    p->algorithm = NULL;
    Alignment* a = (Alignment*) malloc(sizeof(Alignment));
    a->sequence_1 = NULL;
    a->sequence_2 = NULL;
    a->parameters = p;
    a->result = new_result();      
    
    return a;
}

void destroy_alignment(Alignment* a) {
    if (a != NULL) {
        if (a->sequence_1 != NULL) {
            destroy_RNA_Sequence(a->sequence_1);
        }
        if (a->sequence_2 != NULL) {
            destroy_RNA_Sequence(a->sequence_2);
        }
        if (a->result != NULL) {
            destroy_result(a->result);
        }
        if (a->parameters->algorithm != NULL) {
            destroy_RNA_Sequence(a->parameters->algorithm);
        }
        free(a->parameters);
        free(a);
    };    
}

AlignmentMatrix* new_alignment_matrix(unsigned int vector_length,unsigned int seq_length_1,unsigned int seq_length_2){
    AlignmentMatrix* alignment_matrix = (AlignmentMatrix*)malloc(sizeof(AlignmentMatrix));
    alignment_matrix->matrix = (short*)malloc((1+seq_length_1+vector_length-1)*((seq_length_2+vector_length-1)/vector_length)*sizeof(short));
    return alignment_matrix;
}