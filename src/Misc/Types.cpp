#include "Types.hpp"
#include <iostream>
#include <string.h>

using namespace std;

extern "C" String* new_string(char *s) {
    String* str = (String*) malloc(sizeof(String));
    char* new_s = (char*) malloc(strlen(s)+2);
    //pongo en primer lugar el guion
    new_s[0]='-';
    //avanzo el puntero
    new_s++;
    strcpy(new_s, s);
    //retrocedo el puntero para que el string completo tenga el guion
    new_s--;
    str->sequence = new_s;
    str->length = strlen(s);
    //por el guion dado que no lo cuenta el strlen
    str->length++;

    return str;
}

void destroy_string(String* str) {
    if (str != NULL) {
        if (str->sequence != NULL){
            free(str->sequence);
        };
        free(str);
    };
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
        destroy_string(result->sequence_1);
    if(result->sequence_2 != NULL)
        destroy_string(result->sequence_2);
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
            destroy_string(a->sequence_1);
        }
        if (a->sequence_2 != NULL) {
            destroy_string(a->sequence_2);
        }
        if (a->result != NULL) {
            destroy_result(a->result);
        }
        if (a->parameters->algorithm != NULL) {
            destroy_string(a->parameters->algorithm);
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