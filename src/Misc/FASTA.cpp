#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "JSON.hpp"
#include "FASTA.hpp"

using namespace std;

char* parse_FASTA(char* file, Sequence *sequence){
    FILE * fp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;

    fp = fopen(file, "r");
    if (fp == NULL)
        exit(EXIT_FAILURE);
    
    size_t total_length = 0;
     
    while ((read = getline(&line, &len, fp)) != -1) {
        
        if(read>0&&line[read-1]=='\n')
            read--;
        
        if(read > 0 && line[0]!='>')
            total_length += read;
    }

    fseek(fp, 0, SEEK_SET);
    
    char* buffer = (char*)malloc(total_length+1);
    buffer[0] = '-';
    char* buffer_ptr = buffer+1;
    while ((read = getline(&line, &len, fp)) != -1) {
        if(read>0&&line[read-1]=='\n')
            read--;
        
        if(read > 0 && line[0]!='>'){
            strncpy(buffer_ptr, line, read);
            buffer_ptr+=read;
        }
    }
    fclose(fp);
    if (line)
        free(line);

    sequence->sequence = buffer;
    sequence->length = total_length+1;
    
}

void FASTA_to_alignment(Alignment* alignment,char* file1,char* file2){
    
    alignment->sequence_1 = (Sequence*)malloc(sizeof(Sequence));
    alignment->sequence_2 = (Sequence*)malloc(sizeof(Sequence));
    parse_FASTA(file1, alignment->sequence_1);
    parse_FASTA(file2, alignment->sequence_2);

    /*cerr<<(size_t)alignment->sequence_1->sequence<<endl;
    for(int i=0;i<alignment->sequence_1->length;i++)
        cerr<<alignment->sequence_1->sequence[i];
    cerr<<endl;

    cerr<<(size_t)alignment->sequence_2->sequence<<endl;
    for(int i=0;i<alignment->sequence_2->length;i++)
        cerr<<alignment->sequence_2->sequence[i];
    cerr<<endl;*/
}