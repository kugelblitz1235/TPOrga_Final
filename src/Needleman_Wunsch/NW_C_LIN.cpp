#include "NW_C_LIN.hpp"

namespace NW{
namespace C{
namespace LIN{
    void NW(Alignment& alignment, bool debug){
        // Tamaños y punteros a las secuencias
        unsigned int seq1_len = alignment.sequence_1->length;
        unsigned int seq2_len = alignment.sequence_2->length;
        char* seq1 = alignment.sequence_1->sequence;
        char* seq2 = alignment.sequence_2->sequence;
        
        // Matriz de los puntajes que genera el algoritmo NW
        short** scores = (short**)malloc((seq2_len)*sizeof(short*));
        short* tmp_scores = (short*)malloc((seq2_len*seq1_len)*sizeof(short));
        
        // Inicializamos la matriz
        for(unsigned int y = 0;y < seq2_len;y++)
            scores[y] = &tmp_scores[y*seq1_len];
        
        scores[0][0] = 0;
        for(unsigned int y = 1;y < seq2_len;y++){
            scores[y][0] = scores[y-1][0] + alignment.parameters->gap;
        }
        for(unsigned int x = 1;x < seq1_len;x++){
            scores[0][x] = scores[0][x-1] + alignment.parameters->gap;
        }
        
        for(unsigned int y = 1;y < seq2_len;y++){
            for(unsigned int x = 1;x < seq1_len;x++){
                // Moverse por la diagonal, chequeando si se genera un matching
                short score_diag = scores[y-1][x-1] + 
                            alignment.parameters->missmatch*(seq2[y]!= seq1[x]) + 
                            alignment.parameters->match*(seq2[y] == seq1[x]);
                
                // Moverse por la izquierda sumando el puntaje de gap
                short score_left = scores[y][x-1] + alignment.parameters->gap;
                // Moverse por arriba sumando el puntaje de gap
                short score_up = scores[y-1][x] + alignment.parameters->gap;
                // Guardar en la matriz el maximo entre las 3 posibles direcciones
                short best_score = max(score_diag,score_left);
                best_score = max(best_score,score_up);
                scores[y][x] = best_score;
            }
        }
        
        if(debug){// Matriz utilizada para debuggear la matriz de puntajes
            alignment.matrix = (short *) scores;
        }
        
        // Recuperar los 2 strings del mejor alineamiento utilizando backtracking, empezando desde la posicion mas inferior derecha
        backtracking_C(
            (short*)scores,
            alignment,
            1,
            seq1_len-1,seq2_len-1,
            false,
            (score_fun_t)get_score_LIN,
            false
        );

        
        if(!debug){
            free(tmp_scores);
            free(scores);
        }
    }
}
}
}