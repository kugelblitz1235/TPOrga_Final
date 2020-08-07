#include <iostream>
#include "NW_C.hpp"

using namespace std;

void NW_C_LIN(Alignment& alignment){
	
	short** scores = (short**)malloc((alignment.sequence_1->length+1)*sizeof(short*));
	
	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++)
		scores[y] = (short*)malloc((alignment.sequence_2->length+1)*sizeof(short));
	
	scores[0][0] = 0;
	
	for(unsigned int y = 1;y < alignment.sequence_1->length+1;y++){
		scores[y][0] = scores[y-1][0] + alignment.parameters->gap;
	}
	
	for(unsigned int x = 1;x < alignment.sequence_2->length+1;x++){
		scores[0][x] = scores[0][x-1] + alignment.parameters->gap;
	}
	
	for(unsigned int y = 1;y < alignment.sequence_1->length+1;y++){
			
		for(unsigned int x = 1;x < alignment.sequence_2->length+1;x++){
			
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(alignment.sequence_1->sequence[y-1]!= alignment.sequence_2->sequence[x-1]) + 
						  alignment.parameters->match*(alignment.sequence_1->sequence[y-1] == alignment.sequence_2->sequence[x-1]);
			
			short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			
			short best_score = max(score_diag,score_left);
			
			best_score = max(best_score,score_up);
			
			scores[y][x] = best_score;
			
		}
	}
	
	cerr << "Score Matrix" << endl;
	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++){
		for(unsigned int x = 0;x < alignment.sequence_2->length+1;x++){
			cerr << (int)scores[y][x] << " ";
		}
		cerr << endl << endl;
	}	
	
	char* best_sequence_1 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	char* best_sequence_2 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	
	unsigned int length = 0;
	
	unsigned int y = alignment.sequence_1->length;
	unsigned int x = alignment.sequence_2->length;
	
	while(y != 0 || x != 0){
		//cerr << x << " " << y << " " << scores[y][x] << endl;
		if(y > 0 && x > 0){
			short score_diag = scores[y-1][x-1] + 
						  alignment.parameters->missmatch*(alignment.sequence_1->sequence[y-1]!= alignment.sequence_2->sequence[x-1]) + 
						  alignment.parameters->match*(alignment.sequence_1->sequence[y-1] == alignment.sequence_2->sequence[x-1]);
			
			//short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = scores[y-1][x] + alignment.parameters->gap;
			
			if(score_diag == scores[y][x]){
				best_sequence_1[length] = alignment.sequence_1->sequence[y-1];
				best_sequence_2[length] = alignment.sequence_2->sequence[x-1];
				x--;
				y--;
			}else if(score_up == scores[y][x]){
				best_sequence_1[length] = alignment.sequence_1->sequence[y-1];
				best_sequence_2[length] = '_';
				y--;
			}else{
				best_sequence_1[length] = '_';
				best_sequence_2[length] = alignment.sequence_2->sequence[x-1];
				x--;
			}
			
		}else if(x > 0){
			best_sequence_1[length] = '_';
			best_sequence_2[length] = alignment.sequence_2->sequence[x-1];
			x--;
		}else if(y > 0){
			best_sequence_1[length] = alignment.sequence_1->sequence[y-1];
			best_sequence_2[length] = '_';
			y--;
		}
		
		length++;
	}
	
	for(unsigned int i = 0;i < length/2;i++){
		char swap = best_sequence_1[i];
		best_sequence_1[i] = best_sequence_1[length-1-i];
		best_sequence_1[length-1-i] = swap;
		swap = best_sequence_2[i];
		best_sequence_2[i] = best_sequence_2[length-1-i];
		best_sequence_2[length-1-i] = swap;
	}
	
	cerr << "Best sequences" << endl;
	for(unsigned int i = 0;i < length;i++){
		cerr << best_sequence_1[i];
	}cerr << endl;
	for(unsigned int i = 0;i < length;i++){
		cerr << best_sequence_2[i];
	}cerr << endl;

	alignment.result->sequence_1 = new_Sequence_from_string(best_sequence_1);
	alignment.result->sequence_2 = new_Sequence_from_string(best_sequence_2);
	alignment.result->score = scores[alignment.sequence_1->length][alignment.sequence_2->length];

	for(unsigned int y = 0;y < alignment.sequence_1->length+1;y++)
		free(scores[y]);
	
	free(scores);
	
	// alignment.length = length;
	// alignment.sequence_1 = best_sequence_1;
	// alignment.sequence_2 = best_sequence_2;


}

Alignment* alignment_by_NW(std::string implementation, char* sequence_1, char* sequence_2, short gap, short missmatch, short match){
	//creo la estructura vacia
	Alignment* alignment = new_alignment();
	
	//relleno la estructura con los valores correspondientes
	alignment->sequence_1 = new_Sequence_from_string(sequence_1);
	alignment->sequence_2 = new_Sequence_from_string(sequence_2);
	(alignment->parameters)->match = match;
	(alignment->parameters)->missmatch = missmatch;
	(alignment->parameters)->gap = gap;

	if(implementation.compare("C") == 0){
		//ejecuto la implementación en c
		NW_C_LIN(*alignment);
	
	}else if(implementation.compare("LIN") == 0){
		
		//ejecuto el algoritmo en asm lineal
		NWLIN(alignment);

	}
	else{
		throw "No existe la implementación ingresada.";
	}

	//devuelvo la estructura modificada
	return alignment;
}



void NW_C_SSE (Alignment& alignment){
	char* seq1 = alignment.sequence_1->sequence;
	char* seq2 = alignment.sequence_2->sequence;
	unsigned int seq1_len = alignment.sequence_1->length;
	unsigned int seq2_len = alignment.sequence_2->length;
	
	//en este caso hardcodeamos el tamaño del vector
	int vector_len = 4;
	
	int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
	int width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
	int score_matrix_sz = height * width * vector_len; // largo de diagonal
		
	short* score_matrix =  (short*)malloc(score_matrix_sz*sizeof(short));
	
	short* v_aux = (short*)malloc((width-1)*sizeof(short));
	//llenamos el vector auxiliar
	for(int i = 0;i < width-1;i++){
		v_aux[i] = SHRT_MIN;
	}

	int count=0;
	for(int i = 0 ; i < height ; i++){
		unsigned int offset_y = i * width * vector_len;
		for( int j = 0; j < 2 ; j++){
			unsigned int offset_x = j * vector_len;
			//emulamos simd
			for( int k = 0;k < vector_len;k++){
				if( j==1 && k == vector_len-1)
					score_matrix[offset_y + offset_x + k] = 0;
				else
					score_matrix[offset_y + offset_x + k] = SHRT_MIN;
				count++;
			}			
		}
	}
	/******************************************************************************************************/
	//arrays auxiliares para el calculo de los scores
	char* str_row = (char*)malloc((vector_len+1) * sizeof(char)); //vector_len+1 es por debuggeo
	char* str_col = (char*)malloc((vector_len+1) * sizeof(char));
	str_row[vector_len]=0;
	str_col[vector_len]=0;

	for( int i = 0 ; i < height ; i++){
		int offset_y = i * width * vector_len;
		for( int j = 2; j < width ; j++){
			int offset_x = j * vector_len;
			//emulamos simd
			
			//cerr<<"j-vector_len :"<<j<<"-"<<vector_len<<"="<<j-vector_len<<"\n";	
			if(j-vector_len < 0){ //desborde por izquierda
				cerr<<"izq"<<endl;
				//simd : desplazamiento de puntero y levantar datos de memoria
				//levantamos el string con el puntero offseteado para no acceder afuera
				
				int offset_str_row = vector_len - j;
			//	cerr<<offset_str_row<<"\n";	
				
				for(int k = 0;k < vector_len;k++){
			//		cerr<<k<<endl;
					str_row[k] = seq1[offset_y + offset_x + k - vector_len + offset_str_row];
				}
				//simd : shift
				//shifteamos los chars para que quede bien (el sentido es contrario a cuando lo hagamos en simd)
				for(int s = 0;s < offset_str_row; s++){
					for(int k = vector_len-1;k > 0;k--){
						str_row[k+1] = str_row[k];
					}
				}
			}else if(j-vector_len >= width-vector_len){ // desborde por derecha
				cerr<<"der"<<endl;
				//simd : desplazamiento de puntero y levantar datos de memoria
				//levantamos el string con el puntero offseteado para no acceder afuera
				int offset_str_row = j - width-vector_len + 1;
			//	cerr<<offset_str_row<<"\n";	

				for(int k = 0;k < vector_len;k++){
			//		cerr<<k<<endl;
					str_row[k] = seq1[offset_y + offset_x + k - vector_len - offset_str_row];
				}
				//simd : shift
				//shifteamos los chars para que quede bien (el sentido es contrario a cuando lo hagamos en simd)
				for(int s = 0;s < offset_str_row; s++){
					for(int k = 1;k < vector_len;k++){
						str_row[k-1] = str_row[k];
					}
				}

			}else{ //caso feliz
			//	cerr<<"centro"<<endl;
				for(int k = 0;k < vector_len;k++){
					str_row[k] = seq1[offset_y + offset_x + k - vector_len];
				}
			}
			printf("%s\n", str_row);
			
			for(int k = 0;k < vector_len;k++){
				int val = 0;
				if(str_row[k]=='A')
					val=1;
				else if(str_row[k]=='G')
					val=2;
				else if(str_row[k]=='T')
					val=3;
				else if(str_row[k]=='U')
					val=4;
				else
					val=5;
				score_matrix[offset_y + offset_x + k] = val;
				count++;
			}			
		}
	}

	printScoreMatrix(score_matrix, &alignment, (int)vector_len);

	free(str_row);
	free(str_col);
}
