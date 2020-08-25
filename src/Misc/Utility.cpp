#include <iostream>
#include <string>
#include <string.h>
#include <vector>

#include "Utility.hpp"


#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;


void backtracking_C(
	short *score_matrix,
	Alignment& alignment,
	int vector_len,
	unsigned int best_x,unsigned int best_y,
	bool SW,
	score_fun_t score_fun,
	bool debug = false
){
	
	char* seq1 = alignment.sequence_1->sequence;
	char* seq2 = alignment.sequence_2->sequence;

	int seq1_len = alignment.sequence_1->length;
	int seq2_len = alignment.sequence_2->length;
	
	Sequence* best_seq_1 = alignment.result->sequence_1;
	Sequence* best_seq_2 = alignment.result->sequence_2;

	char* best_sequence_1 = (char*)malloc(seq1_len+1+seq2_len+1);
	char* best_sequence_2 = (char*)malloc(seq1_len+1+seq2_len+1);
	
	unsigned int length = 0;
	
	unsigned int x = best_x;
	unsigned int y = best_y;
	
	while(y != 0 || x != 0){
		// cerr<<"step"<<endl;
		// DBG(x);
		// DBG(y);
		if(y > 0 && x > 0){
			if(SW && score_fun(score_matrix, seq1_len,y,x,vector_len) == 0)
				break;

			short score_diag = score_fun(score_matrix, seq1_len,y-1,x-1,vector_len) + 
						  alignment.parameters->missmatch*(seq2[y] != seq1[x]) + 
						  alignment.parameters->match*(seq2[y] == seq1[x]);
			
			//short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = score_fun(score_matrix, seq1_len,y-1,x,vector_len) + alignment.parameters->gap;
			
			if(score_diag ==  score_fun(score_matrix, seq1_len,y,x,vector_len)){
				best_sequence_1[length] = seq2[y];
				best_sequence_2[length] = seq1[x];
				x--;
				y--;
			}else if(score_up == score_fun(score_matrix, seq1_len,y,x,vector_len)){
				best_sequence_1[length] = seq2[y];
				best_sequence_2[length] = '_';
				y--;
			}else{
				best_sequence_1[length] = '_';
				best_sequence_2[length] = seq1[x];
				x--;
			}
		}else if(x > 0){
			best_sequence_1[length] = '_';
			best_sequence_2[length] = seq2[x-1];
			x--;
		}else if(y > 0){
			best_sequence_1[length] = seq1[y-1];
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
	
	// cerr << "Best sequences" << endl;
	best_sequence_1[length]=0;
	best_sequence_2[length]=0;
	DBG(best_sequence_2);
	DBG(best_sequence_1);
	alignment.result->sequence_1 = new_Sequence_from_string(best_sequence_1);
	alignment.result->sequence_2 = new_Sequence_from_string(best_sequence_2);
	
	alignment.result->score = score_fun(score_matrix, seq1_len,best_y,best_x,vector_len);

	//creamos los nuevos strings y destruimos los anteriores
	best_seq_1 = new_Sequence_from_string(best_sequence_1);
	best_seq_2 = new_Sequence_from_string(best_sequence_2);

	destroy_Sequence(alignment.result->sequence_1);
	destroy_Sequence(alignment.result->sequence_2);
	alignment.result->sequence_1 = best_seq_1;
	alignment.result->sequence_2 = best_seq_2;

}


void printStartLine(int row, int vector_len){
	int espacios = (vector_len-1)-(row % vector_len);
	for	(int i = 0; i < espacios; i++){
		cout << "   ";
	}
}

void printEndLine(int row, int vector_len){
	int espacios = (row % vector_len);
	for	(int i = 0; i < espacios; i++){
		cout << "   ";
	}
	cout << endl;
}

void paintLine(vector<vector<char>>& canvasMatrix, int squareHeight, int squareWidth, int x, int y, int length, bool orientation){
	/*
		orientation = false ---> vertical
		orientation = true ----> horizontal
	*/
	
	if(orientation == true){
		for (int j = 0; j < (squareWidth-1)*length+1; j++)
		{
			if(j % (squareWidth-1) == 0)
				canvasMatrix[y][x+j] =  '@';
			else
				canvasMatrix[y][x+j] =  '@';
		}
	}else{
		for (int i = 0; i < (squareHeight-1)*length+1; i++)
		{
			if(i % (squareHeight-1) == 0)
				canvasMatrix[y+i][x] =  '@';
			else
				canvasMatrix[y+i][x] =  '@';
		}
	}
	
}

void paintSquare(vector<vector<char>>& canvasMatrix, int squareHeight, int squareWidth, int x, int y, string data){
	int size = data.size();
	int start_pos = (squareWidth - size)/2;
	int end_pos = start_pos + size;
	
	for (int i = 0; i < squareHeight; i++)
	{
		for (int j = 0; j < squareWidth; j++)
		{
			if (i == 0 || i == squareHeight-1)
			{
				// borde superior e inferior
				canvasMatrix[y+i][x+j] =  '-';
			} else if (j == 0 || j == squareWidth-1)
			{
				// borde lateral
				canvasMatrix[y+i][x+j] = '|';
			} else if (i == squareHeight / 2 && j >= start_pos && j < end_pos)
			{
				// tengo que imprimir el nro
				canvasMatrix[y+i][x+j] = data[j-start_pos];
			} else {
				// espacio vacio
				canvasMatrix[y+i][x+j] = ' ';
			}
		}
		
	}

	canvasMatrix[y][x] = '+';
	canvasMatrix[y][x+squareWidth-1] = '+';
	canvasMatrix[y+squareHeight-1][x] = '+';
	canvasMatrix[y+squareHeight-1][x+squareWidth-1] = '+'; 	
}

void printScoreMatrix2(short* matrix,Alignment* alignment,int vec){
	unsigned int seq1_len = alignment->sequence_1->length;
	unsigned int seq2_len = alignment->sequence_2->length;
	//en este caso hardcodeamos el tamaño del vector
	int vector_len = vec;
	
	int height = ((seq2_len + vector_len - 1)/ vector_len); //cantidad de "franjas" de diagonales
	int width = (1 + seq1_len + vector_len - 1); //cantidad de diagonales por franja
	
	vector<vector<int>> score_matrix(height*vector_len,vector<int>(width,-100000));
	
	for(int i = 0 ; i < height ; i++){
		unsigned int offset_y = i * width * vector_len;
		for( int j = 0; j < width ; j++){
			unsigned int offset_x = j * vector_len;

			for(int k = 0;k < vector_len;k++){
				score_matrix[i*vector_len+k][j] = matrix[offset_y+offset_x+vector_len-1-k];
			}			
		}
	}

	int squareWidth = 8;
	int squareHeight = 3;
	int canvasHeight = height*vector_len*squareHeight;
	int canvasWidth = (width+vector_len-1)*squareWidth;

	vector<vector<char>> canvasMatrix(canvasHeight,vector<char>(canvasWidth,' '));
	for(int i = 0 ; i < height ; i++){
		for(int k=0;k < vector_len; k++){
			for(int j = 0; j < width ; j++){
				int mx = j+(vector_len-1-(i*vector_len+k)%vector_len)+1;
				int my = (i*vector_len+k)+1;
				
				paintSquare(
					canvasMatrix,
					squareHeight,squareWidth,
					mx*(squareWidth-1),
					my*(squareHeight-1),
					to_string(score_matrix[i*vector_len+k][j])
				);	
			}	 
		}
	}

	for(unsigned int i = 0; i < seq1_len; i++){
		paintSquare(
			canvasMatrix,
			squareHeight,squareWidth,
			(i+vector_len+1)*(squareWidth-1),
			0,
			string(1,alignment->sequence_1->sequence[i])
		);	
	}

	for(unsigned int i = 0; i < seq2_len; i++){
		paintSquare(
			canvasMatrix,
			squareHeight,squareWidth,
			0,
			(i+1)*(squareHeight-1),
			string(1,alignment->sequence_2->sequence[i])
		);	
	}

	paintLine(
		canvasMatrix,
		squareHeight,squareWidth,
		(vector_len+1)*(squareWidth-1),(0+1)*(squareHeight-1),
		seq2_len,
		false);
	
	paintLine(
		canvasMatrix,
		squareHeight,squareWidth,
		(vector_len+1)*(squareWidth-1),(0+1)*(squareHeight-1),
		seq1_len,
		true);

	paintLine(
		canvasMatrix,
		squareHeight,squareWidth,
		(vector_len+seq1_len+1)*(squareWidth-1),(0+1)*(squareHeight-1),
		seq2_len,
		false);
		
	paintLine(
		canvasMatrix,
		squareHeight,squareWidth,
		(vector_len+1)*(squareWidth-1),(seq2_len+1)*(squareHeight-1),
		seq1_len,
		true);

	for(auto& row : canvasMatrix){
		for(char c : row){
			cout<<c;
		}cout<<endl;
	}
}

short get_score_SSE(short* score_matrix, int seq_row_len, int y,int x, int vector_len){
	int width = (1 + seq_row_len + vector_len - 1); //cantidad de diagonales por franja

	int franja = y / vector_len;	
	int offset_y = franja * width * vector_len;

	int diagonal = y % vector_len + x + 1;
	int cell = vector_len - 1 - y % vector_len;
	int offset_x = diagonal * vector_len + cell;

	return score_matrix[offset_y + offset_x];
}

short get_score_LIN(short* score_matrix, int seq_row_len, int y,int x, int vector_len){
	short** matrix_ptr = (short**) score_matrix;
	return matrix_ptr[y][x];
}

void print_arr(short* p, int len){
	for(int i=0;i<len;i++)
		cerr<<p[i]<<" ";
	cerr<<endl;
}

void print128_hex(__m128i var)
{
    uint16_t val[8];
    memcpy(val, &var, sizeof(val));
    printf("Reg content: %04X %04X %04X %04X %04X %04X %04X %04X \n", 
           val[0], val[1], val[2], val[3], val[4], val[5], 
           val[6], val[7]);
}

void reg_to_arr(short* p,__m128i reg){
	_mm_storeu_si128((__m128i*)p,reg);
}

void word_to_arr8(short* p,__m128i reg){
	_mm_storeu_si128((__m128i*)p,reg);
}


void word_to_char8(char* p,__m128i reg){
  __m128i a = _mm_packus_epi16(reg,reg);
  _mm_storeu_si64((__m128i*)p,a);
}

__m128i char_to_word8(char* p){
  __m128i a = _mm_loadl_epi64((__m128i*)p);
  __m128i b = _mm_setzero_si128();
	
  return _mm_unpacklo_epi8(a,b);
}

  