#include <iostream>
#include <string>
#include <string.h>

#include "Utility.hpp"


#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;


void backtracking_C(
	short *score_matrix,
	Alignment& alignment,
	int vector_len,
	unsigned int x,unsigned int y,
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

	char* best_sequence_1 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	char* best_sequence_2 = (char*)malloc(alignment.sequence_1->length+1+alignment.sequence_2->length+1);
	
	unsigned int length = 0;
	
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
	// DBG(best_sequence_2);
	// DBG(best_sequence_1);
	alignment.result->sequence_1 = new_Sequence_from_string(best_sequence_1);
	alignment.result->sequence_2 = new_Sequence_from_string(best_sequence_2);
	alignment.result->score = score_fun(score_matrix, seq1_len,alignment.sequence_2->length-1,alignment.sequence_1->length-1,vector_len);

	//creamos los nuevos strings y destruimos los anteriores
	best_seq_1 = new_Sequence_from_string(best_sequence_1);
	best_seq_2 = new_Sequence_from_string(best_sequence_2);

	destroy_Sequence(alignment.result->sequence_1);
	destroy_Sequence(alignment.result->sequence_2);
	alignment.result->sequence_1 = best_seq_1;
	alignment.result->sequence_2 = best_seq_2;

}


void printSpaceLine(int dataSize,int columns){
	for(int x = 0;x < columns*(dataSize+1);x++){
		if(x%(dataSize+1)==0)cout<<"|";
		else cout << " ";
	}
	cout << "|" << endl;
}

void printDivisionLine(int dataSize,int columns){
	for(int x = 0;x < columns*(dataSize+1);x++){
		if(x%(dataSize+1)==0)cout<<"+";
		else cout << "-";
	}
	cout << "+" << endl;
}

void printScoreMatrix(short* matrix,Alignment* alignment,int vec){
	int s1 = alignment->sequence_1->length;
	int s2 = alignment->sequence_2->length;;
	char* s1s = alignment->sequence_1->sequence;
	char* s2s = alignment->sequence_2->sequence;
	/*if(alignment->sequence_1->length >alignment->sequence_2->length){
		s1 = alignment->sequence_1->length;
		s1s = alignment->sequence_1->sequence;
	}else{
		s1 = alignment->sequence_2->length;
		s1s = alignment->sequence_2->sequence;
	}
	
	if(alignment->sequence_1->length <= alignment->sequence_2->length){
		s2 = alignment->sequence_1->length;
		s2s = alignment->sequence_1->sequence;
	}else{
		s2 = alignment->sequence_2->length;
		s2s = alignment->sequence_2->sequence;
	}
	*/
	
	int dataSize = 6;
	
	int rows =  ((s2 + vec - 1) / vec)*vec;
	int columns = s1 + vec + vec - 1;
	
	printDivisionLine(dataSize,columns+1);
	for(int i = 0;i < dataSize/4;i++)
		printSpaceLine(dataSize,columns+1);
		
	for(int x = 0;x < columns;x++){
		cout << "|";
		string number = "";
		if(vec <= x && x < columns-vec+1)
			number = string(1,s1s[x-vec]);
		
		for(int i = 0;i < dataSize/2;i++)
			cout << " ";
		cout << number;
		for(int i = 0;i < (dataSize+1)/2;i++)
			cout << " ";
	}
	cout << "|" << endl;
	for(int i = 0;i < dataSize/4;i++)
		printSpaceLine(dataSize,columns+1);

	for(int y = 0;y < rows;y++){
		printDivisionLine(dataSize,columns+1);
		for(int i = 0;i < dataSize/4;i++)
			printSpaceLine(dataSize,columns+1);
		
		cout << "|";
		string number = "";
		if(y < s2)
			number = string(1,s2s[y]);
		int space = dataSize-number.size();
		
		for(int i = 0;i < space/2;i++)
			cout << " ";
		cout << number;
		for(int i = 0;i < (space+1)/2;i++)
			cout << " ";
			
		for(int x = 0;x < columns;x++){
			cout << "|";
			string number;
			if(x < (vec-1-(y%vec)) || (columns-(y%vec) - 1)< x){
				number = " ";
			}else{
				int xx = x-vec;
				int yy = y;
				int index = 0;
				index += (yy/vec)*(s1+vec)*vec;
				index += vec*(xx+4);
				index += (1-vec)*(vec-1-yy%vec);
				number = to_string(matrix[index]);
			}
			
			int space = dataSize-number.size();
			
			for(int i = 0;i < space/2;i++)
				cout << " ";
			cout << number;
			for(int i = 0;i < (space+1)/2;i++)
				cout << " ";
		}
		cout << "|" << endl;
		for(int i = 0;i < dataSize/4;i++)
			printSpaceLine(dataSize,columns+1);
	}
	printDivisionLine(dataSize,columns+1);
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

void print_arr(short* p, unsigned int len){
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

  