// #include <iostream>
// #include <fstream>
// #include <string>
// #include <string.h>
// #include <vector>

#include "Utility.hpp"


#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;


extern "C" void backtracking_C(
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
	
	// Sequence* best_seq_1 = alignment.result->sequence_1;
	// Sequence* best_seq_2 = alignment.result->sequence_2;

	char* best_sequence_1 = (char*)malloc(seq1_len+1+seq2_len+1);
	char* best_sequence_2 = (char*)malloc(seq1_len+1+seq2_len+1);
	
	unsigned int length = 0;
	
	unsigned int x = best_x;
	unsigned int y = best_y;
	
	if(debug)printf("Comienza backtracking: x = %d, y = %d\n", x, y);
	while(y != 0 || x != 0){
		// cerr<<"step"<<endl;
		// DBG(x);
		// DBG(y);
		if(SW && score_fun(score_matrix, seq1_len,y,x,vector_len) == 0) {
			// printf("Termina por score = 0\n");
			break;
		}
		if(y > 0 && x > 0){

			short score_diag = score_fun(score_matrix, seq1_len,y-1,x-1,vector_len) + 
						  alignment.parameters->missmatch*(seq2[y] != seq1[x]) + 
						  alignment.parameters->match*(seq2[y] == seq1[x]);
			
			//short score_left = scores[y][x-1] + alignment.parameters->gap;
			
			short score_up = score_fun(score_matrix, seq1_len,y-1,x,vector_len) + alignment.parameters->gap;
			
			if(score_diag ==  score_fun(score_matrix, seq1_len,y,x,vector_len)){
				// printf("Diagonal\n");
				best_sequence_1[length] = seq1[x];
				best_sequence_2[length] = seq2[y];
				x--;
				y--;
			}else if(score_up == score_fun(score_matrix, seq1_len,y,x,vector_len)){
				// printf("Arriba\n");
				best_sequence_1[length] = '-';
				best_sequence_2[length] = seq2[y];
				y--;
			}else{
				// printf("Izquierda\n");
				best_sequence_1[length] = seq1[x];
				best_sequence_2[length] = '-';
				x--;
			}
		}else if(x > 0){
			// printf("Tope superior\nIzquierda\n");
			best_sequence_1[length] = seq1[x];
			best_sequence_2[length] = '-';
			x--;
		}else if(y > 0){
			// printf("Tope Izquierda\nArriba\n");
			best_sequence_1[length] = '-';
			best_sequence_2[length] = seq2[y];
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
	// DBG(best_sequence_1);
	// DBG(best_sequence_2);
	alignment.result->sequence_1 = new_Sequence_from_string(best_sequence_1);
	alignment.result->sequence_2 = new_Sequence_from_string(best_sequence_2);
	
	alignment.result->score = score_fun(score_matrix, seq1_len,best_y,best_x,vector_len);

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

void printScoreMatrix(short* matrix, Alignment* alignment, int vector_len, ofstream &ofs){
	unsigned int seq1_len = alignment->sequence_1->length;
	unsigned int seq2_len = alignment->sequence_2->length;
	//en este caso hardcodeamos el tamaño del vector
	// int vector_len = vec;
	
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
			ofs << c;
		}
		ofs << endl;
	}
}

extern "C" short get_score_SSE(short* score_matrix, unsigned int seq_row_len, int y,int x, int vector_len){
	int width = (1 + seq_row_len + vector_len - 1); //cantidad de diagonales por franja

	int franja = y / vector_len;	
	int offset_y = franja * width * vector_len;

	int diagonal = y % vector_len + x + 1;
	int cell = vector_len - 1 - y % vector_len;
	int offset_x = diagonal * vector_len + cell;

	return score_matrix[offset_y + offset_x];
}

short get_score_LIN(short* score_matrix, unsigned int seq_row_len, int y,int x, int vector_len){
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

bool valid_subsequence(Sequence* res_seq,Sequence* original_seq) {

	char* res_seq_clean = (char*) malloc(res_seq->length + 1);
	int res_seq_clean_len = 0;
	for(unsigned int i = 0; i < res_seq->length; i++){
		if (res_seq->sequence[i] != '-') {
			res_seq_clean[res_seq_clean_len] = res_seq->sequence[i];
			res_seq_clean_len++;
		}
	}
	res_seq_clean[res_seq_clean_len] = 0;

	return strstr(original_seq->sequence, res_seq_clean) != NULL;
}

bool valid_score(Alignment& alignment) {
	Parameters* params = alignment.parameters;
	int len = min(alignment.result->sequence_1->length, alignment.result->sequence_2->length);
	int score = 0;
	for (int i = 1; i < len; i++) {
		char seq1_char = alignment.result->sequence_1->sequence[i];
		char seq2_char = alignment.result->sequence_2->sequence[i];
		if(seq1_char == '-' || seq2_char == '-')
			score += params->gap;
		else{
			if(seq1_char == seq2_char)
				score += params->match;
			else
				score += params->missmatch;
		}
	}
	return score == alignment.result->score;
}

bool valid_alignment(Alignment& alignment){
	Result* result = alignment.result;
	
	bool valid_seq1 = valid_subsequence(result->sequence_1,alignment.sequence_1);
	bool valid_seq2 = valid_subsequence(result->sequence_2,alignment.sequence_2);
	bool valid_length = (result->sequence_1->length == result->sequence_2->length);
	bool valid_scr = valid_score(alignment);

	return valid_seq1 && valid_seq2 && valid_length && valid_scr;
}

bool check_scr_matrix_manual(short ** valid_matrix, Alignment* alignment, score_fun_t score_fun, int vector_len, bool verbose) {
    unsigned int seq1_len = alignment->sequence_1->length;
    unsigned int seq2_len = alignment->sequence_2->length;

	short (*matrix)[seq1_len] = (short (*)[seq1_len]) valid_matrix;

    bool matrix_ok = true;
    bool position_ok;
    for(unsigned int i = 0 ; i < seq2_len ; i++){
        for(unsigned int j = 0 ; j < seq1_len ; j++){
			// DBG(i);
			// DBG(j);
            position_ok = matrix[i][j] == score_fun(alignment->matrix,seq1_len,i,j,vector_len);
            // cout << "llegue" << endl;
			matrix_ok &= position_ok;
            if (!position_ok && verbose){
                printf("Score matrices differ at: \n");
                DBG(i);
                DBG(j);
                DBG(matrix[i][j]);
                DBG(score_fun(alignment->matrix,seq1_len,i,j,vector_len));
            }
        }
    }

    return matrix_ok;
}

extern "C" void print_registers(){
	long unsigned int rax,rbx,rcx,rdx,rdi,rsi,rbp,rsp,r8,r9,r10,r11,r12,r13,r14,r15;
	asm("mov %%rax, %0": "=r" (rax));
	asm("mov %%rbx, %0": "=r" (rbx));
	asm("mov %%rcx, %0": "=r" (rcx));
	asm("mov %%rdx, %0": "=r" (rdx));
	asm("mov %%rdi, %0": "=r" (rdi));
	asm("mov %%rsi, %0": "=r" (rsi));
	asm("mov %%rbp, %0": "=r" (rbp));
	asm("mov %%rsp, %0": "=r" (rsp));
	asm("mov %%r8, %0": "=r" (r8));
	asm("mov %%r9, %0": "=r" (r9));
	asm("mov %%r10, %0": "=r" (r10));
	asm("mov %%r11, %0": "=r" (r11));
	asm("mov %%r12, %0": "=r" (r12));
	asm("mov %%r13, %0": "=r" (r13));
	asm("mov %%r14, %0": "=r" (r14));
	asm("mov %%r15, %0": "=r" (r15));
	printf("---------------------------- Register dump ----------------------------\n");
	printf("rax: 0x%016lx, %ld\t\t\tr8:  0x%016lx, %ld\n",rax,rax,r8,r8);
	printf("rbx: 0x%016lx, %ld\t\t\tr9:  0x%016lx, %ld\n",rbx,rbx,r9,r9);
	printf("rcx: 0x%016lx, %ld\t\t\tr10: 0x%016lx, %ld\n",rcx,rcx,r10,r10);
	printf("rdx: 0x%016lx, %ld\t\t\tr11: 0x%016lx, %ld\n",rdx,rdx,r11,r11);
	printf("rdi: 0x%016lx, %ld\t\t\tr12: 0x%016lx, %ld\n",rdi,rdi,r12,r12);
	printf("rsi: 0x%016lx, %ld\t\t\tr13: 0x%016lx, %ld\n",rsi,rsi,r13,r13);
	printf("rbp: 0x%016lx\t\t\t\tr14: 0x%016lx, %ld\n",rbp,r14,r14);
	printf("rsp: 0x%016lx\t\t\t\tr15: 0x%016lx, %ld\n",rsp,r15,r15);
	printf("-----------------------------------------------------------------------\n");
}

extern "C" void print_xmm(short *stack, unsigned int n){
	printf("---------------------------- Register dump ----------------------------\n");
	for (unsigned int i = 0; i < n; i++){
		printf("xmm%02d: ", i);
		for (unsigned int j = 0; j < 8; j++) {
			printf("0x%04hX ", stack[i*8+j]);
		}
		cout << endl;
	}
	printf("-----------------------------------------------------------------------\n");
}