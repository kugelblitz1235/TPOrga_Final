#include <iostream>
#include <string>
#include <string.h>

#include "Utility.hpp"

using namespace std;


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
	int s1,s2;
	char* s1s;
	char* s2s;
	if(alignment->sequence_1->length >alignment->sequence_2->length){
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
