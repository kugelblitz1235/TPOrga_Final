#include <iostream>
#include <bitset>

#include "NW_C.hpp"

#define DBG(x) cerr << #x << " = " << (x) <<"\n"

using namespace std;

namespace NW{

	Alignment* get_alignment(std::string implementation, char* sequence_1, char* sequence_2, short gap, short missmatch, short match){
		// Creo la estructura vacia
		Alignment* alignment = new_alignment();
		
		// Relleno la estructura con los valores correspondientes
		alignment->sequence_1 = new_Sequence_from_string(sequence_1);
		alignment->sequence_2 = new_Sequence_from_string(sequence_2);
		(alignment->parameters)->match = match;
		(alignment->parameters)->missmatch = missmatch;
		(alignment->parameters)->gap = gap;

		if(implementation.compare("C_LIN") == 0)NW::C::LIN::NW(*alignment, true);
		else if(implementation.compare("C_SSE") == 0)NW::C::SSE::NW(*alignment, true);
		else if(implementation.compare("C_AVX512") == 0)NW::C::AVX512::NW(*alignment, true);
		else if(implementation.compare("C_AVX") == 0)NW::C::AVX::NW(*alignment, true);
		else if(implementation.compare("ASM_LIN") == 0)NW_ASM_LIN(alignment);
		else if(implementation.compare("ASM_SSE") == 0)NW_ASM_SSE(alignment, true);
		else if(implementation.compare("ASM_AVX") == 0)NW_ASM_AVX(alignment, true);
		else if(implementation.compare("ASM_AVX512") == 0)NW_ASM_AVX512(alignment, true);
		else throw "No existe la implementación ingresada.";
		

		//devuelvo la estructura modificada
		return alignment;
	}
}
