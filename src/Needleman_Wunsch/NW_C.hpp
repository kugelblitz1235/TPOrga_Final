#ifndef NW_C_H
#define NW_C_H
#include "../Misc/Types.hpp"
#include "../Misc/Utility.hpp"
#include "../Misc/AlignAlgo.hpp"
#include <string>
#include <limits.h>
#include <emmintrin.h>
#include <immintrin.h>

namespace NW{
    namespace C {
        //Version lineal del algoritmo utilizando la lógica de paralelización
        void with_logic_SSE (Alignment& alignment, int vector_len, bool debug);
        
        namespace LIN{
            //Version lineal del algoritmo de alineamiento global Needleman-Wunsch
            void NW (Alignment& alignment, bool debug);
        }
        
        namespace SSE{
            //Version paralelizada, con el set de instrucciones SSE, del algoritmo de alineamiento global Needleman-Wunsch
            void NW (Alignment& alignment, bool debug);
        }

        namespace AVX{
            //Version paralelizada, con el set de instrucciones AVX, del algoritmo de alineamiento global Needleman-Wunsch
            void NW (Alignment& alignment, bool debug);
        }

        namespace AVX512{
            //Version paralelizada, con el set de instrucciones AVX512, del algoritmo de alineamiento global Needleman-Wunsch
            void NW (Alignment& alignment, bool debug);
        }
    }
    //Método utilizado para seleccionar la versión del algoritmo a ejecutar
    Alignment* get_alignment(std::string implementation, char * sequence_1, char* sequence_2, short gap, short missmatch, short match);
}

#endif // NW_C_Ha