#include <unistd.h>
#include <iostream>
#include <string>
#include <map>
#include <emmintrin.h>

#include "Misc/Types.hpp"
#include "Smith_Waterman/SW_C_LIN.hpp"
#include "Smith_Waterman/SW_C_SIMDlogic.hpp"
#include "Smith_Waterman/SW_C_SSE.hpp"
#include "Smith_Waterman/SW_C_AVX.hpp"
#include "Smith_Waterman/SW_C_AVX512.hpp"
#include "Needleman_Wunsch/NW_C_LIN.hpp"
#include "Needleman_Wunsch/NW_C_SIMDlogic.hpp"
#include "Needleman_Wunsch/NW_C_SSE.hpp"
#include "Needleman_Wunsch/NW_C_AVX.hpp"
#include "Needleman_Wunsch/NW_C_AVX512.hpp"
#include "Misc/JSON.hpp"
#include "Misc/AlignAlgo.hpp"
#include "Misc/Utility.hpp"
#include "Misc/FASTA.hpp"

using namespace std;

map<string, string> algorithms = {
  {"SW_C_LIN", "Smith-Waterman implementado en C de forma lineal."},
  {"SW_C_SSE", "Smith-Waterman implementado en C usando instrucciones SSE."},
  {"SW_C_AVX", "Smith-Waterman implementado en C usando instrucciones AVX."},
  {"SW_C_AVX512", "Smith-Waterman implementado en C usando instrucciones AVX512."},
  {"SW_ASM_LIN", "Smith-Waterman implementado en ASM de forma lineal."},
  {"SW_ASM_SSE", "Smith-Waterman implementado en ASM usando instrucciones SSE."},
  {"SW_ASM_AVX", "Smith-Waterman implementado en ASM usando instrucciones AVX."},
  {"SW_ASM_AVX512", "Smith-Waterman implementado en ASM usando instrucciones AVX512."},
  {"NW_C_LIN", "Needlman-Wunsch implementado en C de forma lineal."},
  {"NW_C_SSE", "Needlman-Wunsch implementado en C usando instrucciones SSE."},
  {"NW_C_AVX", "Needlman-Wunsch implementado en C usando instrucciones AVX."},
  {"NW_C_AVX512", "Needlman-Wunsch implementado en C usando instrucciones AVX512."},
  {"NW_ASM_LIN", "Needlman-Wunsch implementado en ASM de forma lineal."},
  {"NW_ASM_SSE", "Needlman-Wunsch implementado en ASM usando instrucciones SSE."},
  {"NW_ASM_AVX", "Needlman-Wunsch implementado en ASM usando instrucciones AVX."},
  {"NW_ASM_AVX512", "Needlman-Wunsch implementado en ASM usando instrucciones AVX512."},
};

//funcion que selecciona el algoritmo, la implementacion y recibe los parametros para ejecutarla
void align_sequences(Alignment &alignment, string s1, string s2){
  FASTA_to_alignment(&alignment, s1.c_str(), s2.c_str());
  int min_length = min(alignment.sequence_1->length, alignment.sequence_2->length);
  string implementation = *(alignment.parameters->algorithm);

  if(implementation.find("SSE") != implementation.end() && min_len < 8){
    cerr < "Los algoritmos SSE necesitan strings de largo mínimo 8" << endl;
    exit(1);
  }
  else if(implementation.find("AVX512") != implementation.end() && min_len < 32){
    cerr < "Los algoritmos AVX512 necesitan strings de largo mínimo 32" << endl;
    exit(1);
  }
  else if(implementation.find("AVX") != implementation.end() && min_len < 16){
    cerr < "Los algoritmos AVX necesitan strings de largo mínimo 16" << endl;
    exit(1);
  }


  if(implementation.compare("SW_C_LIN") == 0)SW::C::LIN::SW(alignment, false);
  else if(implementation.compare("SW_C_SSE") == 0)SW::C::SSE::SW(alignment, false);
  else if(implementation.compare("SW_C_AVX") == 0)SW::C::AVX::SW(alignment, false);
  else if(implementation.compare("SW_C_AVX512") == 0)SW::C::AVX512::SW(alignment, false);
  else if(implementation.compare("SW_ASM_LIN") == 0)SW_ASM_LIN(&alignment);
  else if(implementation.compare("SW_ASM_SSE") == 0)SW_ASM_SSE(&alignment, false);
  else if(implementation.compare("SW_ASM_AVX") == 0)SW_ASM_AVX(&alignment, false);
  else if(implementation.compare("SW_ASM_AVX512") == 0)SW_ASM_AVX512(&alignment, false);
  else if(implementation.compare("NW_C_LIN") == 0)NW::C::LIN::NW(alignment, false);
  else if(implementation.compare("NW_C_SSE") == 0)NW::C::SSE::NW(alignment, false);
  else if(implementation.compare("NW_C_AVX512") == 0)NW::C::AVX512::NW(alignment, false);
  else if(implementation.compare("NW_C_AVX") == 0)NW::C::AVX::NW(alignment, false);
  else if(implementation.compare("NW_ASM_LIN") == 0)NW_ASM_LIN(&alignment);
  else if(implementation.compare("NW_ASM_SSE") == 0)NW_ASM_SSE(&alignment, false);
  else if(implementation.compare("NW_ASM_AVX") == 0)NW_ASM_AVX(&alignment, false);
  else if(implementation.compare("NW_ASM_AVX512") == 0)NW_ASM_AVX512(&alignment, false);
  else {
    throw "No existe la implementación ingresada.";
    exit(1);
  }
}

void imprimir_ayuda() {
  cerr << "Especificar algoritmo con -a" << endl;
  cerr << "Especificar puntaje de match con -p" << endl;
  cerr << "Especificar puntaje de mismatch con -q" << endl;
  cerr << "Especificar puntaje de gap con -r" << endl;
  cerr << "Los algoritmos existentes son: " << endl;
		for (auto& alg_desc: algorithms) cerr << "- " << alg_desc.first << "\t" << alg_desc.second << endl;

  exit(1);
}

int main (int argc, char **argv)
{
  Alignment a = *new_alignment();

  // Argument handling
  int c;

  opterr = 0;

  string algorithm;
  string s1, s2;
  
  while ((c = getopt(argc, argv, "a:s:t:p:q:r:h")) != -1)
    switch (c)
      {
      case 'a':
        algorithm = string(optarg);
        a.parameters->algorithm = &algorithm;
        break;
      case 's':
        s1 = string(optarg);
        break;
      case 't':
        s2 = string(optarg);
        break;
      case 'p':
        a.parameters->match = atoi(optarg);
        break;
      case 'q':
        a.parameters->missmatch = atoi(optarg);
        break;
      case 'r':
        a.parameters->gap = atoi(optarg);
        break;
      case 'h':
        imprimir_ayuda();
        break;
      case '?':
        switch(optopt)
        {
        case 'a':
          fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
          break;
        case 's':
          fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
          break;
        case 't':
          fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
          break;
        case 'p':
          fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
          break;
        case 'q':
          fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
          break;
        case 'r':
          fprintf (stderr, "La opción -%c requiere un argumento.\n", optopt);
          break;
        default:
          fprintf (stderr, "La opción `-%c' no existe. Imprimir ayuda con -h.\n", optopt);
          break;
        }
        return 1;
      default:
        imprimir_ayuda();
        exit(1);
      }

  if (a.parameters->algorithm == NULL | s1 == "" | s2 == "") {
    cerr << "Faltan argumentos obligatorios. Especificar un algoritmo y dos secuencias (-a, -s, -t).\nImprimir ayuda con -h.\n" << endl;
    exit(1);
  }

  if(algorithms.find(*(a.parameters->algorithm)) == algorithms.end()) {
    cout << "Algoritmo no encontrado: " << *(a.parameters->algorithm) << endl;
		cerr << "Los algoritmos existentes son: " << endl;
		for (auto& alg_desc: algorithms) cerr << "- " << alg_desc.first << "\t" << alg_desc.second << endl;
		exit(1);
  }

  align_sequences(a, s1, s2);

  printf("Alineamiento terminado:\nAlgoritmo: %s\n", (*a.parameters->algorithm).c_str());
  printf("Parámetros:\n\tmatch: %d\n\tmismatch: %d\n\tgap: %d\n", a.parameters->match, a.parameters->missmatch, a.parameters->gap);
  printf("Alineamiento:\n%s\n%s\n", a.result->sequence_1->sequence, a.result->sequence_2->sequence);
  printf("Puntaje:\n%d\n", a.result->score);

  return 0;
}