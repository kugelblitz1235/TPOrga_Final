#include <unistd.h>
#include <iostream>
#include <string>
#include <emmintrin.h>

#include "Misc/Types.hpp"
#include "Needleman_Wunsch/NW_C.hpp"
#include "Smith_Waterman/SW_C.hpp"
#include "Misc/JSON.hpp"
#include "Misc/AlignAlgo.hpp"
#include "Misc/Utility.hpp"
#include "Misc/FASTA.hpp"

using namespace std;

//funcion que selecciona el algoritmo, la implementacion y recibe los parametros para ejecutarla
Alignment* align_sequences (std::string algorithm, std::string implementation, char* sequence_1, char* sequence_2, short gap, short missmatch, short match){
    try
    {
          //los cambio a mayusculas para una mejor comparacion
          for (auto & c: algorithm) c = toupper(c);
          for (auto & c: implementation) c = toupper(c);

          if(algorithm.compare("NW") == 0){
            return NW::get_alignment(implementation, sequence_1, sequence_2, gap, missmatch, match);
          }
          else if (algorithm.compare("SW") == 0)
          {
            return SW::get_alignment(implementation, sequence_1, sequence_2, gap, missmatch, match);
          }
          else{
            throw "Los parámetros ingresados son inválidos.";
          }

    }
    catch(char * ex){
      std::cout<<ex;
    }
	return NULL;
}

int main (int argc, char **argv)
{
  // Alignment* a = new_alignment();

  // // Argument handling
  // int c;

  // opterr = 0;

  // while ((c = getopt(argc, argv, "a:s:t:p:q:r:")) != -1)
  //   switch (c)
  //     {
  //     case 'a':
  //       a->parameters->algorithm = new_string(optarg);
  //       break;
  //     case 's':
  //       a->sequence_1 = new_string(optarg);
  //       break;
  //     case 't':
  //       a->sequence_2 = new_string(optarg);
  //       break;
  //     case 'p':
  //       a->parameters->match = atoi(optarg);
  //       break;
  //     case 'q':
  //       a->parameters->missmatch = atoi(optarg);
  //       break;
  //     case 'e':
  //       a->parameters->gap = atoi(optarg);
  //       break;
  //     case '?':
  //       switch(optopt)
  //       {
  //       case 'a':
  //         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
  //         break;
  //       case 's':
  //         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
  //         break;
  //       case 't':
  //         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
  //         break;
  //       case 'p':
  //         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
  //         break;
  //       case 'q':
  //         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
  //         break;
  //       case 'r':
  //         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
  //         break;
  //       default:
  //         fprintf (stderr, "Unknown option `-%c'.\n", optopt);
  //         break;
  //       }
  //       return 1;
  //     default:
  //       abort ();
  //     }

  // if (a->parameters->algorithm == NULL | a->sequence_1 == NULL | a->sequence_1 == NULL)
  // {
  //     cerr << "Missing arguments. Please specify an algorithm and two sequences (-a, -s, -t)." << endl;
  //     exit(1);
  // }

  return 0;
}