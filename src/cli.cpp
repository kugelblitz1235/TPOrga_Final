#include <unistd.h>
#include <iostream>
#include <string.h>
#include <string>
#include "Misc/Types.hpp"
#include "Needleman_Wunsch/NW_C.hpp"
#include "Smith_Waterman/SW_C.hpp"
#include "Misc/JSON.hpp"
#include "Misc/AlignAlgo.hpp"
#include "Misc/Utility.hpp"

using namespace std;

//funcion que selecciona el algoritmo, la implementacion y recibe los parametros para ejecutarla
Alignment* align_sequences (std::string algorithm, std::string implementation, char* sequence_1, char* sequence_2, short gap, short missmatch, short match){
	
  try
  {
    
      //los cambio a mayusculas para una mejor comparacion
      for (auto & c: algorithm) c = toupper(c);
      for (auto & c: implementation) c = toupper(c);

      if(algorithm.compare("NW") == 0){
        
        return alignment_by_NW(implementation, sequence_1, sequence_2, gap, missmatch, match);

      }
      else if (algorithm.compare("SW") == 0)
      {
        return alignment_by_SW(implementation, sequence_1, sequence_2, gap, missmatch, match);
      }
      else{
        throw "Los parámetros ingresados son inválidos.";
      }

  }catch(char * ex){
    std::cout<<ex;
  }
	
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

  //Alignment* alignment = alignment_by_SW((char*) "TGTTACGG",(char*) "GGTTGACTA", -2,-3,3 );
/*Alignment* alignment = align_sequences("NW", "LIN",(char*) "GCATGCU",(char*) "GATTACA", -1,-1,1 );
  
  Sequence* res_seq_1 = alignment->result->sequence_1;
  cout<<res_seq_1->length<< endl;
  for(unsigned int i = 0;i < res_seq_1->length;i++){
      cout << res_seq_1->sequence[i];
  }cout << endl;
  
  Sequence* res_seq_2 = alignment->result->sequence_2;
  for(unsigned int i = 0;i < res_seq_2->length;i++){
      cout << res_seq_2->sequence[i];
  }cout << endl;
*/
  Alignment* alignment = new_alignment();
  alignment->sequence_1 = new_Sequence_from_string((char*) "GCATGCU");
	alignment->sequence_2 = new_Sequence_from_string((char*) "GATTACCA");
  NW_C_SSE(*alignment);
  
  /*
   cerr << (long long)matrix << endl;
   for(int i = 0;i < 9+1;i++){
     for(int x = 0;x < 8+1;x++){
       cerr << matrix[i*(8+1)+x] << " ";
     }cerr << endl;
   }cerr << endl;
  
   char file[] = "./test.json";
   save_object_as_JSON(a, file);
  
   destroy_alignment(a);
*/
//cd src && make clean && make && cd .. && ./cli
/*
  short scores[100000];
	for(int i = 0;i < 100000;i++)
	scores[i] = i-1000/(i+1);

  printScoreMatrix(scores,alignment,4);
*/
  return 0;
}
