#include <unistd.h>
#include <iostream>
#include <string.h>
#include "Types.hpp"
#include "NW_C.hpp"
#include "SW_C.hpp"
#include "JSON.hpp"
#include "AlignAlgo.hpp"

using namespace std;

int
main (int argc, char **argv)
{
  Alignment* a = new_alignment();

  // Argument handling
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "a:s:t:p:q:r:")) != -1)
    switch (c)
      {
      case 'a':
        a->parameters->algorithm = new_string(optarg);
        break;
      case 's':
        a->sequence_1 = new_string(optarg);
        break;
      case 't':
        a->sequence_2 = new_string(optarg);
        break;
      case 'p':
        a->parameters->match = atoi(optarg);
        break;
      case 'q':
        a->parameters->missmatch = atoi(optarg);
        break;
      case 'e':
        a->parameters->gap = atoi(optarg);
        break;
      case '?':
        switch(optopt)
        {
        case 'a':
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          break;
        case 's':
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          break;
        case 't':
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          break;
        case 'p':
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          break;
        case 'q':
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          break;
        case 'r':
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          break;
        default:
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          break;
        }
        return 1;
      default:
        abort ();
      }

  if (a->parameters->algorithm == NULL | a->sequence_1 == NULL | a->sequence_1 == NULL)
  {
      cerr << "Missing arguments. Please specify an algorithm and two sequences (-a, -s, -t)." << endl;
      exit(1);
  }

  Alignment* alignment = new_alignment();
  char* seq1 = "TGTTACGG";
  char* seq2 = "GGTTGACTA";
  short gap_pen=-2;
  short missmatch_pen=-3;
  short match_score=3;
  short* matrix = (short*)SWLIN(seq2,seq1,match_score,missmatch_pen,gap_pen,9,8,alignment);
  //cerr << (long long)matrix << endl;
   
  for(int i = 0;i < 9+1;i++){
    for(int x = 0;x < 8+1;x++){
      cerr << matrix[i*(8+1)+x] << " ";
    }cerr << endl;
  }cerr << endl;
  
 // char file[] = "./test.json";
  // save_object_as_JSON(a, file);
  
  destroy_alignment(a);

  return 0;
}
