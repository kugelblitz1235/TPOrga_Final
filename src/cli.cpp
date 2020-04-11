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
  //char* seq1 = "TGTTACGG";
  //char* seq2 = "GGTTGACTA";
  char* seq1 = "GCATGCU";
  char* seq2 = "GATTACA";
  short gap_pen=-1;
  short missmatch_pen=-1;
  short match_score=1;
  alignment->sequence_1=new_string(seq2);
  alignment->sequence_2=new_string(seq1);
  (alignment->parameters)->match=match_score;
  (alignment->parameters)->missmatch=missmatch_pen;
  (alignment->parameters)->gap=gap_pen;
  //SWLIN(alignment);
  NWLIN(alignment);
  String* res_seq_1 = alignment->result->sequence_1;
  cout<<res_seq_1->length<< endl;
  for(int i = 0;i < res_seq_1->length;i++){
      cout << res_seq_1->sequence[i];
  }cout << endl;
  
  String* res_seq_2 = alignment->result->sequence_2;
  for(int i = 0;i < res_seq_2->length;i++){
      cout << res_seq_2->sequence[i];
  }cout << endl;
  
  //cerr << (long long)matrix << endl;
  //for(int i = 0;i < 7+1;i++){
  //  for(int x = 0;x < 7+1;x++){
  //    cerr << matrix[i*(7+1)+x] << " ";
  //  }cerr << endl;
  //}cerr << endl;
  
 // char file[] = "./test.json";
  // save_object_as_JSON(a, file);
  
  destroy_alignment(a);

  return 0;
}
