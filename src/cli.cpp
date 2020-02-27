#include <unistd.h>
#include <iostream>
#include <cstring>
#include "Types.hpp"
#include "NW_C.hpp"
#include "SW_C.hpp"
#include "JSON.hpp"

using namespace std;

int
main (int argc, char **argv)
{
  // Default parameters
  Sequence s1;
  s1.length = 0;
  s1.sequence = NULL;
  Sequence s2;
  s2.length = 0;
  s2.sequence = NULL;
  Parameters p;
  p.match = 2;
  p.gap = -1;
  p.missmatch = -1;
  p.algorithm = NULL;
  Sequence rs1;
  rs1.length = 0;
  rs1.sequence = "TTTT";
  Sequence rs2;
  rs2.length = 0;
  rs2.sequence = "ATCG";
  Result r;
  r.sequence_1 = &rs1;
  r.sequence_2 = &rs2;
  r.score = 0;  
  Alignment a;
  a.sequence_1 = &s1;
  a.sequence_2 = &s2;
  a.parameters = &p;
  a.result = &r;  

  // Argument handling
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "a:s:t:p:q:r:")) != -1)
    switch (c)
      {
      case 'a':
        a.parameters->algorithm = optarg;
        break;
      case 's':
        a.sequence_1->sequence = optarg;
        a.sequence_1->length = strlen(optarg);
        break;
      case 't':
        a.sequence_2->sequence = optarg;
        a.sequence_2->length = strlen(optarg);
        break;
      case 'p':
        a.parameters->match = atoi(optarg);
        break;
      case 'q':
        a.parameters->missmatch = atoi(optarg);
        break;
      case 'e':
        a.parameters->gap = atoi(optarg);
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

  if (a.parameters->algorithm == NULL | a.sequence_1->sequence == NULL | a.sequence_1->sequence == NULL)
  {
      cerr << "Missing arguments." << endl;
      exit(1);
  }

  char file[] = "./test.json";
  save_object_as_JSON(&a, file);

  return 0;
}
