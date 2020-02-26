#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int
main (int argc, char **argv)
{
  char *algorithm = NULL;
  char *source = NULL;
  char *target = NULL;
  short match = 2;
  short missmatch = -1;
  short gap = -1;
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "a:s:t:p:q:r:")) != -1)
    switch (c)
      {
      case 'a':
        algorithm = optarg;
        break;
      case 's':
        source = optarg;
        break;
      case 't':
        target = optarg;
        break;
      case 'p':
        match = atoi(optarg);
        break;
      case 'q':
        missmatch = atoi(optarg);
        break;
      case 'e':
        gap = atoi(optarg);
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

  if (algorithm == NULL | source == NULL | target == NULL)
  {
      fprintf (stderr, "Missing arguments.\n");
      exit(1);
  }

  printf ("algorithm = %s, source = %s, target = %s, match = %i, missmatch = %i, gap = %i\n",
          algorithm, source, target, match, missmatch, gap);

  return 0;
}
