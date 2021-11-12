#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int edge (int a, int b) {
  assert (a != b);
  assert (a >  0);
  assert (b >  0);
  int min, max;
  min = a; max = b;
  if (min > max) { min = b; max = a; }

  return (max - 2) * (max - 1) / 2 + min; }

int main (int argc, char** argv) {


  FILE* file = fopen (argv[1], "r");
  int   size = atoi  (argv[2]);

  int cnf = 0;
  int nVar = 0, nCls = 0;
  int tmp = fscanf (file, " p cnf %i %i ", &nVar, &nCls);

  if (tmp == 2) {
    cnf = 1;
    printf ("p cnf %i %i\n", nVar, nCls);
  }

  int length = 0;
  int lit;
  while (1) {
    tmp = fscanf (file, " %i ", &lit);
    if (tmp == EOF) break;
    if ((length == 0) && (cnf == 0)) {
      length++;
      printf ("%i\t", lit); }
    else {
      if (lit != 0) {
        int var = abs(lit);
        int a = 1;
        int b = 1;
        int plus = var;
        for (int i = size - 1; i > 0; i--) {
          if (plus <= i) { b = a + plus; i = 0; }
          else { a++; plus -= i; }
        }
        if (lit < 0) printf("-");
//        printf ("%i %i %i\n", var, a, b);
        printf ("%i ", edge (a, b));
        length++;
      }
      if (lit == 0) {
        printf ("0\n");
        length = 0;
      }
    }
  }

}
