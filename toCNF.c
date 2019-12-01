#include <stdio.h>
#include <stdlib.h>

#define ALLOC	100000

#define POS(c,e)     ((c) * 2 * nEdge + (e) + 1)
#define NEG(c,e)     ((c) * 2 * nEdge + (e) + 1 + nEdge)
#define CANON(g)     ((g) + canonStart + 1)
#define KILLS(c,g)   ((c) + (g) * nCls + killsStart + 1)
#define SPLIT(c,m,l) ((c) * (size0 + size1) + (l) * (size0) + (m) + splitStart + 1)
#define ORDER(c,e)   ((c - nEdge + 1) * nEdge + (e) + orderStart + 1)

//#define BREAK
//#define SORT
//#define FIXCANON

int canonStart, killsStart, splitStart, orderStart;
int nNode, nEdge, nGraph, nClass;
int half, size0, size1;
int *set, seed;

void setUnits (int current, int depth, int max) {
  int i;
  for (i = current; i < nEdge; i++) {
    if (seed >= 0) set[depth] = i;
    if (depth + 1 == max) seed--;
    setUnits (i + 1, depth + 1, max); } }

void printSplit (int cls, int size, int offset, int part) {
   int i, j, e;
   for (i = 0; i < size; i++) {
      for (e = offset, j = 1; j < size; e++, j = j << 1) {
        if (i & j) printf("%i ", POS(cls,e));
        else       printf("%i ", NEG(cls,e)); }
      printf("%i 0\n", SPLIT(cls,i,part));
      for (e = offset, j = 1; j < size; e++, j = j << 1) {
        if (i & j) printf("-%i -%i 0\n", SPLIT(cls,i,part), POS(cls,e));
        else       printf("-%i -%i 0\n", SPLIT(cls,i,part), NEG(cls,e)); } } }

int main (int argc, char** argv) {
  int i, j, k, g, c, e, in;
  int nCls;
  int *mask, *eqcl;
  nNode = nGraph = nClass = nEdge  = 0;

  eqcl = (int*) malloc (sizeof(int) * ALLOC);
  mask = (int*) malloc (sizeof(int) * ALLOC);

  // parse map
  FILE* input;
  input = fopen(argv[1], "r");

  int first = 1;
  while (1) {
    int nat, tmp, cClass, cMask;
    tmp = fscanf (input, " %i ", &nat);
    if (tmp == EOF) break;

    if (first == 1) {
      cMask  = 0;
      cClass = nat;
      if (nat > nClass) nClass = nat; }

    else if (nat > 0) {
      cMask |= 1 << (nat - 1);
      if (nat > nEdge) nEdge = nat; }

    first = 0;
    if (nat == 0) {
      eqcl[nGraph] = cClass;
      mask[nGraph] = cMask;
      nGraph++;
      first = 1; } }

  if (nEdge ==  3) nNode = 3;
  if (nEdge ==  6) nNode = 4;
  if (nEdge == 10) nNode = 5;
  if (nEdge == 15) nNode = 6;

  if ((nEdge ==  6) && (nClass !=   11)) { printf("ERROR: not all classes present\n"); exit (0); }
  if ((nEdge == 10) && (nClass !=   34)) { printf("ERROR: not all classes present\n"); exit (0); }
  if ((nEdge == 15) && (nClass !=  156)) { printf("ERROR: not all classes present\n"); exit (0); }
  if ((nEdge == 21) && (nClass != 1044)) { printf("ERROR: not all classes present\n"); exit (0); }

  half  = (nEdge + 1) / 2;
  size0 = (1 << ((nEdge + 1)/2));
  size1 = (1 << ((nEdge    )/2));

  nCls = atoi (argv[2]);

  canonStart = nCls * 2 * nEdge;
  killsStart = canonStart + nGraph;
  splitStart = killsStart + nGraph * nCls;
//#ifdef SORT
  orderStart = splitStart + nCls * (size0 + size1);
//#else
//  orderStart = splitStart;
//#endif

  int totalCls = nGraph + nClass + nGraph * nCls * 3 + nCls * (size0 * (half + 1) + size1 * (nEdge - half + 1));

  for (i = 1; i <= nClass; i++) {
    for (in = 0, g = 0; g < nGraph; g++)
      if (eqcl[g] == i) {
        if (in == 4) { totalCls+= 6; orderStart++; in = 2; }
        in++; }
    for (j = 1; j <= 3; j++)
      if (in > j) for (k = 0; k < j; k++) totalCls++; }

  int totalVar = orderStart;
#ifdef SORT
  int totalVar += (nCls - nEdge + 1) * nEdge;
#endif

  if (argc > 3) {
    totalCls += nCls - nEdge + 1; }

#ifdef BREAK
  totalCls += (nEdge - 1) * (nEdge + 1);
#endif
#ifdef FIXCANON
  totalCls += nEdge - 3;
#endif
  printf ("p cnf %i %i\n", totalVar, totalCls);

  if (argc > 3) {
    seed = abs(atoi(argv[3]));
    int max  = nCls - nEdge + 1;
    set = (int *) malloc (sizeof(int) * max);
    int mod = 1;
    for (i = 1; i <= max; i++) {
      mod *= nEdge - i + 1;
      mod /= i; }
    seed = seed % mod;
    printf("c seed = %i \n", seed);
    setUnits (0, 0, max);
    for (i = 0; i < max; i++)
      printf("%i 0\n", POS(nEdge-1+i,set[i])); }

  // if no clause KILLS a graph, it must be the CANON; nGraph.
  for (g = 0; g < nGraph; g++) {
    for (c = 0; c < nCls; c++)
      printf ("%i ", KILLS(c,g));
    printf ("%i 0\n", CANON(g)); }

  // at least one graph of each isomorphism class must be CANON; nClass.
  for (i = 1; i <= nClass; i++) {
    for (g = 0; g < nGraph; g++)
      if (eqcl[g] == i) printf ("%i ", CANON(g));
    printf ("0\n"); }

  // define SPLIT variables.
  for (c = 0; c < nCls; c++) {
    printSplit (c, size0,    0, 0);
    printSplit (c, size1, half, 1); }

  // link SPLIT with CANON and KILLS variables; nGraph * nCls * 3.
  for (g = 0; g < nGraph; g++) {
    int mask0 = mask[g] & ((1 << half) - 1);
    int mask1 = mask[g] >> half;
    for (c = 0; c < nCls; c++) {
      printf ("-%i -%i -%i 0\n", SPLIT(c,mask0,0), SPLIT(c,mask1,1), CANON(g));
      printf ("%i -%i 0\n", SPLIT(c,mask0,0), KILLS(c,g));
      printf ("%i -%i 0\n", SPLIT(c,mask1,1), KILLS(c,g)); } }

  // at most one CANON per isomorphism class
  int next = splitStart + nCls * (size0 + size1) + 1;
  int set[4];
  for (i = 1; i <= nClass; i++) {
    for (in = 0, g = 0; g < nGraph; g++) {
      if (eqcl[g] == i) {
        if (in == 4) {
          printf ("%i %i 0\n", -set[0], -set[1]); printf ("%i %i 0\n", -set[0], -set[2]);
          printf ("%i %i 0\n", -set[1], -set[2]); printf ("%i %i 0\n", -set[0], -next  );
          printf ("%i %i 0\n", -set[1], -next  ); printf ("%i %i 0\n", -set[2], -next  );
          set[0] = -next++; set[1] = set[3]; in = 2; }
        set[in++] = CANON(g); } }
    for (j = 1; j <= 3; j++)
      if (in > j) for (k = 0; k < j; k++)
        printf ("%i %i 0\n", -set[j], -set[k]); }

#ifdef BREAK
  for (i = 1; i < nEdge; i++) {
    printf ("%i 0\n", NEG(i-1,i));
    printf ("-%i 0\n", POS(i-1,i));
    for (j = 0; j < nEdge; j++)
      if (j != i) printf ("-%i 0\n", NEG(i-1,j)); }
#endif

#ifdef FIXCANON
  int allone = (1 << nEdge) - 1;
  for (g = 0; g < nGraph; g++)
    for (e = 0; e < nEdge; e++)
      if (mask[g] != allone && e != 0 && e != 1 && e != (2 * nNode - 3))
        if ((mask[g] | (1 << e)) == allone)
          printf("-%i 0\n", CANON(g));
#endif
/*
  for (g = 0; g < nGraph; g++) {
    int count = 0;
    for (e = 0; e < nEdge; e++)
      if (mask[g] & (1 << e)) count++;
    if ((count > 4) && !(mask[g] & 1))
      printf("-%i 0\n", CANON(g)); }
*/
#ifdef SORT
  c = nEdge - 1;
  printf("-%i %i 0\n", ORDER(c,0), POS(c,0));
  printf("%i -%i 0\n", ORDER(c,0), POS(c,0));

  for (c = nEdge; c < nCls; c++) {
    printf("-%i -%i %i 0\n", ORDER(c-1,0), ORDER(c,0), POS(c,0));
    printf("-%i %i -%i 0\n", ORDER(c-1,0), ORDER(c,0), POS(c,0)); }

  for (c = nEdge - 1; c < nCls; c++)
    for (e = 1; e < nEdge; e++)
      printf("%i -%i %i 0\n", ORDER(c,e-1), POS(c,e), ORDER(c,e));

  for (c = nEdge - 1; c < nCls; c++)
    for (e = 0; e < nEdge - 1; e++) {
      printf("-%i %i 0\n", ORDER(c,e), ORDER(c,e+1));
      printf("%i -%i %i 0\n", ORDER(c,e), POS(c,e+1), ORDER(c,e+1)); }

  for (c = nEdge - 1; c < nCls - 1; c++)
    for (e = 0; e < nEdge; e++)
      printf("%i -%i 0\n", ORDER(c,e), POS(c+1,e));
#endif
}
