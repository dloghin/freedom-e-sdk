/**
 * Compute e using Euler's formula
 */
#define PFDEBUG
#define FIXEDN

#ifdef FIXEDN
#define N 20
#else
#include<stdlib.h>
#endif // FIXEDN

#include<stdint.h>

#ifdef PFDEBUG
#include <stdio.h>
#endif

int main(int argc, char** argv) {
  int n, i;
  float e, fact;

#ifdef FIXEDN
  n = N;
#else
  if (argc < 2) {
#ifdef PFDEBUG
    printf("Usage: %s <iterations>\n", argv[0]);
#endif // PFDEBUG
    return -1;
  }
  n = atoi(argv[1]);
#endif // FIXEDN

  e = 2.0;
  fact = 1.0;
  for (i = 2; i < n; i++) {
    fact = fact / i;
    e = e + fact;
  }
#ifdef PFDEBUG
  printf("Euler's number with %i iterations: ", n);
  printf("%9.8f\n", e);
#endif
	return 0;
}
