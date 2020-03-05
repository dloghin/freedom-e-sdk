/**
 * Compute e using Euler's formula
 */
#define PFDEBUG

#include<stdlib.h>
#include<stdint.h>

#ifdef PFDEBUG
#include <stdio.h>
#endif

int main(int argc, char** argv) {
  int n, i;
  float e, fact;

  if (argc < 2) {
#ifdef PFDEBUG
    printf("Usage: %s <iterations>\n", argv[0]);
#endif
    return -1;
  }

  n = atoi(argv[1]);
  e = 2.0;
  fact = 1.0;
  for (i = 2; i < n; i++) {
    fact = fact / i;
    e = e + fact;
  }
#ifdef PFDEBUG
  printf("%9.8f\n", e);
#endif
	return 0;
}
