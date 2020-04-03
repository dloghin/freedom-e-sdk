#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "../common/posit.h"
#include "../common/perf.h"

#define PFDEBUG

#define N 182

#ifndef WITH_MALLOC
static element_t sa[N * N], sb[N * N], sc[N * N];
#endif

void init_matrix(element_t** a, element_t** b, element_t** c, int n) {
  int i, j;

#ifdef WITH_MALLOC
  *a = (element_t*)malloc(n * n * sizeof(float));
  *b = (element_t*)malloc(n * n * sizeof(float));
  *c = (element_t*)malloc(n * n * sizeof(float));
#else
  *a = &(sa[0]);
  *b = &(sb[0]);
  *c = &(sc[0]);
#endif

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      (*a)[i*n+j] = zero;
      (*b)[i*n+j] = pi;
  }

  for (i = 0; i < n; i++)
    (*a)[i*n+i] = one;
}

void mm(element_t* a, element_t* b, element_t* c, int n, int m, int p) {
  int i, j, k;
  float r;
  for (i = 0; i < n; i++)
    for (j = 0; j < p; j++) {
      r = zero;
      for (k = 0; k < m; k++)
        r = r + a[i*m+k] * b[k*p+j];
      c[i*p+j] = r;
    }
}

#ifdef PFDEBUG
void print_matrix(element_t* a, int n, int m) {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++)
#if defined(__x86_64__)
      printf("%3.2f ", a[i*n+j]);
#else
      printf("%x ", a[i*n+j]);
#endif      
    printf("\n");
  }
  printf("\n");
}
#endif

int main() {
  element_t *a, *b, *c;
  // init constants  
  init_constants();
  // init matrix
  init_matrix(&a, &b, &c, N);
  // multiply
  unsigned long long startc = read_cycles();
  mm(a, b, c, N, N, N);
  unsigned long long endc = read_cycles();
  endc = endc - startc;
#ifdef PFDEBUG
  printf("Cycles %lu\n", endc);
//  print_matrix(a, N, N);
//  print_matrix(b, N, N);
//	print_matrix(c, N, N);
#endif
#ifdef WITH_MALLOC
  free(a);
  free(b);
  free(c);
#endif
  return 0;
}
