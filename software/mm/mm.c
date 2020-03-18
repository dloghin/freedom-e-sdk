#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#if !defined(__x86_64__)
#include <metal/cpu.h>
#endif

#define PFDEBUG
// #define WITH_POSIT_32

#define N 64

typedef float element_t;

#ifndef WITH_MALLOC
static element_t sa[N * N], sb[N * N], sc[N * N];
#endif

// constants
element_t zero, one, two, pi;
#if (defined WITH_POSIT_32)       // posit 32,3
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x40000000;
uint32_t posit_two = 0x44000000;
uint32_t posit_pi = 0x46487eca;
#elif (defined WITH_POSIT_16)   // posit 16,2
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00004000;
uint32_t posit_two = 0x00004800;
#elif (defined WITH_POSIT_8)    // posit 8,1
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00000040;
uint32_t posit_two = 0x00000050;
#else
uint32_t fp32_zero = 0x00000000;
uint32_t fp32_one = 0x3f800000;
uint32_t fp32_two = 0x40000000;
uint32_t fp32_pi = 0x40490fd9;
#endif /* WITH_POSIT */

void init() {
#if (defined WITH_POSIT_8 || defined WITH_POSIT_16 || defined WITH_POSIT_32)
  *((uint32_t*)&zero) = posit_zero;
  *((uint32_t*)&one) = posit_one;
  *((uint32_t*)&two) = posit_two;
  *((uint32_t*)&pi) = posit_pi;
#else
  *((uint32_t*)&zero) = fp32_zero;
  *((uint32_t*)&one) = fp32_one;
  *((uint32_t*)&two) = fp32_two;
  *((uint32_t*)&pi) = fp32_pi;
#endif /* WITH_POSIT */
}

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

unsigned long long read_cycles() {
#if defined(__x86_64__)
  return 0;
#else
  struct metal_cpu *mycpu = metal_cpu_get(0);
  return metal_cpu_get_timer(mycpu);
#endif
}

int main() {
  element_t *a, *b, *c;
  // init constants  
  init();
  // init matrix
  init_matrix(&a, &b, &c, N);
  // multiply
  unsigned long long startc = read_cycles();
  mm(a, b, c, N, N, N);
  unsigned long long endc = read_cycles();
  endc = endc - startc;
#ifdef PFDEBUG
  printf("Cycles %lu\n", endc);
  print_matrix(a, N, N);
  print_matrix(b, N, N);
  print_matrix(c, N, N);
#endif
#ifdef WITH_MALLOC
  free(a);
  free(b);
  free(c);
#endif
  return 0;
}
