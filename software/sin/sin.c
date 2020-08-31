/*
 * Approximation of sin(1) using the series:
 * sin(x) = x - x^3/3! + x^5/5! - x^7/7! + ...
 */

#include <stdint.h>
#include "../common/posit.h"
#include "../common/perf.h"

//#define PFDEBUG

#ifdef PFDEBUG
#include <stdio.h>
#endif

int main() {
	int index;
	int n = 10;
	init_constants();
	element_t sin = one;
	element_t fact = one;
	element_t sign = one;
	element_t i = two;

	unsigned long long startc = read_cycles();
	for(index = 1; index <= n; index++){
		sign = sign * minus_one;
		fact = fact * i * (i + one);
		i = i + two;
		sin = sin + sign * one / fact;
	}
	unsigned long long endc = read_cycles();
	endc = endc - startc;

#ifdef PFDEBUG
	printf("Cycles %lu\n", endc);
	printf("Sin(1) %x\n", sin);
#if defined(__x86_64__)
	printf("Result: %e\n", sin);
#endif	// __x86_64__
#endif
	return 0;
}
