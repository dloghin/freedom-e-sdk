#include <stdint.h>
#include <stdlib.h>

#define DIM 100600000

void init(uint64_t* a, uint64_t size, uint64_t val) {
	uint64_t i;
	for (i = 0; i < size; i++)
		a[i] = val;
}

uint64_t sum(uint64_t* a, uint64_t size) {
	uint64_t sum = 0;
	uint64_t i;
	for (i = 0; i < size; i++)
		sum += a[i];
	return sum;
}

int main() { 
	uint64_t* a = (uint64_t*)malloc(DIM * sizeof(uint64_t));
	init(a, DIM, 1);
	uint64_t s1 = sum(a, DIM);
	init(a, DIM, 2);
	uint64_t s2 = sum(a, DIM);
	free(a);
	return (int)(s1+s2);
}
