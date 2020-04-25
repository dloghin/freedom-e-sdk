#include "common2.h"

#include "/home/dumi/git/riscv/freedom-e-sdk-dloghin/software/common/perf.h"

#if (defined WITH_POSIT_32)
// posit(32,3)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x40000000;
uint32_t posit_two = 0x44000000;
uint32_t posit_three = 0x46000000;
uint32_t posit_four = 0x48000000;
uint32_t posit_minus_one = 0xc0000000;
uint32_t posit_e = 0x456fc2a3;
#elif (defined WITH_POSIT_16)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00004000;
uint32_t posit_two = 0x00004800;
uint32_t posit_three = 0x00004c00;
uint32_t posit_four = 0x00005000;
uint32_t posit_minus_one = 0x0000c000;
uint32_t posit_e = 0x00004ae0;
#elif (defined WITH_POSIT_8)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00000040;
uint32_t posit_two = 0x00000048;
uint32_t posit_three = 0x0000004c;
uint32_t posit_four = 0x00000050;
uint32_t posit_minus_one = 0x000000c0;
uint32_t posit_e = 0x00000055;
#else
uint32_t fp32_zero = 0x00000000;
uint32_t fp32_one = 0x3f800000;
uint32_t fp32_two = 0x40000000;
uint32_t fp32_three = 0x40400000;
uint32_t fp32_four = 0x40800000;
uint32_t fp32_minus_one = 0xbf800000;
uint32_t fp32_e = 0x402df854;
#endif /* WITH_POSIT */

// constants
element_t minus_one = -1.0;
element_t zero = 0.0;
element_t one = 1.0;
element_t two = 2.0;
element_t three = 3.0;
element_t four = 4.0;
element_t e = 2.7182818;

// data
element_t data[51200];

#ifdef WDEBUG
void write_data(element_t* data, size_t count, char* file) {
	FILE* f = fopen(file, "wb");
	if (!f) {
		printf("Error opening file %s\n", file);
		return;
	}
	size_t n = fwrite(data, sizeof(element_t), count, f);
	if (n != count) {
		printf("Wrote %ld elements, expected %ld elements!\n", n, count);
	}
	fclose(f);
}
#endif /* PDEBUG */

void init() {
#if (defined WITH_POSIT_8 || defined WITH_POSIT_16 || defined WITH_POSIT_32)
  *((uint32_t*)&zero) = posit_zero;
  *((uint32_t*)&one) = posit_one;
  *((uint32_t*)&minus_one) = posit_minus_one;
  *((uint32_t*)&e) = posit_e;
#else
  *((uint32_t*)&zero) = fp32_zero;
  *((uint32_t*)&one) = fp32_one;
  *((uint32_t*)&minus_one) = fp32_minus_one;
  *((uint32_t*)&e) = fp32_e;
#endif /* WITH_POSIT */
}

void get_max(element_t* a, size_t n) {
	element_t max = a[0];
	int idx = 0;
	for (size_t i = 1; i < n; i++) {
		if (a[i] > max) {
			max = a[i];
			idx = i;
		}
	}
	printf("%i %5.4f\n", idx, max);
}

int main() {
	size_t n;
	element_t* output;
	init();
	unsigned long long startc = read_cycles();
	forward(&output, &n);
	unsigned long long endc = read_cycles();
	get_max(output, n);
	return 0;
}
