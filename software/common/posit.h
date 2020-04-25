#pragma once

// #define WITH_POSIT_32
// #define WITH_POSIT_16
// #define WITH_POSIT_8

typedef float element_t;
typedef int int_t;

// constants
element_t zero, one, minus_one, two, seven, pi, eps, minus_eps, hundred;
#if (defined WITH_POSIT_32)       // posit 32,3
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x40000000;
uint32_t posit_two = 0x44000000;
uint32_t posit_pi = 0x46487eca;
uint32_t posit_eps = 0xd000001;
uint32_t posit_minus_eps = 0xf2ffffff;
uint32_t posit_hundred = 0x5a400000;
uint32_t posit_seven = 0x4b000000;
#elif (defined WITH_POSIT_16)   // posit 16,2
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00004000;
uint32_t posit_two = 0x00004800;
uint32_t posit_pi = 0x00004c90;
uint32_t posit_eps = 0x00000701;
uint32_t posit_minus_eps = 0x0000f8ff;
uint32_t posit_hundred = 0x00006a40;
uint32_t posit_seven = 0x5600;
#elif (defined WITH_POSIT_8)    // posit 8,1
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00000040;
uint32_t posit_two = 0x00000050;
uint32_t posit_pi = 0x00000059;
uint32_t posit_eps = 0x00000003;		// 0.001
uint32_t posit_minus_eps = 0x000000fd;	// -0.001
uint32_t posit_hundred = 0x0000007a;
uint32_t posit_seven = 0x66;
#else
uint32_t fp32_zero = 0x00000000;
uint32_t fp32_one = 0x3f800000;
uint32_t fp32_minus_one = 0xbf800000;
uint32_t fp32_two = 0x40000000;
uint32_t fp32_pi = 0x40490fd9;
uint32_t fp32_eps = 0x358637bd; 		// 0.000001
uint32_t fp32_minus_eps = 0xb58637bd; 	// -0.000001
uint32_t fp32_hundred = 0x42c80000;
uint32_t fp32_seven = 0x40e00000;		// 7.0
#endif /* WITH_POSIT */

void init_constants() {
#if (defined WITH_POSIT_8 || defined WITH_POSIT_16 || defined WITH_POSIT_32)
  *((uint32_t*)&zero) = posit_zero;
  *((uint32_t*)&one) = posit_one;
  *((uint32_t*)&two) = posit_two;
  *((uint32_t*)&pi) = posit_pi;
  *((uint32_t*)&eps) = posit_eps;
  *((uint32_t*)&minus_eps) = posit_minus_eps;
  *((uint32_t*)&hundred) = posit_hundred;
  *((uint32_t*)&seven) = posit_seven;
#else
  *((uint32_t*)&zero) = fp32_zero;
  *((uint32_t*)&one) = fp32_one;
  *((uint32_t*)&minus_one) = fp32_minus_one;
  *((uint32_t*)&two) = fp32_two;
  *((uint32_t*)&pi) = fp32_pi;
  *((uint32_t*)&eps) = fp32_eps;
  *((uint32_t*)&minus_eps) = fp32_minus_eps;
  *((uint32_t*)&hundred) = fp32_hundred;
  *((uint32_t*)&seven) = fp32_seven;
#endif /* WITH_POSIT */
}
