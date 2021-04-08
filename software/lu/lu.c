/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 3.0 structured OpenMP C versions - LU

  This benchmark is an OpenMP C version of the NPB LU code.
  
  The OpenMP C 2.3 versions are derived by RWCP from the serial Fortran versions 
  in "NPB 2.3-serial" developed by NAS. 3.0 translation is performed by the UVSQ.

  Permission to use, copy, distribute and modify this software for any
  purpose with or without fee is hereby granted.
  This software is provided "as is" without express or implied warranty.

  Information on OpenMP activities at RWCP is available at:

           http://pdplab.trc.rwcp.or.jp/pdperf/Omni/
  
  Information on NAS Parallel Benchmarks 2.3 is available at:
  
           http://www.nas.nasa.gov/NAS/NPB/

--------------------------------------------------------------------*/
/*--------------------------------------------------------------------

  Authors: S. Weeratunga
           V. Venkatakrishnan
           E. Barszcz
           M. Yarrow

  OpenMP C version: S. Satoh

  3.0 structure translation: M. Popov
  
--------------------------------------------------------------------*/

// https://github.com/benchmark-subsetting/NPB3.0-omp-C

/* global variables */
#include <stdint.h>

#include "../common/perf.h"

#include "applu.h"

//#define WITH_POSIT_32
//#define PFDEBUG

#ifdef PFDEBUG
#include <stdio.h>
#endif

// constants
element_t zero, one, two, three, four, five, six, ten, hundred, c_epsilon, c_e__8;

#if defined WITH_POSIT_32
/*
// posit(32,2)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x40000000;
uint32_t posit_two = 0x48000000;
 */

// posit(32,3)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x40000000;
uint32_t posit_two = 0x44000000;
uint32_t posit_three = 0x46000000;
uint32_t posit_four = 0x48000000;
uint32_t posit_five = 0x49000000;
uint32_t posit_six = 0x4a000000;
uint32_t posit_ten = 0x4d000000;
uint32_t posit_hundred = 0x5a400000;
uint32_t posit_thousand = 0x63e80000;
uint32_t posit_01 = 0x251eb851;
uint32_t posit_001 = 0x1c0c49ba;
uint32_t posit_0001 = 0x1546dc5d;
uint32_t posit_00001 = 0xf4f8b58;
uint32_t posit_e__8 = 0x6abcc77;
#elif (defined WITH_POSIT_16)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00004000;
uint32_t posit_two = 0x00004800;
uint32_t posit_three = 0x00004c00;
uint32_t posit_four = 0x00005000;
uint32_t posit_five = 0x00005200;
uint32_t posit_six = 0x00005400;
uint32_t posit_ten = 0x00005a00;
uint32_t posit_hundred = 0x00006a40;
uint32_t posit_thousand = 0x000073e8;
uint32_t posit_1 = 0x000024cc;
uint32_t posit_01 = 0x0000151e;
uint32_t posit_001 = 0x00000c0c;
uint32_t posit_0001 = 0x000006a3;
uint32_t posit_00001 = 0x000003a7;
uint32_t posit_e__8 = 0xaa;
#elif (defined WITH_POSIT_8)
uint32_t posit_zero = 0x00000000;
uint32_t posit_one = 0x00000040;
uint32_t posit_two = 0x00000050;
uint32_t posit_three = 0x00000058;
uint32_t posit_four = 0x00000060;
uint32_t posit_five = 0x00000062;
uint32_t posit_six = 0x00000064;
uint32_t posit_ten = 0x0000006a;
uint32_t posit_hundred = 0x00000079;
uint32_t posit_thousand = 0x0000007d;
uint32_t posit_1 = 0x00000014;
uint32_t posit_01 = 0x00000006;
uint32_t posit_001 = 0x00000002;
uint32_t posit_e__8 = 0x00000002;
#else
uint32_t fp32_zero = 0x00000000;
uint32_t fp32_one = 0x3f800000;
uint32_t fp32_two = 0x40000000;
uint32_t fp32_three = 0x40400000;
uint32_t fp32_four = 0x40800000;
uint32_t fp32_five = 0x40a00000;
uint32_t fp32_six = 0x40c00000;
uint32_t fp32_ten = 0x41200000;
uint32_t fp32_hundred = 0x42c80000;
uint32_t fp32_01 = 0x3c23d70a;
uint32_t fp32_001 = 0x3a83126f;
uint32_t fp32_0001 = 0x38d1b717;
uint32_t fp32_00001 = 0x3727c5ac;
uint32_t fp32_e__8 = 0x322bcc77;
#endif /* WITH_POSIT */

void init_constants() {
#if (defined WITH_POSIT_8 || defined WITH_POSIT_16 || defined WITH_POSIT_32)
	*((uint32_t*)&zero) = posit_zero;
	*((uint32_t*)&one) = posit_one;
	*((uint32_t*)&two) = posit_two;
	*((uint32_t*)&three) = posit_three;
	*((uint32_t*)&four) = posit_four;
	*((uint32_t*)&five) = posit_five;
	*((uint32_t*)&six) = posit_six;
	*((uint32_t*)&ten) = posit_ten;
	*((uint32_t*)&hundred) = posit_hundred;
	*((uint32_t*)&c_epsilon) = posit_1;
	*((uint32_t*)&c_e__8) = posit_e__8;
#else
	*((uint32_t*)&zero) = fp32_zero;
	*((uint32_t*)&one) = fp32_one;
	*((uint32_t*)&two) = fp32_two;
	*((uint32_t*)&three) = fp32_three;
	*((uint32_t*)&four) = fp32_four;
	*((uint32_t*)&five) = fp32_five;
	*((uint32_t*)&six) = fp32_six;
	*((uint32_t*)&ten) = fp32_ten;
	*((uint32_t*)&hundred) = fp32_hundred;
	*((uint32_t*)&c_epsilon) = fp32_01;
	*((uint32_t*)&c_e__8) = fp32_e__8;
#endif /* WITH_POSIT */
}

/*--------------------------------------------------------------------
c   parameters which can be overridden in runtime config file
c   isiz1,isiz2,isiz3 give the maximum size
c   ipr = 1 to print out verbose information
c   omega = 2.0 is correct for all classes
c   tolrsd is tolerance levels for steady state residuals
c-------------------------------------------------------------------*/

#define IPR_DEFAULT	0
#define	OMEGA_DEFAULT	(one + two/ten)

#define	TOLRSD1_DEF	c_e__8
#define	TOLRSD2_DEF	c_e__8
#define	TOLRSD3_DEF	c_e__8
#define	TOLRSD4_DEF	c_e__8
#define	TOLRSD5_DEF	c_e__8

#define	C1		(one + four/ten)
#define	C2		(four/ten)
#define	C3		(one/ten)
#define	C4		one
#define	C5		(one + four/ten)

#define	DT_DEFAULT	(one/two)

float sqrt_asm(float x) {
	float res;
#if defined(__x86_64__)
	asm("fsqrt" : "=t" (res) : "0" (x));
#else
	asm("fsqrt.s %0,%1\n\t" : "=f" (res) : "f" (x) : "cc");
#endif
	return res;
}

element_t my_fabs(element_t x) {
	if (x < zero)
		return -x;
	return x;
}

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define	pow2(a) ((a)*(a))

/* function declarations */
static void blts (int nx, int ny, int nz, int k,
		  element_t omega,
		  element_t v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		  element_t ldz[ISIZ1][ISIZ2][5][5],
		  element_t ldy[ISIZ1][ISIZ2][5][5],
		  element_t ldx[ISIZ1][ISIZ2][5][5],
		  element_t d[ISIZ1][ISIZ2][5][5],
		  int ist, int iend, int jst, int jend,
		  int nx0, int ny0 );
static void buts(int nx, int ny, int nz, int k,
		 element_t omega,
		 element_t v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		 element_t tv[ISIZ1][ISIZ2][5],
		 element_t d[ISIZ1][ISIZ2][5][5],
		 element_t udx[ISIZ1][ISIZ2][5][5],
		 element_t udy[ISIZ1][ISIZ2][5][5],
		 element_t udz[ISIZ1][ISIZ2][5][5],
		 int ist, int iend, int jst, int jend,
		 int nx0, int ny0 );
static int domain(void);
static void erhs(void);
static void error(void);
static void exact( int i, int j, int k, element_t u000ijk[5] );
static void jacld(int k);
static void jacu(int k);
static void l2norm (int nx0, int ny0, int nz0,
		    int ist, int iend,
		    int jst, int jend,
		    element_t v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		    element_t sum[5]);
static void pintgr(void);
static int read_input(void);
static void rhs(void);
static void setbv(void);
static void setcoeff(void);
static void setiv(void);
static void ssor(void);
static void verify(element_t xcr[5], element_t xce[5], element_t xci,
		   char *class, boolean *verified);

/*--------------------------------------------------------------------
      program applu
--------------------------------------------------------------------*/

int main(int argc, char **argv) {

/*--------------------------------------------------------------------
c
c   driver for the performance evaluation of the solver for
c   five coupled parabolic/elliptic partial differential equations.
c
--------------------------------------------------------------------*/

  char class;
  boolean verified;
  element_t mflops;
  int nthreads = 1;

  init_constants();

#ifdef PFDEBUG
	printf("\n\n NAS Parallel Benchmarks 3.0 structured OpenMP C version"
			" - SP Benchmark\n\n");
#endif

  unsigned long long startc = read_cycles();

/*--------------------------------------------------------------------
c   read input data
--------------------------------------------------------------------*/
  if (read_input())
	  return -1;

/*--------------------------------------------------------------------
c   set up domain sizes
--------------------------------------------------------------------*/
  if (domain())
	  return -1;

/*--------------------------------------------------------------------
c   set up coefficients
--------------------------------------------------------------------*/
  setcoeff();

/*--------------------------------------------------------------------
c   set the boundary values for dependent variables
--------------------------------------------------------------------*/
  setbv();

/*--------------------------------------------------------------------
c   set the initial values for dependent variables
--------------------------------------------------------------------*/
  setiv();

/*--------------------------------------------------------------------
c   compute the forcing term based on prescribed exact solution
--------------------------------------------------------------------*/
  erhs();
  
#pragma omp parallel
{  
  
#if defined(_OPENMP)  
#pragma omp master
  nthreads = omp_get_num_threads();
#endif /* _OPENMP */  
}

/*--------------------------------------------------------------------
c   perform the SSOR iterations
--------------------------------------------------------------------*/
  ssor();

/*--------------------------------------------------------------------
c   compute the solution error
--------------------------------------------------------------------*/
  error();

/*--------------------------------------------------------------------
c   compute the surface integral
--------------------------------------------------------------------*/
  pintgr();

/*--------------------------------------------------------------------
c   verification test
--------------------------------------------------------------------*/
  verify ( rsdnm, errnm, frc, &class, &verified );

  unsigned long long endc = read_cycles();
  endc = endc - startc;
#ifdef PFDEBUG
  printf("Cycles %lu\n", endc);
#endif

#ifdef PFDEBUG
	if (verified)
		printf("Verification: SUCCESS\n");
	else
		printf("Verification: FAILED\n");
#endif

	return verified;
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void blts (int nx, int ny, int nz, int k,
		  element_t omega,
/*--------------------------------------------------------------------
c   To improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.  Only needed in v.
--------------------------------------------------------------------*/
		  element_t v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		  element_t ldz[ISIZ1][ISIZ2][5][5],
		  element_t ldy[ISIZ1][ISIZ2][5][5],
		  element_t ldx[ISIZ1][ISIZ2][5][5],
		  element_t d[ISIZ1][ISIZ2][5][5],
		  int ist, int iend, int jst, int jend,
		  int nx0, int ny0 ) {
/*--------------------------------------------------------------------
c
c   compute the regular-sparse, block lower triangular solution:
c
c                     v <-- ( L-inv ) * v
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, m;
  element_t tmp, tmp1;
  element_t tmat[5][5];

#pragma omp for nowait schedule(static)
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (m = 0; m < 5; m++) {
	v[i][j][k][m] = v[i][j][k][m]
	  - omega * (  ldz[i][j][m][0] * v[i][j][k-1][0]
		       + ldz[i][j][m][1] * v[i][j][k-1][1]
		       + ldz[i][j][m][2] * v[i][j][k-1][2]
		       + ldz[i][j][m][3] * v[i][j][k-1][3]
		       + ldz[i][j][m][4] * v[i][j][k-1][4]  );
      }
    }
  }

#pragma omp for nowait schedule(static)
  for (i = ist; i <= iend; i++) {
    
#if defined(_OPENMP)      
    if (i != ist) {
	while (flag[i-1] == 0) {
#pragma omp flush(flag)
	    ;
	}
    }
    if (i != iend) {
	while (flag[i] == 1) {
#pragma omp flush(flag)
	    ;
	}
    }
#endif /* _OPENMP */
    
    for (j = jst; j <= jend; j++) {
      for (m = 0; m < 5; m++) {

	v[i][j][k][m] = v[i][j][k][m]
	  - omega * ( ldy[i][j][m][0] * v[i][j-1][k][0]
		      + ldx[i][j][m][0] * v[i-1][j][k][0]
		      + ldy[i][j][m][1] * v[i][j-1][k][1]
		      + ldx[i][j][m][1] * v[i-1][j][k][1]
		      + ldy[i][j][m][2] * v[i][j-1][k][2]
		      + ldx[i][j][m][2] * v[i-1][j][k][2]
		      + ldy[i][j][m][3] * v[i][j-1][k][3]
		      + ldx[i][j][m][3] * v[i-1][j][k][3]
		      + ldy[i][j][m][4] * v[i][j-1][k][4]
		      + ldx[i][j][m][4] * v[i-1][j][k][4] );
      }
       
/*--------------------------------------------------------------------
c   diagonal block inversion
c
c   forward elimination
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	tmat[m][0] = d[i][j][m][0];
	tmat[m][1] = d[i][j][m][1];
	tmat[m][2] = d[i][j][m][2];
	tmat[m][3] = d[i][j][m][3];
	tmat[m][4] = d[i][j][m][4];
      }

      tmp1 = one / tmat[0][0];
      tmp = tmp1 * tmat[1][0];
      tmat[1][1] =  tmat[1][1]
	- tmp * tmat[0][1];
      tmat[1][2] =  tmat[1][2]
	- tmp * tmat[0][2];
      tmat[1][3] =  tmat[1][3]
	- tmp * tmat[0][3];
      tmat[1][4] =  tmat[1][4]
	- tmp * tmat[0][4];
      v[i][j][k][1] = v[i][j][k][1]
	- v[i][j][k][0] * tmp;

      tmp = tmp1 * tmat[2][0];
      tmat[2][1] =  tmat[2][1]
	- tmp * tmat[0][1];
      tmat[2][2] =  tmat[2][2]
	- tmp * tmat[0][2];
      tmat[2][3] =  tmat[2][3]
	- tmp * tmat[0][3];
      tmat[2][4] =  tmat[2][4]
	- tmp * tmat[0][4];
      v[i][j][k][2] = v[i][j][k][2]
	- v[i][j][k][0] * tmp;

      tmp = tmp1 * tmat[3][0];
      tmat[3][1] =  tmat[3][1]
	- tmp * tmat[0][1];
      tmat[3][2] =  tmat[3][2]
	- tmp * tmat[0][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[0][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[0][4];
      v[i][j][k][3] = v[i][j][k][3]
	- v[i][j][k][0] * tmp;

      tmp = tmp1 * tmat[4][0];
      tmat[4][1] =  tmat[4][1]
	- tmp * tmat[0][1];
      tmat[4][2] =  tmat[4][2]
	- tmp * tmat[0][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[0][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[0][4];
      v[i][j][k][4] = v[i][j][k][4]
	- v[i][j][k][0] * tmp;

      tmp1 = one / tmat[ 1][1];
      tmp = tmp1 * tmat[ 2][1];
      tmat[2][2] =  tmat[2][2]
	- tmp * tmat[1][2];
      tmat[2][3] =  tmat[2][3]
	- tmp * tmat[1][3];
      tmat[2][4] =  tmat[2][4]
	- tmp * tmat[1][4];
      v[i][j][k][2] = v[i][j][k][2]
	- v[i][j][k][1] * tmp;

      tmp = tmp1 * tmat[3][1];
      tmat[3][2] =  tmat[3][2]
	- tmp * tmat[1][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[1][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[1][4];
      v[i][j][k][3] = v[i][j][k][3]
	- v[i][j][k][1] * tmp;

      tmp = tmp1 * tmat[4][1];
      tmat[4][2] =  tmat[4][2]
	- tmp * tmat[1][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[1][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[1][4];
      v[i][j][k][4] = v[i][j][k][4]
	- v[i][j][k][1] * tmp;

      tmp1 = one / tmat[2][2];
      tmp = tmp1 * tmat[3][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[2][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[2][4];
      v[i][j][k][3] = v[i][j][k][3]
        - v[i][j][k][2] * tmp;

      tmp = tmp1 * tmat[4][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[2][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[2][4];
      v[i][j][k][4] = v[i][j][k][4]
	- v[i][j][k][2] * tmp;

      tmp1 = one / tmat[3][3];
      tmp = tmp1 * tmat[4][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[3][4];
      v[i][j][k][4] = v[i][j][k][4]
	- v[i][j][k][3] * tmp;

/*--------------------------------------------------------------------
c   back substitution
--------------------------------------------------------------------*/
      v[i][j][k][4] = v[i][j][k][4]
	/ tmat[4][4];

      v[i][j][k][3] = v[i][j][k][3]
	- tmat[3][4] * v[i][j][k][4];
      v[i][j][k][3] = v[i][j][k][3]
	/ tmat[3][3];

      v[i][j][k][2] = v[i][j][k][2]
	- tmat[2][3] * v[i][j][k][3]
	- tmat[2][4] * v[i][j][k][4];
      v[i][j][k][2] = v[i][j][k][2]
	/ tmat[2][2];

      v[i][j][k][1] = v[i][j][k][1]
	- tmat[1][2] * v[i][j][k][2]
	- tmat[1][3] * v[i][j][k][3]
	- tmat[1][4] * v[i][j][k][4];
      v[i][j][k][1] = v[i][j][k][1]
	/ tmat[1][1];

      v[i][j][k][0] = v[i][j][k][0]
	- tmat[0][1] * v[i][j][k][1]
	- tmat[0][2] * v[i][j][k][2]
	- tmat[0][3] * v[i][j][k][3]
	- tmat[0][4] * v[i][j][k][4];
      v[i][j][k][0] = v[i][j][k][0]
	/ tmat[0][0];
    }
    
#if defined(_OPENMP)    
    if (i != ist) flag[i-1] = 0;
    if (i != iend) flag[i] = 1;
#pragma omp flush(flag)    
#endif /* _OPENMP */    
  }
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void buts(int nx, int ny, int nz, int k,
		 element_t omega,
/*--------------------------------------------------------------------
c   To improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.  Only needed in v.
--------------------------------------------------------------------*/
		 element_t v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		 element_t tv[ISIZ1][ISIZ2][5],
		 element_t d[ISIZ1][ISIZ2][5][5],
		 element_t udx[ISIZ1][ISIZ2][5][5],
		 element_t udy[ISIZ1][ISIZ2][5][5],
		 element_t udz[ISIZ1][ISIZ2][5][5],
		 int ist, int iend, int jst, int jend,
		 int nx0, int ny0 ) {
/*--------------------------------------------------------------------
c
c   compute the regular-sparse, block upper triangular solution:
c
c                     v <-- ( U-inv ) * v
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, m;
  element_t tmp, tmp1;
  element_t tmat[5][5];

#pragma omp for nowait schedule(static)
  for (i = iend; i >= ist; i--) {
    for (j = jend; j >= jst; j--) {
      for (m = 0; m < 5; m++) {
	tv[i][j][m] = 
	  omega * (  udz[i][j][m][0] * v[i][j][k+1][0]
		     + udz[i][j][m][1] * v[i][j][k+1][1]
		     + udz[i][j][m][2] * v[i][j][k+1][2]
		     + udz[i][j][m][3] * v[i][j][k+1][3]
		     + udz[i][j][m][4] * v[i][j][k+1][4] );
      }
    }
  }

#pragma omp for nowait schedule(static)
  for (i = iend; i >= ist; i--) {
#if defined(_OPENMP)      
    if (i != iend) {
      while (flag[i+1] == 0) {
#pragma omp flush(flag)
	;
      }
    }
    if (i != ist) {
      while (flag[i] == 1) {
#pragma omp flush(flag)
	;
      }
    }
#endif /* _OPENMP */
    
    for (j = jend; j >= jst; j--) {
      for (m = 0; m < 5; m++) {
	tv[i][j][m] = tv[i][j][m]
	  + omega * ( udy[i][j][m][0] * v[i][j+1][k][0]
		      + udx[i][j][m][0] * v[i+1][j][k][0]
		      + udy[i][j][m][1] * v[i][j+1][k][1]
		      + udx[i][j][m][1] * v[i+1][j][k][1]
		      + udy[i][j][m][2] * v[i][j+1][k][2]
		      + udx[i][j][m][2] * v[i+1][j][k][2]
		      + udy[i][j][m][3] * v[i][j+1][k][3]
		      + udx[i][j][m][3] * v[i+1][j][k][3]
		      + udy[i][j][m][4] * v[i][j+1][k][4]
		      + udx[i][j][m][4] * v[i+1][j][k][4] );
      }

/*--------------------------------------------------------------------
c   diagonal block inversion
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	tmat[m][0] = d[i][j][m][0];
	tmat[m][1] = d[i][j][m][1];
	tmat[m][2] = d[i][j][m][2];
	tmat[m][3] = d[i][j][m][3];
	tmat[m][4] = d[i][j][m][4];
      }

      tmp1 = one / tmat[0][0];
      tmp = tmp1 * tmat[1][0];
      tmat[1][1] =  tmat[1][1]
	- tmp * tmat[0][1];
      tmat[1][2] =  tmat[1][2]
	- tmp * tmat[0][2];
      tmat[1][3] =  tmat[1][3]
	- tmp * tmat[0][3];
      tmat[1][4] =  tmat[1][4]
	- tmp * tmat[0][4];
      tv[i][j][1] = tv[i][j][1]
	- tv[i][j][0] * tmp;

      tmp = tmp1 * tmat[2][0];
      tmat[2][1] =  tmat[2][1]
	- tmp * tmat[0][1];
      tmat[2][2] =  tmat[2][2]
	- tmp * tmat[0][2];
      tmat[2][3] =  tmat[2][3]
	- tmp * tmat[0][3];
      tmat[2][4] =  tmat[2][4]
	- tmp * tmat[0][4];
      tv[i][j][2] = tv[i][j][2]
	- tv[i][j][0] * tmp;

      tmp = tmp1 * tmat[3][0];
      tmat[3][1] =  tmat[3][1]
	- tmp * tmat[0][1];
      tmat[3][2] =  tmat[3][2]
	- tmp * tmat[0][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[0][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[0][4];
      tv[i][j][3] = tv[i][j][3]
	- tv[i][j][0] * tmp;

      tmp = tmp1 * tmat[4][0];
      tmat[4][1] =  tmat[4][1]
	- tmp * tmat[0][1];
      tmat[4][2] =  tmat[4][2]
	- tmp * tmat[0][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[0][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[0][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][0] * tmp;

      tmp1 = one / tmat[1][1];
      tmp = tmp1 * tmat[2][1];
      tmat[2][2] =  tmat[2][2]
	- tmp * tmat[1][2];
      tmat[2][3] =  tmat[2][3]
	- tmp * tmat[1][3];
      tmat[2][4] =  tmat[2][4]
	- tmp * tmat[1][4];
      tv[i][j][2] = tv[i][j][2]
	- tv[i][j][1] * tmp;

      tmp = tmp1 * tmat[3][1];
      tmat[3][2] =  tmat[3][2]
	- tmp * tmat[1][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[1][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[1][4];
      tv[i][j][3] = tv[i][j][3]
	- tv[i][j][1] * tmp;

      tmp = tmp1 * tmat[4][1];
      tmat[4][2] =  tmat[4][2]
	- tmp * tmat[1][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[1][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[1][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][1] * tmp;

      tmp1 = one / tmat[2][2];
      tmp = tmp1 * tmat[3][2];
      tmat[3][3] =  tmat[3][3]
	- tmp * tmat[2][3];
      tmat[3][4] =  tmat[3][4]
	- tmp * tmat[2][4];
      tv[i][j][3] = tv[i][j][3]
	- tv[i][j][2] * tmp;

      tmp = tmp1 * tmat[4][2];
      tmat[4][3] =  tmat[4][3]
	- tmp * tmat[2][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[2][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][2] * tmp;

      tmp1 = one / tmat[3][3];
      tmp = tmp1 * tmat[4][3];
      tmat[4][4] =  tmat[4][4]
	- tmp * tmat[3][4];
      tv[i][j][4] = tv[i][j][4]
	- tv[i][j][3] * tmp;

/*--------------------------------------------------------------------
c   back substitution
--------------------------------------------------------------------*/
      tv[i][j][4] = tv[i][j][4]
	/ tmat[4][4];

      tv[i][j][3] = tv[i][j][3]
	- tmat[3][4] * tv[i][j][4];
      tv[i][j][3] = tv[i][j][3]
	/ tmat[3][3];

      tv[i][j][2] = tv[i][j][2]
	- tmat[2][3] * tv[i][j][3]
	- tmat[2][4] * tv[i][j][4];
      tv[i][j][2] = tv[i][j][2]
	/ tmat[2][2];

      tv[i][j][1] = tv[i][j][1]
	- tmat[1][2] * tv[i][j][2]
	- tmat[1][3] * tv[i][j][3]
	- tmat[1][4] * tv[i][j][4];
      tv[i][j][1] = tv[i][j][1]
	/ tmat[1][1];

      tv[i][j][0] = tv[i][j][0]
	- tmat[0][1] * tv[i][j][1]
	- tmat[0][2] * tv[i][j][2]
	- tmat[0][3] * tv[i][j][3]
	- tmat[0][4] * tv[i][j][4];
      tv[i][j][0] = tv[i][j][0]
	/ tmat[0][0];

      v[i][j][k][0] = v[i][j][k][0] - tv[i][j][0];
      v[i][j][k][1] = v[i][j][k][1] - tv[i][j][1];
      v[i][j][k][2] = v[i][j][k][2] - tv[i][j][2];
      v[i][j][k][3] = v[i][j][k][3] - tv[i][j][3];
      v[i][j][k][4] = v[i][j][k][4] - tv[i][j][4];
    }
    
#if defined(_OPENMP)    
    if (i != iend) flag[i+1] = 0;
    if (i != ist) flag[i] = 1;
#pragma omp flush(flag)
#endif /* _OPENMP */    
  }
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static int domain(void) {

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/

  nx = nx0;
  ny = ny0;
  nz = nz0;

/*--------------------------------------------------------------------
c   check the sub-domain size
--------------------------------------------------------------------*/
  if ( nx < 4 || ny < 4 || nz < 4 ) {
#ifdef PFDEBUG
    printf("     SUBDOMAIN SIZE IS TOO SMALL - \n"
	   "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n"
	   "     SO THAT NX, NY AND NZ ARE GREATER THAN OR EQUAL\n"
	   "     TO 4 THEY ARE CURRENTLY%3d%3d%3d\n", nx, ny, nz);
#endif
    return -1;
  }

  if ( nx > ISIZ1 || ny > ISIZ2 || nz > ISIZ3 ) {
#ifdef PFDEBUG
    printf("     SUBDOMAIN SIZE IS TOO LARGE - \n"
	   "     ADJUST PROBLEM SIZE OR NUMBER OF PROCESSORS\n"
	   "     SO THAT NX, NY AND NZ ARE LESS THAN OR EQUAL TO \n"
	   "     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY.  THEY ARE\n"
	   "     CURRENTLY%4d%4d%4d\n", nx, ny, nz);
#endif
    return -1;
  }

/*--------------------------------------------------------------------
c   set up the start and end in i and j extents for all processors
--------------------------------------------------------------------*/
  ist = 1;
  iend = nx - 2;

  jst = 1;
  jend = ny - 2;

  return 0;
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void erhs(void) {

#pragma omp parallel
{

/*--------------------------------------------------------------------
c
c   compute the right hand side based on exact solution
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int iglob, jglob;
  int L1, L2;
  int ist1, iend1;
  int jst1, jend1;
  element_t  dsspm;
  element_t  xi, eta, zeta;
  element_t  q;
  element_t  u21, u31, u41;
  element_t  tmp;
  element_t  u21i, u31i, u41i, u51i;
  element_t  u21j, u31j, u41j, u51j;
  element_t  u21k, u31k, u41k, u51k;
  element_t  u21im1, u31im1, u41im1, u51im1;
  element_t  u21jm1, u31jm1, u41jm1, u51jm1;
  element_t  u21km1, u31km1, u41km1, u51km1;

  dsspm = dssp;

#pragma omp for
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      for (k = 0; k < nz; k++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = zero;
	}
      }
    }
  }

#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    xi = ( (element_t)(iglob) ) / ( nx0 - 1 );
    for (j = 0; j < ny; j++) {
      jglob = j;
      eta = ( (element_t)(jglob) ) / ( ny0 - 1 );
      for (k = 0; k < nz; k++) {
	zeta = ( (element_t)(k) ) / ( nz - 1 );
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  ce[m][0]
	    + ce[m][1] * xi
	    + ce[m][2] * eta
	    + ce[m][3] * zeta
	    + ce[m][4] * xi * xi
	    + ce[m][5] * eta * eta
	    + ce[m][6] * zeta * zeta
	    + ce[m][7] * xi * xi * xi
	    + ce[m][8] * eta * eta * eta
	    + ce[m][9] * zeta * zeta * zeta
	    + ce[m][10] * xi * xi * xi * xi
	    + ce[m][11] * eta * eta * eta * eta
	    + ce[m][12] * zeta * zeta * zeta * zeta;
	}
      }
    }
  }

/*--------------------------------------------------------------------
c   xi-direction flux differences
--------------------------------------------------------------------*/

  L1 = 0;
  L2 = nx-1;

#pragma omp for
  for (i = L1; i <= L2; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k < nz - 1; k++) {
	flux[i][j][k][0] = rsd[i][j][k][1];
	u21 = rsd[i][j][k][1] / rsd[i][j][k][0];
	q = one/two * (  rsd[i][j][k][1] * rsd[i][j][k][1]
		      + rsd[i][j][k][2] * rsd[i][j][k][2]
		      + rsd[i][j][k][3] * rsd[i][j][k][3] )
	  / rsd[i][j][k][0];
	flux[i][j][k][1] = rsd[i][j][k][1] * u21 + C2 * 
	  ( rsd[i][j][k][4] - q );
	flux[i][j][k][2] = rsd[i][j][k][2] * u21;
	flux[i][j][k][3] = rsd[i][j][k][3] * u21;
	flux[i][j][k][4] = ( C1 * rsd[i][j][k][4] - C2 * q ) * u21;
      }
    }
  }

#pragma omp for
  for (j = jst; j <= jend; j++) {
    for (k = 1; k <= nz - 2; k++) {
      for (i = ist; i <= iend; i++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] =  frct[i][j][k][m]
	    - tx2 * ( flux[i+1][j][k][m] - flux[i-1][j][k][m] );
	}
      }
      for (i = ist; i <= L2; i++) {
	tmp = one / rsd[i][j][k][0];

	u21i = tmp * rsd[i][j][k][1];
	u31i = tmp * rsd[i][j][k][2];
	u41i = tmp * rsd[i][j][k][3];
	u51i = tmp * rsd[i][j][k][4];

	tmp = one / rsd[i-1][j][k][0];

	u21im1 = tmp * rsd[i-1][j][k][1];
	u31im1 = tmp * rsd[i-1][j][k][2];
	u41im1 = tmp * rsd[i-1][j][k][3];
	u51im1 = tmp * rsd[i-1][j][k][4];

	flux[i][j][k][1] = (four/three) * tx3 *
	  ( u21i - u21im1 );
	flux[i][j][k][2] = tx3 * ( u31i - u31im1 );
	flux[i][j][k][3] = tx3 * ( u41i - u41im1 );
	flux[i][j][k][4] = (one/two) * ( one - C1*C5 )
	  * tx3 * ( ( u21i * u21i + u31i * u31i + u41i * u41i )
		    - ( u21im1*u21im1 + u31im1*u31im1 + u41im1*u41im1 ) )
	  + (one/six)
	  * tx3 * ( u21i*u21i - u21im1*u21im1 )
	  + C1 * C5 * tx3 * ( u51i - u51im1 );
      }

      for (i = ist; i <= iend; i++) {
	frct[i][j][k][0] = frct[i][j][k][0]
	  + dx1 * tx1 * (            rsd[i-1][j][k][0]
				     - two * rsd[i][j][k][0]
				     +       	    rsd[i+1][j][k][0] );
	frct[i][j][k][1] = frct[i][j][k][1]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][1] - flux[i][j][k][1] )
	  + dx2 * tx1 * (            rsd[i-1][j][k][1]
				     - two * rsd[i][j][k][1]
				     +           rsd[i+1][j][k][1] );
	frct[i][j][k][2] = frct[i][j][k][2]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][2] - flux[i][j][k][2] )
	  + dx3 * tx1 * (            rsd[i-1][j][k][2]
				     - two * rsd[i][j][k][2]
				     +           rsd[i+1][j][k][2] );
	frct[i][j][k][3] = frct[i][j][k][3]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][3] - flux[i][j][k][3] )
	  + dx4 * tx1 * (            rsd[i-1][j][k][3]
				     - two * rsd[i][j][k][3]
				     +           rsd[i+1][j][k][3] );
	frct[i][j][k][4] = frct[i][j][k][4]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][4] - flux[i][j][k][4] )
	  + dx5 * tx1 * (            rsd[i-1][j][k][4]
				     - two * rsd[i][j][k][4]
				     +           rsd[i+1][j][k][4] );
      }

/*--------------------------------------------------------------------
c   Fourth-order dissipation
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	frct[1][j][k][m] = frct[1][j][k][m]
	  - dsspm * ( + five * rsd[1][j][k][m]
		      - four * rsd[2][j][k][m]
		      +           rsd[3][j][k][m] );
	frct[2][j][k][m] = frct[2][j][k][m]
	  - dsspm * ( - four * rsd[1][j][k][m]
		      + six * rsd[2][j][k][m]
		      - four * rsd[3][j][k][m]
		      +           rsd[4][j][k][m] );
      }

      ist1 = 3;
      iend1 = nx - 4;
      for (i = ist1; i <=iend1; i++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = frct[i][j][k][m]
	    - dsspm * (            rsd[i-2][j][k][m]
				   - four * rsd[i-1][j][k][m]
				   + six * rsd[i][j][k][m]
				   - four * rsd[i+1][j][k][m]
				   +           rsd[i+2][j][k][m] );
	}
      }

      for (m = 0; m < 5; m++) {
	frct[nx-3][j][k][m] = frct[nx-3][j][k][m]
	  - dsspm * (             rsd[nx-5][j][k][m]
				  - four * rsd[nx-4][j][k][m]
				  + six * rsd[nx-3][j][k][m]
				  - four * rsd[nx-2][j][k][m]  );
	frct[nx-2][j][k][m] = frct[nx-2][j][k][m]
	  - dsspm * (             rsd[nx-4][j][k][m]
				  - four * rsd[nx-3][j][k][m]
				  + five * rsd[nx-2][j][k][m] );
      }
    }
  }

/*--------------------------------------------------------------------
c   eta-direction flux differences
--------------------------------------------------------------------*/

  L1 = 0;
  L2 = ny-1;

#pragma omp for
  for (i = ist; i <= iend; i++) {
    for (j = L1; j <= L2; j++) {
      for (k = 1; k <= nz - 2; k++) {
	flux[i][j][k][0] = rsd[i][j][k][2];
	u31 = rsd[i][j][k][2] / rsd[i][j][k][0];
	q = (one/two) * (  rsd[i][j][k][1] * rsd[i][j][k][1]
		      + rsd[i][j][k][2] * rsd[i][j][k][2]
		      + rsd[i][j][k][3] * rsd[i][j][k][3] )
	  / rsd[i][j][k][0];
	flux[i][j][k][1] = rsd[i][j][k][1] * u31;
	flux[i][j][k][2] = rsd[i][j][k][2] * u31 + C2 * 
	  ( rsd[i][j][k][4] - q );
	flux[i][j][k][3] = rsd[i][j][k][3] * u31;
	flux[i][j][k][4] = ( C1 * rsd[i][j][k][4] - C2 * q ) * u31;
      }
    }
  }

#pragma omp for
  for (i = ist; i <= iend; i++) {
    for (k = 1; k <= nz - 2; k++) {
      for (j = jst; j <= jend; j++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] =  frct[i][j][k][m]
	    - ty2 * ( flux[i][j+1][k][m] - flux[i][j-1][k][m] );
	}
      }
      for (j = jst; j <= L2; j++) {
	tmp = one / rsd[i][j][k][0];

	u21j = tmp * rsd[i][j][k][1];
	u31j = tmp * rsd[i][j][k][2];
	u41j = tmp * rsd[i][j][k][3];
	u51j = tmp * rsd[i][j][k][4];

	tmp = one / rsd[i][j-1][k][0];

	u21jm1 = tmp * rsd[i][j-1][k][1];
	u31jm1 = tmp * rsd[i][j-1][k][2];
	u41jm1 = tmp * rsd[i][j-1][k][3];
	u51jm1 = tmp * rsd[i][j-1][k][4];

	flux[i][j][k][1] = ty3 * ( u21j - u21jm1 );
	flux[i][j][k][2] = (four/three) * ty3 *
	  ( u31j - u31jm1 );
	flux[i][j][k][3] = ty3 * ( u41j - u41jm1 );
	flux[i][j][k][4] = (one/two) * ( one - C1*C5 )
	  * ty3 * ( ( u21j  *u21j + u31j  *u31j + u41j  *u41j )
		    - ( u21jm1*u21jm1 + u31jm1*u31jm1 + u41jm1*u41jm1 ) )
	  + (one/six)
	  * ty3 * ( u31j*u31j - u31jm1*u31jm1 )
	  + C1 * C5 * ty3 * ( u51j - u51jm1 );
      }

      for (j = jst; j <= jend; j++) {
	frct[i][j][k][0] = frct[i][j][k][0]
	  + dy1 * ty1 * (            rsd[i][j-1][k][0]
				     - two * rsd[i][j][k][0]
				     +           rsd[i][j+1][k][0] );
	frct[i][j][k][1] = frct[i][j][k][1]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][1] - flux[i][j][k][1] )
	  + dy2 * ty1 * (            rsd[i][j-1][k][1]
				     - two * rsd[i][j][k][1]
				     +           rsd[i][j+1][k][1] );
	frct[i][j][k][2] = frct[i][j][k][2]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][2] - flux[i][j][k][2] )
	  + dy3 * ty1 * (            rsd[i][j-1][k][2]
				     - two * rsd[i][j][k][2]
				     +           rsd[i][j+1][k][2] );
	frct[i][j][k][3] = frct[i][j][k][3]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][3] - flux[i][j][k][3] )
	  + dy4 * ty1 * (            rsd[i][j-1][k][3]
				     - two * rsd[i][j][k][3]
				     +           rsd[i][j+1][k][3] );
	frct[i][j][k][4] = frct[i][j][k][4]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][4] - flux[i][j][k][4] )
	  + dy5 * ty1 * (            rsd[i][j-1][k][4]
				     - two * rsd[i][j][k][4]
				     +           rsd[i][j+1][k][4] );
      }

/*--------------------------------------------------------------------
c   fourth-order dissipation
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	frct[i][1][k][m] = frct[i][1][k][m]
	  - dsspm * ( + five * rsd[i][1][k][m]
		      - four * rsd[i][2][k][m]
		      +           rsd[i][3][k][m] );
	frct[i][2][k][m] = frct[i][2][k][m]
	  - dsspm * ( - four * rsd[i][1][k][m]
		      + six * rsd[i][2][k][m]
		      - four * rsd[i][3][k][m]
		      +           rsd[i][4][k][m] );
      }

      jst1 = 3;
      jend1 = ny - 4;

      for (j = jst1; j <= jend1; j++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = frct[i][j][k][m]
	    - dsspm * (            rsd[i][j-2][k][m]
				   - four * rsd[i][j-1][k][m]
				   + six * rsd[i][j][k][m]
				   - four * rsd[i][j+1][k][m]
				   +           rsd[i][j+2][k][m] );
	}
      }

      for (m = 0; m < 5; m++) {
	frct[i][ny-3][k][m] = frct[i][ny-3][k][m]
	  - dsspm * (             rsd[i][ny-5][k][m]
				  - four * rsd[i][ny-4][k][m]
				  + six * rsd[i][ny-3][k][m]
				  - four * rsd[i][ny-2][k][m]  );
	frct[i][ny-2][k][m] = frct[i][ny-2][k][m]
	  - dsspm * (             rsd[i][ny-4][k][m]
				  - four * rsd[i][ny-3][k][m]
				  + five * rsd[i][ny-2][k][m]  );
      }

    }
  }

/*--------------------------------------------------------------------
c   zeta-direction flux differences
--------------------------------------------------------------------*/
#pragma omp for
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 0; k <= nz-1; k++) {
	flux[i][j][k][0] = rsd[i][j][k][3];
	u41 = rsd[i][j][k][3] / rsd[i][j][k][0];
	q = (one/two) * (  rsd[i][j][k][1] * rsd[i][j][k][1]
		      + rsd[i][j][k][2] * rsd[i][j][k][2]
		      + rsd[i][j][k][3] * rsd[i][j][k][3] )
	  / rsd[i][j][k][0];
	flux[i][j][k][1] = rsd[i][j][k][1] * u41;
	flux[i][j][k][2] = rsd[i][j][k][2] * u41;
	flux[i][j][k][3] = rsd[i][j][k][3] * u41 + C2 * 
	  ( rsd[i][j][k][4] - q );
	flux[i][j][k][4] = ( C1 * rsd[i][j][k][4] - C2 * q ) * u41;
      }

      for (k = 1; k <= nz - 2; k++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] =  frct[i][j][k][m]
	    - tz2 * ( flux[i][j][k+1][m] - flux[i][j][k-1][m] );
	}
      }
      for (k = 1; k <= nz-1; k++) {
	tmp = one / rsd[i][j][k][0];

	u21k = tmp * rsd[i][j][k][1];
	u31k = tmp * rsd[i][j][k][2];
	u41k = tmp * rsd[i][j][k][3];
	u51k = tmp * rsd[i][j][k][4];

	tmp = one / rsd[i][j][k-1][0];

	u21km1 = tmp * rsd[i][j][k-1][1];
	u31km1 = tmp * rsd[i][j][k-1][2];
	u41km1 = tmp * rsd[i][j][k-1][3];
	u51km1 = tmp * rsd[i][j][k-1][4];

	flux[i][j][k][1] = tz3 * ( u21k - u21km1 );
	flux[i][j][k][2] = tz3 * ( u31k - u31km1 );
	flux[i][j][k][3] = (four/three) * tz3 * ( u41k
					    - u41km1 );
	flux[i][j][k][4] = (one/two) * ( one - C1*C5 )
	  * tz3 * ( ( u21k  *u21k + u31k  *u31k + u41k  *u41k )
		    - ( u21km1*u21km1 + u31km1*u31km1 + u41km1*u41km1 ) )
	  + (one/six)
	  * tz3 * ( u41k*u41k - u41km1*u41km1 )
	  + C1 * C5 * tz3 * ( u51k - u51km1 );
      }

      for (k = 1; k <= nz - 2; k++) {
	frct[i][j][k][0] = frct[i][j][k][0]
	  + dz1 * tz1 * (            rsd[i][j][k+1][0]
				     - two * rsd[i][j][k][0]
				     +           rsd[i][j][k-1][0] );
	frct[i][j][k][1] = frct[i][j][k][1]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][1] - flux[i][j][k][1] )
	  + dz2 * tz1 * (            rsd[i][j][k+1][1]
				     - two * rsd[i][j][k][1]
				     +           rsd[i][j][k-1][1] );
	frct[i][j][k][2] = frct[i][j][k][2]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][2] - flux[i][j][k][2] )
	  + dz3 * tz1 * (            rsd[i][j][k+1][2]
				     - two * rsd[i][j][k][2]
				     +           rsd[i][j][k-1][2] );
	frct[i][j][k][3] = frct[i][j][k][3]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][3] - flux[i][j][k][3] )
	  + dz4 * tz1 * (            rsd[i][j][k+1][3]
				     - two * rsd[i][j][k][3]
				     +           rsd[i][j][k-1][3] );
	frct[i][j][k][4] = frct[i][j][k][4]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][4] - flux[i][j][k][4] )
	  + dz5 * tz1 * (            rsd[i][j][k+1][4]
				     - two * rsd[i][j][k][4]
				     +           rsd[i][j][k-1][4] );
      }

/*--------------------------------------------------------------------
c   fourth-order dissipation
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	frct[i][j][1][m] = frct[i][j][1][m]
	  - dsspm * ( + five * rsd[i][j][1][m]
		      - four * rsd[i][j][2][m]
		      +           rsd[i][j][3][m] );
	frct[i][j][2][m] = frct[i][j][2][m]
	  - dsspm * (- four * rsd[i][j][1][m]
		     + six * rsd[i][j][2][m]
		     - four * rsd[i][j][3][m]
		     +           rsd[i][j][4][m] );
      }

      for (k = 3; k <= nz - 4; k++) {
	for (m = 0; m < 5; m++) {
	  frct[i][j][k][m] = frct[i][j][k][m]
	    - dsspm * (           rsd[i][j][k-2][m]
				  - four * rsd[i][j][k-1][m]
				  + six * rsd[i][j][k][m]
				  - four * rsd[i][j][k+1][m]
				  +           rsd[i][j][k+2][m] );
	}
      }

      for (m = 0; m < 5; m++) {
	frct[i][j][nz-3][m] = frct[i][j][nz-3][m]
	  - dsspm * (            rsd[i][j][nz-5][m]
				 - four * rsd[i][j][nz-4][m]
				 + six * rsd[i][j][nz-3][m]
				 - four * rsd[i][j][nz-2][m]  );
        frct[i][j][nz-2][m] = frct[i][j][nz-2][m]
	  - dsspm * (             rsd[i][j][nz-4][m]
				  - four * rsd[i][j][nz-3][m]
				  + five * rsd[i][j][nz-2][m]  );
      }
    }
  }
}
}
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void error(void) {

/*--------------------------------------------------------------------
c
c   compute the solution error
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int iglob, jglob;
  element_t  tmp;
  element_t  u000ijk[5];

  for (m = 0; m < 5; m++) {
    errnm[m] = zero;
  }
  
  for (i = ist; i <= iend; i++) {
    iglob = i;
    for (j = jst; j <= jend; j++) {
      jglob = j;
      for (k = 1; k <= nz-2; k++) {
	exact( iglob, jglob, k, u000ijk );
	for (m = 0; m < 5; m++) {
	  tmp = ( u000ijk[m] - u[i][j][k][m] );
	  errnm[m] = errnm[m] + tmp *tmp;
	}
      }
    }
  }

  for (m = 0; m < 5; m++) {
    errnm[m] = sqrt_asm ( errnm[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
  }
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void exact( int i, int j, int k, element_t u000ijk[5] ) {

/*--------------------------------------------------------------------
c
c   compute the exact solution at (i,j,k)
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int m;
  element_t xi, eta, zeta;

  xi  = ((element_t)i) / (nx0 - 1);
  eta  = ((element_t)j) / (ny0 - 1);
  zeta = ((element_t)k) / (nz - 1);

  for (m = 0; m < 5; m++) {
    u000ijk[m] =  ce[m][0]
      + ce[m][1] * xi
      + ce[m][2] * eta
      + ce[m][3] * zeta
      + ce[m][4] * xi * xi
      + ce[m][5] * eta * eta
      + ce[m][6] * zeta * zeta
      + ce[m][7] * xi * xi * xi
      + ce[m][8] * eta * eta * eta
      + ce[m][9] * zeta * zeta * zeta
      + ce[m][10] * xi * xi * xi * xi
      + ce[m][11] * eta * eta * eta * eta
      + ce[m][12] * zeta * zeta * zeta * zeta;
  }
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void jacld(int k) {

/*--------------------------------------------------------------------
c   compute the lower triangular part of the jacobian matrix
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j;
  element_t  r43;
  element_t  c1345;
  element_t  c34;
  element_t  tmp1, tmp2, tmp3;

  r43 = ( four / three );
  c1345 = C1 * C3 * C4 * C5;
  c34 = C3 * C4;

#pragma omp for nowait schedule(static)
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {

/*--------------------------------------------------------------------
c   form the block daigonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      d[i][j][0][0] =  one
	+ dt * two * (   tx1 * dx1
			 + ty1 * dy1
			 + tz1 * dz1 );
      d[i][j][0][1] =  zero;
      d[i][j][0][2] =  zero;
      d[i][j][0][3] =  zero;
      d[i][j][0][4] =  zero;

      d[i][j][1][0] =  dt * two
	* (  tx1 * ( - r43 * c34 * tmp2 * u[i][j][k][1] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][1] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][1] ) );
      d[i][j][1][1] =  one
	+ dt * two
	* (  tx1 * r43 * c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * two * (   tx1 * dx2
			 + ty1 * dy2
			 + tz1 * dz2  );
      d[i][j][1][2] = zero;
      d[i][j][1][3] = zero;
      d[i][j][1][4] = zero;

      d[i][j][2][0] = dt * two
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][2] )
	     + ty1 * ( - r43 * c34 * tmp2 * u[i][j][k][2] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][2] ) );
      d[i][j][2][1] = zero;
      d[i][j][2][2] = one
	+ dt * two
	* (  tx1 *       c34 * tmp1
	     + ty1 * r43 * c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * two * (  tx1 * dx3
			+ ty1 * dy3
			+ tz1 * dz3 );
      d[i][j][2][3] = zero;
      d[i][j][2][4] = zero;

      d[i][j][3][0] = dt * two
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + tz1 * ( - r43 * c34 * tmp2 * u[i][j][k][3] ) );
      d[i][j][3][1] = zero;
      d[i][j][3][2] = zero;
      d[i][j][3][3] = one
	+ dt * two
	* (  tx1 *       c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 * r43 * c34 * tmp1 )
	+ dt * two * (  tx1 * dx4
			+ ty1 * dy4
			+ tz1 * dz4 );
      d[i][j][3][4] = zero;

      d[i][j][4][0] = dt * two
	* ( tx1 * ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		    - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + ty1 * ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		      - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + tz1 * ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		      - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] ) );
      d[i][j][4][1] = dt * two
	* ( tx1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + ty1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + tz1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1] );
      d[i][j][4][2] = dt * two
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2]
	    + ty1 * ( r43*c34 -c1345 ) * tmp2 * u[i][j][k][2]
	    + tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2] );
      d[i][j][4][3] = dt * two
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + ty1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][3] );
      d[i][j][4][4] = one
	+ dt * two * ( tx1 * c1345 * tmp1
		       + ty1 * c1345 * tmp1
		       + tz1 * c1345 * tmp1 )
        + dt * two * (  tx1 * dx5
			+  ty1 * dy5
			+  tz1 * dz5 );

/*--------------------------------------------------------------------
c   form the first block sub-diagonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i][j][k-1][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      a[i][j][0][0] = - dt * tz1 * dz1;
      a[i][j][0][1] =   zero;
      a[i][j][0][2] =   zero;
      a[i][j][0][3] = - dt * tz2;
      a[i][j][0][4] =   zero;

      a[i][j][1][0] = - dt * tz2
	* ( - ( u[i][j][k-1][1]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k-1][1] );
      a[i][j][1][1] = - dt * tz2 * ( u[i][j][k-1][3] * tmp1 )
	- dt * tz1 * c34 * tmp1
	- dt * tz1 * dz2 ;
      a[i][j][1][2] = zero;
      a[i][j][1][3] = - dt * tz2 * ( u[i][j][k-1][1] * tmp1 );
      a[i][j][1][4] = zero;

      a[i][j][2][0] = - dt * tz2
	* ( - ( u[i][j][k-1][2]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k-1][2] );
      a[i][j][2][1] = zero;
      a[i][j][2][2] = - dt * tz2 * ( u[i][j][k-1][3] * tmp1 )
	- dt * tz1 * ( c34 * tmp1 )
	- dt * tz1 * dz3;
      a[i][j][2][3] = - dt * tz2 * ( u[i][j][k-1][2] * tmp1 );
      a[i][j][2][4] = zero;

      a[i][j][3][0] = - dt * tz2
	* ( - ( u[i][j][k-1][3] * tmp1 ) *( u[i][j][k-1][3] * tmp1 )
	    + (one/two) * C2
	    * ( ( u[i][j][k-1][1] * u[i][j][k-1][1]
		  + u[i][j][k-1][2] * u[i][j][k-1][2]
		  + u[i][j][k-1][3] * u[i][j][k-1][3] ) * tmp2 ) )
	- dt * tz1 * ( - r43 * c34 * tmp2 * u[i][j][k-1][3] );
      a[i][j][3][1] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][1] * tmp1 ) );
      a[i][j][3][2] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][2] * tmp1 ) );
      a[i][j][3][3] = - dt * tz2 * ( two - C2 )
	* ( u[i][j][k-1][3] * tmp1 )
	- dt * tz1 * ( r43 * c34 * tmp1 )
	- dt * tz1 * dz4;
      a[i][j][3][4] = - dt * tz2 * C2;

      a[i][j][4][0] = - dt * tz2
	* ( ( C2 * (  u[i][j][k-1][1] * u[i][j][k-1][1]
                      + u[i][j][k-1][2] * u[i][j][k-1][2]
                      + u[i][j][k-1][3] * u[i][j][k-1][3] ) * tmp2
	      - C1 * ( u[i][j][k-1][4] * tmp1 ) )
	    * ( u[i][j][k-1][3] * tmp1 ) )
	- dt * tz1
	* ( - ( c34 - c1345 ) * tmp3 * (u[i][j][k-1][1]*u[i][j][k-1][1])
	    - ( c34 - c1345 ) * tmp3 * (u[i][j][k-1][2]*u[i][j][k-1][2])
	    - ( r43*c34 - c1345 )* tmp3 * (u[i][j][k-1][3]*u[i][j][k-1][3])
	    - c1345 * tmp2 * u[i][j][k-1][4] );
      a[i][j][4][1] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][1]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k-1][1];
      a[i][j][4][2] = - dt * tz2
	* ( - C2 * ( u[i][j][k-1][2]*u[i][j][k-1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k-1][2];
      a[i][j][4][3] = - dt * tz2
	* ( C1 * ( u[i][j][k-1][4] * tmp1 )
            - (one/two) * C2
            * ( (  u[i][j][k-1][1]*u[i][j][k-1][1]
		   + u[i][j][k-1][2]*u[i][j][k-1][2]
		   + three*u[i][j][k-1][3]*u[i][j][k-1][3] ) * tmp2 ) )
	- dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k-1][3];
      a[i][j][4][4] = - dt * tz2
	* ( C1 * ( u[i][j][k-1][3] * tmp1 ) )
	- dt * tz1 * c1345 * tmp1
	- dt * tz1 * dz5;

/*--------------------------------------------------------------------
c   form the second block sub-diagonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i][j-1][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      b[i][j][0][0] = - dt * ty1 * dy1;
      b[i][j][0][1] =   zero;
      b[i][j][0][2] = - dt * ty2;
      b[i][j][0][3] =   zero;
      b[i][j][0][4] =   zero;

      b[i][j][1][0] = - dt * ty2
	* ( - ( u[i][j-1][k][1]*u[i][j-1][k][2] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j-1][k][1] );
      b[i][j][1][1] = - dt * ty2 * ( u[i][j-1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy2;
      b[i][j][1][2] = - dt * ty2 * ( u[i][j-1][k][1] * tmp1 );
      b[i][j][1][3] = zero;
      b[i][j][1][4] = zero;

      b[i][j][2][0] = - dt * ty2
	* ( - ( u[i][j-1][k][2] * tmp1 ) *( u[i][j-1][k][2] * tmp1 )
	    + one/two * C2 * ( (  u[i][j-1][k][1] * u[i][j-1][k][1]
			       + u[i][j-1][k][2] * u[i][j-1][k][2]
			       + u[i][j-1][k][3] * u[i][j-1][k][3] )
			    * tmp2 ) )
	- dt * ty1 * ( - r43 * c34 * tmp2 * u[i][j-1][k][2] );
      b[i][j][2][1] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][1] * tmp1 ) );
      b[i][j][2][2] = - dt * ty2 * ( ( two - C2 )
				  * ( u[i][j-1][k][2] * tmp1 ) )
	- dt * ty1 * ( r43 * c34 * tmp1 )
	- dt * ty1 * dy3;
      b[i][j][2][3] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][3] * tmp1 ) );
      b[i][j][2][4] = - dt * ty2 * C2;

      b[i][j][3][0] = - dt * ty2
	* ( - ( u[i][j-1][k][2]*u[i][j-1][k][3] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j-1][k][3] );
      b[i][j][3][1] = zero;
      b[i][j][3][2] = - dt * ty2 * ( u[i][j-1][k][3] * tmp1 );
      b[i][j][3][3] = - dt * ty2 * ( u[i][j-1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy4;
      b[i][j][3][4] = zero;

      b[i][j][4][0] = - dt * ty2
	* ( ( C2 * (  u[i][j-1][k][1] * u[i][j-1][k][1]
		      + u[i][j-1][k][2] * u[i][j-1][k][2]
		      + u[i][j-1][k][3] * u[i][j-1][k][3] ) * tmp2
	      - C1 * ( u[i][j-1][k][4] * tmp1 ) )
	    * ( u[i][j-1][k][2] * tmp1 ) )
	- dt * ty1
	* ( - (     c34 - c1345 )*tmp3*(pow2(u[i][j-1][k][1]))
	    - ( r43*c34 - c1345 )*tmp3*(pow2(u[i][j-1][k][2]))
	    - (     c34 - c1345 )*tmp3*(pow2(u[i][j-1][k][3]))
	    - c1345*tmp2*u[i][j-1][k][4] );
      b[i][j][4][1] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][1]*u[i][j-1][k][2] ) * tmp2 )
	- dt * ty1
	* ( c34 - c1345 ) * tmp2 * u[i][j-1][k][1];
      b[i][j][4][2] = - dt * ty2
	* ( C1 * ( u[i][j-1][k][4] * tmp1 )
	    - (one/two) * C2
	    * ( (  u[i][j-1][k][1]*u[i][j-1][k][1]
                   + three * u[i][j-1][k][2]*u[i][j-1][k][2]
		   + u[i][j-1][k][3]*u[i][j-1][k][3] ) * tmp2 ) )
	- dt * ty1
	* ( r43*c34 - c1345 ) * tmp2 * u[i][j-1][k][2];
      b[i][j][4][3] = - dt * ty2
	* ( - C2 * ( u[i][j-1][k][2]*u[i][j-1][k][3] ) * tmp2 )
	- dt * ty1 * ( c34 - c1345 ) * tmp2 * u[i][j-1][k][3];
      b[i][j][4][4] = - dt * ty2
	* ( C1 * ( u[i][j-1][k][2] * tmp1 ) )
	- dt * ty1 * c1345 * tmp1
	- dt * ty1 * dy5;

/*--------------------------------------------------------------------
c   form the third block sub-diagonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i-1][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      c[i][j][0][0] = - dt * tx1 * dx1;
      c[i][j][0][1] = - dt * tx2;
      c[i][j][0][2] =   zero;
      c[i][j][0][3] =   zero;
      c[i][j][0][4] =   zero;

      c[i][j][1][0] = - dt * tx2
	* ( - ( u[i-1][j][k][1] * tmp1 ) *( u[i-1][j][k][1] * tmp1 )
	    + C2 * (one/two) * (  u[i-1][j][k][1] * u[i-1][j][k][1]
                             + u[i-1][j][k][2] * u[i-1][j][k][2]
                             + u[i-1][j][k][3] * u[i-1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - r43 * c34 * tmp2 * u[i-1][j][k][1] );
      c[i][j][1][1] = - dt * tx2
	* ( ( two - C2 ) * ( u[i-1][j][k][1] * tmp1 ) )
	- dt * tx1 * ( r43 * c34 * tmp1 )
	- dt * tx1 * dx2;
      c[i][j][1][2] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][2] * tmp1 ) );
      c[i][j][1][3] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][3] * tmp1 ) );
      c[i][j][1][4] = - dt * tx2 * C2;

      c[i][j][2][0] = - dt * tx2
	* ( - ( u[i-1][j][k][1] * u[i-1][j][k][2] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i-1][j][k][2] );
      c[i][j][2][1] = - dt * tx2 * ( u[i-1][j][k][2] * tmp1 );
      c[i][j][2][2] = - dt * tx2 * ( u[i-1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx3;
      c[i][j][2][3] = zero;
      c[i][j][2][4] = zero;

      c[i][j][3][0] = - dt * tx2
	* ( - ( u[i-1][j][k][1]*u[i-1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i-1][j][k][3] );
      c[i][j][3][1] = - dt * tx2 * ( u[i-1][j][k][3] * tmp1 );
      c[i][j][3][2] = zero;
      c[i][j][3][3] = - dt * tx2 * ( u[i-1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx4;
      c[i][j][3][4] = zero;

      c[i][j][4][0] = - dt * tx2
	* ( ( C2 * (  u[i-1][j][k][1] * u[i-1][j][k][1]
		      + u[i-1][j][k][2] * u[i-1][j][k][2]
		      + u[i-1][j][k][3] * u[i-1][j][k][3] ) * tmp2
	      - C1 * ( u[i-1][j][k][4] * tmp1 ) )
	    * ( u[i-1][j][k][1] * tmp1 ) )
	- dt * tx1
	* ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i-1][j][k][1]) )
	    - (     c34 - c1345 ) * tmp3 * ( pow2(u[i-1][j][k][2]) )
	    - (     c34 - c1345 ) * tmp3 * ( pow2(u[i-1][j][k][3]) )
	    - c1345 * tmp2 * u[i-1][j][k][4] );
      c[i][j][4][1] = - dt * tx2
	* ( C1 * ( u[i-1][j][k][4] * tmp1 )
	    - (one/two) * C2
	    * ( (  three*u[i-1][j][k][1]*u[i-1][j][k][1]
		   + u[i-1][j][k][2]*u[i-1][j][k][2]
		   + u[i-1][j][k][3]*u[i-1][j][k][3] ) * tmp2 ) )
	- dt * tx1
	* ( r43*c34 - c1345 ) * tmp2 * u[i-1][j][k][1];
      c[i][j][4][2] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][2]*u[i-1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i-1][j][k][2];
      c[i][j][4][3] = - dt * tx2
	* ( - C2 * ( u[i-1][j][k][3]*u[i-1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i-1][j][k][3];
      c[i][j][4][4] = - dt * tx2
	* ( C1 * ( u[i-1][j][k][1] * tmp1 ) )
	- dt * tx1 * c1345 * tmp1
	- dt * tx1 * dx5;
    }
  }
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void jacu(int k) {

/*--------------------------------------------------------------------
c   compute the upper triangular part of the jacobian matrix
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j;
  element_t  r43;
  element_t  c1345;
  element_t  c34;
  element_t  tmp1, tmp2, tmp3;

  r43 = ( four/ three );
  c1345 = C1 * C3 * C4 * C5;
  c34 = C3 * C4;

#pragma omp for nowait schedule(static)
#if defined(_OPENMP)  
  for (i = iend; i >= ist; i--) {
      for (j = jend; j >= jst; j--) {
#else	  
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
#endif	

/*--------------------------------------------------------------------
c   form the block daigonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      d[i][j][0][0] =  one
	+ dt * two * (   tx1 * dx1
			 + ty1 * dy1
			 + tz1 * dz1 );
      d[i][j][0][1] =  zero;
      d[i][j][0][2] =  zero;
      d[i][j][0][3] =  zero;
      d[i][j][0][4] =  zero;

      d[i][j][1][0] =  dt * two
	* (  tx1 * ( - r43 * c34 * tmp2 * u[i][j][k][1] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][1] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][1] ) );
      d[i][j][1][1] =  one
	+ dt * two
	* (  tx1 * r43 * c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * two * (   tx1 * dx2
			 + ty1 * dy2
			 + tz1 * dz2  );
      d[i][j][1][2] = zero;
      d[i][j][1][3] = zero;
      d[i][j][1][4] = zero;

      d[i][j][2][0] = dt * two
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][2] )
	     + ty1 * ( - r43 * c34 * tmp2 * u[i][j][k][2] )
	     + tz1 * ( -       c34 * tmp2 * u[i][j][k][2] ) );
      d[i][j][2][1] = zero;
      d[i][j][2][2] = one
	+ dt * two
	* (  tx1 *       c34 * tmp1
	     + ty1 * r43 * c34 * tmp1
	     + tz1 *       c34 * tmp1 )
	+ dt * two * (  tx1 * dx3
			+ ty1 * dy3
			+ tz1 * dz3 );
      d[i][j][2][3] = zero;
      d[i][j][2][4] = zero;

      d[i][j][3][0] = dt * two
	* (  tx1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + ty1 * ( -       c34 * tmp2 * u[i][j][k][3] )
	     + tz1 * ( - r43 * c34 * tmp2 * u[i][j][k][3] ) );
      d[i][j][3][1] = zero;
      d[i][j][3][2] = zero;
      d[i][j][3][3] = one
	+ dt * two
	* (  tx1 *       c34 * tmp1
	     + ty1 *       c34 * tmp1
	     + tz1 * r43 * c34 * tmp1 )
	+ dt * two * (  tx1 * dx4
			+ ty1 * dy4
			+ tz1 * dz4 );
      d[i][j][3][4] = zero;

      d[i][j][4][0] = dt * two
	* ( tx1 * ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		    - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + ty1 * ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		      - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] )
	    + tz1 * ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][1]) )
		      - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][2]) )
		      - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k][3]) )
		      - ( c1345 ) * tmp2 * u[i][j][k][4] ) );
      d[i][j][4][1] = dt * two
	* ( tx1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + ty1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1]
	    + tz1 * (     c34 - c1345 ) * tmp2 * u[i][j][k][1] );
      d[i][j][4][2] = dt * two
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2]
	    + ty1 * ( r43*c34 -c1345 ) * tmp2 * u[i][j][k][2]
	    + tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][2] );
      d[i][j][4][3] = dt * two
	* ( tx1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + ty1 * ( c34 - c1345 ) * tmp2 * u[i][j][k][3]
	    + tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k][3] );
      d[i][j][4][4] = one
        + dt * two * ( tx1 * c1345 * tmp1
		       + ty1 * c1345 * tmp1
		       + tz1 * c1345 * tmp1 )
        + dt * two * (  tx1 * dx5
			+  ty1 * dy5
			+  tz1 * dz5 );

/*--------------------------------------------------------------------
c   form the first block sub-diagonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i+1][j][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      a[i][j][0][0] = - dt * tx1 * dx1;
      a[i][j][0][1] =   dt * tx2;
      a[i][j][0][2] =   zero;
      a[i][j][0][3] =   zero;
      a[i][j][0][4] =   zero;

      a[i][j][1][0] =  dt * tx2
	* ( - ( u[i+1][j][k][1] * tmp1 ) *( u[i+1][j][k][1] * tmp1 )
	    + C2 * (one/two) * (  u[i+1][j][k][1] * u[i+1][j][k][1]
                             + u[i+1][j][k][2] * u[i+1][j][k][2]
                             + u[i+1][j][k][3] * u[i+1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - r43 * c34 * tmp2 * u[i+1][j][k][1] );
      a[i][j][1][1] =  dt * tx2
	* ( ( two - C2 ) * ( u[i+1][j][k][1] * tmp1 ) )
	- dt * tx1 * ( r43 * c34 * tmp1 )
	- dt * tx1 * dx2;
      a[i][j][1][2] =  dt * tx2
	* ( - C2 * ( u[i+1][j][k][2] * tmp1 ) );
      a[i][j][1][3] =  dt * tx2
	* ( - C2 * ( u[i+1][j][k][3] * tmp1 ) );
      a[i][j][1][4] =  dt * tx2 * C2 ;

      a[i][j][2][0] =  dt * tx2
	* ( - ( u[i+1][j][k][1] * u[i+1][j][k][2] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i+1][j][k][2] );
      a[i][j][2][1] =  dt * tx2 * ( u[i+1][j][k][2] * tmp1 );
      a[i][j][2][2] =  dt * tx2 * ( u[i+1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx3;
      a[i][j][2][3] = zero;
      a[i][j][2][4] = zero;

      a[i][j][3][0] = dt * tx2
	* ( - ( u[i+1][j][k][1]*u[i+1][j][k][3] ) * tmp2 )
	- dt * tx1 * ( - c34 * tmp2 * u[i+1][j][k][3] );
      a[i][j][3][1] = dt * tx2 * ( u[i+1][j][k][3] * tmp1 );
      a[i][j][3][2] = zero;
      a[i][j][3][3] = dt * tx2 * ( u[i+1][j][k][1] * tmp1 )
	- dt * tx1 * ( c34 * tmp1 )
	- dt * tx1 * dx4;
      a[i][j][3][4] = zero;

      a[i][j][4][0] = dt * tx2
	* ( ( C2 * (  u[i+1][j][k][1] * u[i+1][j][k][1]
		      + u[i+1][j][k][2] * u[i+1][j][k][2]
		      + u[i+1][j][k][3] * u[i+1][j][k][3] ) * tmp2
	      - C1 * ( u[i+1][j][k][4] * tmp1 ) )
	    * ( u[i+1][j][k][1] * tmp1 ) )
	- dt * tx1
	* ( - ( r43*c34 - c1345 ) * tmp3 * ( pow2(u[i+1][j][k][1]) )
	    - (     c34 - c1345 ) * tmp3 * ( pow2(u[i+1][j][k][2]) )
	    - (     c34 - c1345 ) * tmp3 * ( pow2(u[i+1][j][k][3]) )
	    - c1345 * tmp2 * u[i+1][j][k][4] );
      a[i][j][4][1] = dt * tx2
	* ( C1 * ( u[i+1][j][k][4] * tmp1 )
	    - (one/two) * C2
	    * ( (  three*u[i+1][j][k][1]*u[i+1][j][k][1]
		   + u[i+1][j][k][2]*u[i+1][j][k][2]
		   + u[i+1][j][k][3]*u[i+1][j][k][3] ) * tmp2 ) )
	- dt * tx1
	* ( r43*c34 - c1345 ) * tmp2 * u[i+1][j][k][1];
      a[i][j][4][2] = dt * tx2
	* ( - C2 * ( u[i+1][j][k][2]*u[i+1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i+1][j][k][2];
      a[i][j][4][3] = dt * tx2
	* ( - C2 * ( u[i+1][j][k][3]*u[i+1][j][k][1] ) * tmp2 )
	- dt * tx1
	* (  c34 - c1345 ) * tmp2 * u[i+1][j][k][3];
      a[i][j][4][4] = dt * tx2
	* ( C1 * ( u[i+1][j][k][1] * tmp1 ) )
	- dt * tx1 * c1345 * tmp1
	- dt * tx1 * dx5;

/*--------------------------------------------------------------------
c   form the second block sub-diagonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i][j+1][k][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      b[i][j][0][0] = - dt * ty1 * dy1;
      b[i][j][0][1] =   zero;
      b[i][j][0][2] =  dt * ty2;
      b[i][j][0][3] =   zero;
      b[i][j][0][4] =   zero;

      b[i][j][1][0] =  dt * ty2
	* ( - ( u[i][j+1][k][1]*u[i][j+1][k][2] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j+1][k][1] );
      b[i][j][1][1] =  dt * ty2 * ( u[i][j+1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy2;
      b[i][j][1][2] =  dt * ty2 * ( u[i][j+1][k][1] * tmp1 );
      b[i][j][1][3] = zero;
      b[i][j][1][4] = zero;

      b[i][j][2][0] =  dt * ty2
	* ( - ( u[i][j+1][k][2] * tmp1 ) *( u[i][j+1][k][2] * tmp1 )
	    + (one/two) * C2 * ( (  u[i][j+1][k][1] * u[i][j+1][k][1]
			       + u[i][j+1][k][2] * u[i][j+1][k][2]
			       + u[i][j+1][k][3] * u[i][j+1][k][3] )
			    * tmp2 ) )
	- dt * ty1 * ( - r43 * c34 * tmp2 * u[i][j+1][k][2] );
      b[i][j][2][1] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][1] * tmp1 ) );
      b[i][j][2][2] =  dt * ty2 * ( ( two - C2 )
				 * ( u[i][j+1][k][2] * tmp1 ) )
	- dt * ty1 * ( r43 * c34 * tmp1 )
	- dt * ty1 * dy3;
      b[i][j][2][3] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][3] * tmp1 ) );
      b[i][j][2][4] =  dt * ty2 * C2;

      b[i][j][3][0] =  dt * ty2
	* ( - ( u[i][j+1][k][2]*u[i][j+1][k][3] ) * tmp2 )
	- dt * ty1 * ( - c34 * tmp2 * u[i][j+1][k][3] );
      b[i][j][3][1] = zero;
      b[i][j][3][2] =  dt * ty2 * ( u[i][j+1][k][3] * tmp1 );
      b[i][j][3][3] =  dt * ty2 * ( u[i][j+1][k][2] * tmp1 )
	- dt * ty1 * ( c34 * tmp1 )
	- dt * ty1 * dy4;
      b[i][j][3][4] = zero;

      b[i][j][4][0] =  dt * ty2
	* ( ( C2 * (  u[i][j+1][k][1] * u[i][j+1][k][1]
		      + u[i][j+1][k][2] * u[i][j+1][k][2]
		      + u[i][j+1][k][3] * u[i][j+1][k][3] ) * tmp2
	      - C1 * ( u[i][j+1][k][4] * tmp1 ) )
	    * ( u[i][j+1][k][2] * tmp1 ) )
	- dt * ty1
	* ( - (     c34 - c1345 )*tmp3*( pow2(u[i][j+1][k][1]) )
	    - ( r43*c34 - c1345 )*tmp3*( pow2(u[i][j+1][k][2]) )
	    - (     c34 - c1345 )*tmp3*( pow2(u[i][j+1][k][3]) )
	    - c1345*tmp2*u[i][j+1][k][4] );
      b[i][j][4][1] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][1]*u[i][j+1][k][2] ) * tmp2 )
	- dt * ty1
	* ( c34 - c1345 ) * tmp2 * u[i][j+1][k][1];
      b[i][j][4][2] =  dt * ty2
	* ( C1 * ( u[i][j+1][k][4] * tmp1 )
	    - (one/two) * C2
	    * ( (  u[i][j+1][k][1]*u[i][j+1][k][1]
		   + three * u[i][j+1][k][2]*u[i][j+1][k][2]
		   + u[i][j+1][k][3]*u[i][j+1][k][3] ) * tmp2 ) )
	- dt * ty1
	* ( r43*c34 - c1345 ) * tmp2 * u[i][j+1][k][2];
      b[i][j][4][3] =  dt * ty2
	* ( - C2 * ( u[i][j+1][k][2]*u[i][j+1][k][3] ) * tmp2 )
	- dt * ty1 * ( c34 - c1345 ) * tmp2 * u[i][j+1][k][3];
      b[i][j][4][4] =  dt * ty2
	* ( C1 * ( u[i][j+1][k][2] * tmp1 ) )
	- dt * ty1 * c1345 * tmp1
	- dt * ty1 * dy5;

/*--------------------------------------------------------------------
c   form the third block sub-diagonal
--------------------------------------------------------------------*/
      tmp1 = one / u[i][j][k+1][0];
      tmp2 = tmp1 * tmp1;
      tmp3 = tmp1 * tmp2;

      c[i][j][0][0] = - dt * tz1 * dz1;
      c[i][j][0][1] =   zero;
      c[i][j][0][2] =   zero;
      c[i][j][0][3] = dt * tz2;
      c[i][j][0][4] =   zero;

      c[i][j][1][0] = dt * tz2
	* ( - ( u[i][j][k+1][1]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k+1][1] );
      c[i][j][1][1] = dt * tz2 * ( u[i][j][k+1][3] * tmp1 )
	- dt * tz1 * c34 * tmp1
	- dt * tz1 * dz2 ;
      c[i][j][1][2] = zero;
      c[i][j][1][3] = dt * tz2 * ( u[i][j][k+1][1] * tmp1 );
      c[i][j][1][4] = zero;

      c[i][j][2][0] = dt * tz2
	* ( - ( u[i][j][k+1][2]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( - c34 * tmp2 * u[i][j][k+1][2] );
      c[i][j][2][1] = zero;
      c[i][j][2][2] = dt * tz2 * ( u[i][j][k+1][3] * tmp1 )
	- dt * tz1 * ( c34 * tmp1 )
	- dt * tz1 * dz3;
      c[i][j][2][3] = dt * tz2 * ( u[i][j][k+1][2] * tmp1 );
      c[i][j][2][4] = zero;

      c[i][j][3][0] = dt * tz2
	* ( - ( u[i][j][k+1][3] * tmp1 ) *( u[i][j][k+1][3] * tmp1 )
	    + (one/two) * C2
	    * ( ( u[i][j][k+1][1] * u[i][j][k+1][1]
		  + u[i][j][k+1][2] * u[i][j][k+1][2]
		  + u[i][j][k+1][3] * u[i][j][k+1][3] ) * tmp2 ) )
	- dt * tz1 * ( - r43 * c34 * tmp2 * u[i][j][k+1][3] );
      c[i][j][3][1] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][1] * tmp1 ) );
      c[i][j][3][2] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][2] * tmp1 ) );
      c[i][j][3][3] = dt * tz2 * ( two - C2 )
	* ( u[i][j][k+1][3] * tmp1 )
	- dt * tz1 * ( r43 * c34 * tmp1 )
	- dt * tz1 * dz4;
      c[i][j][3][4] = dt * tz2 * C2;

      c[i][j][4][0] = dt * tz2
	* ( ( C2 * (  u[i][j][k+1][1] * u[i][j][k+1][1]
                      + u[i][j][k+1][2] * u[i][j][k+1][2]
                      + u[i][j][k+1][3] * u[i][j][k+1][3] ) * tmp2
	      - C1 * ( u[i][j][k+1][4] * tmp1 ) )
	    * ( u[i][j][k+1][3] * tmp1 ) )
	- dt * tz1
	* ( - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k+1][1]) )
	    - ( c34 - c1345 ) * tmp3 * ( pow2(u[i][j][k+1][2]) )
	    - ( r43*c34 - c1345 )* tmp3 * ( pow2(u[i][j][k+1][3]) )
	    - c1345 * tmp2 * u[i][j][k+1][4] );
      c[i][j][4][1] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][1]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k+1][1];
      c[i][j][4][2] = dt * tz2
	* ( - C2 * ( u[i][j][k+1][2]*u[i][j][k+1][3] ) * tmp2 )
	- dt * tz1 * ( c34 - c1345 ) * tmp2 * u[i][j][k+1][2];
      c[i][j][4][3] = dt * tz2
	* ( C1 * ( u[i][j][k+1][4] * tmp1 )
            - one/two * C2
            * ( (  u[i][j][k+1][1]*u[i][j][k+1][1]
		   + u[i][j][k+1][2]*u[i][j][k+1][2]
		   + three*u[i][j][k+1][3]*u[i][j][k+1][3] ) * tmp2 ) )
	- dt * tz1 * ( r43*c34 - c1345 ) * tmp2 * u[i][j][k+1][3];
      c[i][j][4][4] = dt * tz2
	* ( C1 * ( u[i][j][k+1][3] * tmp1 ) )
	- dt * tz1 * c1345 * tmp1
	- dt * tz1 * dz5;
    }
  }
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void l2norm (int nx0, int ny0, int nz0,
		    int ist, int iend,
		    int jst, int jend,
/*--------------------------------------------------------------------
c   To improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.  Only needed in v.
--------------------------------------------------------------------*/
		    element_t v[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5],
		    element_t sum[5]) {

#pragma omp parallel 
{

/*--------------------------------------------------------------------
c   to compute the l2-norm of vector v.
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  element_t sum0=zero, sum1=zero, sum2=zero, sum3=zero, sum4=zero;

#pragma omp single  
  for (m = 0; m < 5; m++) {
    sum[m] = zero;
  }

#pragma omp for nowait
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz0-2; k++) {
	  sum0 = sum0 + v[i][j][k][0] * v[i][j][k][0];
	  sum1 = sum1 + v[i][j][k][1] * v[i][j][k][1];
	  sum2 = sum2 + v[i][j][k][2] * v[i][j][k][2];
	  sum3 = sum3 + v[i][j][k][3] * v[i][j][k][3];
	  sum4 = sum4 + v[i][j][k][4] * v[i][j][k][4];
      }
    }
  }

#pragma omp critical
  {
      sum[0] += sum0;
      sum[1] += sum1;
      sum[2] += sum2;
      sum[3] += sum3;
      sum[4] += sum4;
  }
#pragma omp barrier  
  
#pragma omp single  
  for (m = 0;  m < 5; m++) {
    sum[m] = sqrt_asm ( sum[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
  }
}
}
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void pintgr(void) {

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k;
  int ibeg, ifin, ifin1;
  int jbeg, jfin, jfin1;
  int iglob, iglob1, iglob2;
  int jglob, jglob1, jglob2;
  element_t phi1[ISIZ2+2][ISIZ3+2];	/* phi1(0:isiz2+1,0:isiz3+1) */
  element_t phi2[ISIZ2+2][ISIZ3+2];	/* phi2(0:isiz2+1,0:isiz3+1) */
  element_t  frc1, frc2, frc3;

/*--------------------------------------------------------------------
c   set up the sub-domains for integeration in each processor
--------------------------------------------------------------------*/
  ibeg = nx;
  ifin = 0;
  iglob1 = -1;
  iglob2 = nx-1;
  if (iglob1 >= ii1 && iglob2 < ii2+nx) ibeg = 0;
  if (iglob1 >= ii1-nx && iglob2 <= ii2) ifin = nx;
  if (ii1 >= iglob1 && ii1 <= iglob2) ibeg = ii1;
  if (ii2 >= iglob1 && ii2 <= iglob2) ifin = ii2;
  jbeg = ny;
  jfin = -1;
  jglob1 = 0;
  jglob2 = ny-1;
  if (jglob1 >= ji1 && jglob2 < ji2+ny) jbeg = 0;
  if (jglob1 > ji1-ny && jglob2 <= ji2) jfin = ny;
  if (ji1 >= jglob1 && ji1 <= jglob2) jbeg = ji1;
  if (ji2 >= jglob1 && ji2 <= jglob2) jfin = ji2;
  ifin1 = ifin;
  jfin1 = jfin;
  if (ifin1 == ii2) ifin1 = ifin -1;
  if (jfin1 == ji2) jfin1 = jfin -1;

/*--------------------------------------------------------------------
c   initialize
--------------------------------------------------------------------*/
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = zero;
      phi2[i][k] = zero;
    }
  }
  for (i = ibeg; i <= ifin; i++) {
    iglob = i;
    for (j = jbeg; j <= jfin; j++) {
      jglob = j;

      k = ki1;

      phi1[i][j] = C2*(  u[i][j][k][4]
			- (one/two) * (  pow2(u[i][j][k][1])
				    + pow2(u[i][j][k][2])
				    + pow2(u[i][j][k][3]) )
			/ u[i][j][k][0] );

      k = ki2;

      phi2[i][j] = C2*(  u[i][j][k][4]
			- (one/two) * (  pow2(u[i][j][k][1])
				    + pow2(u[i][j][k][2])
				    + pow2(u[i][j][k][3]) )
			/ u[i][j][k][0] );
    }
  }

  frc1 = zero;

  for (i = ibeg; i <= ifin1; i++) {
    for (j = jbeg; j <= jfin1; j++) {
      frc1 = frc1 + (  phi1[i][j]
		       + phi1[i+1][j]
		       + phi1[i][j+1]
		       + phi1[i+1][j+1]
		       + phi2[i][j]
		       + phi2[i+1][j]
		       + phi2[i][j+1]
		       + phi2[i+1][j+1] );
    }
  }

  frc1 = dxi * deta * frc1;

/*--------------------------------------------------------------------
c   initialize
--------------------------------------------------------------------*/
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = zero;
      phi2[i][k] = zero;
    }
  }
  jglob = jbeg;
  if (jglob == ji1) {
    for (i = ibeg; i <= ifin; i++) {
      iglob = i;
      for (k = ki1; k <= ki2; k++) {
	phi1[i][k] = C2*(  u[i][jbeg][k][4]
			  - (one/two) * (  pow2(u[i][jbeg][k][1])
				      + pow2(u[i][jbeg][k][2])
				      + pow2(u[i][jbeg][k][3]) )
			  / u[i][jbeg][k][0] );
      }
    }
  }

  jglob = jfin;
  if (jglob == ji2) {
    for (i = ibeg; i <= ifin; i++) {
      iglob = i;
      for (k = ki1; k <= ki2; k++) {
	phi2[i][k] = C2*(  u[i][jfin][k][4]
			  - (one/two) * (  pow2(u[i][jfin][k][1])
				      + pow2(u[i][jfin][k][2])
				      + pow2(u[i][jfin][k][3]) )
			  / u[i][jfin][k][0] );
      }
    }
  }


  frc2 = zero;
  for (i = ibeg; i <= ifin1; i++) {
    for (k = ki1; k <= ki2-1; k++) {
      frc2 = frc2 + (  phi1[i][k]
		       + phi1[i+1][k]
		       + phi1[i][k+1]
		       + phi1[i+1][k+1]
		       + phi2[i][k]
		       + phi2[i+1][k]
		       + phi2[i][k+1]
		       + phi2[i+1][k+1] );
    }
  }


  frc2 = dxi * dzeta * frc2;

/*--------------------------------------------------------------------
c   initialize
--------------------------------------------------------------------*/
  for (i = 0; i <= ISIZ2+1; i++) {
    for (k = 0; k <= ISIZ3+1; k++) {
      phi1[i][k] = zero;
      phi2[i][k] = zero;
    }
  }
  iglob = ibeg;
  if (iglob == ii1) {
    for (j = jbeg; j <= jfin; j++) {
      jglob = j;
      for (k = ki1; k <= ki2; k++) {
	phi1[j][k] = C2*(  u[ibeg][j][k][4]
			  - (one/two) * (  pow2(u[ibeg][j][k][1])
				      + pow2(u[ibeg][j][k][2])
				      + pow2(u[ibeg][j][k][3]) )
			  / u[ibeg][j][k][0] );
      }
    }
  }

  iglob = ifin;
  if (iglob == ii2) {
    for (j = jbeg; j <= jfin; j++) {
      jglob = j;
      for (k = ki1; k <= ki2; k++) {
	phi2[j][k] = C2*(  u[ifin][j][k][4]
			  - (one/two) * (  pow2(u[ifin][j][k][1])
				      + pow2(u[ifin][j][k][2])
				      + pow2(u[ifin][j][k][3]) )
			  / u[ifin][j][k][0] );
      }
    }
  }

  frc3 = zero;

  for (j = jbeg; j <= jfin1; j++) {
    for (k = ki1; k <= ki2-1; k++) {
      frc3 = frc3 + (  phi1[j][k]
		       + phi1[j+1][k]
		       + phi1[j][k+1]
		       + phi1[j+1][k+1]
		       + phi2[j][k]
		       + phi2[j+1][k]
		       + phi2[j][k+1]
		       + phi2[j+1][k+1] );
    }
  }

  frc3 = deta * dzeta * frc3;
  frc = (one/four) * ( frc1 + frc2 + frc3 );
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static int read_input(void) {
    
/*--------------------------------------------------------------------
c    if input file does not exist, it uses defaults
c       ipr = 1 for detailed progress output
c       inorm = how often the norm is printed (once every inorm iterations)
c       itmax = number of pseudo time steps
c       dt = time step
c       omega 1 over-relaxation factor for SSOR
c       tolrsd = steady state residual tolerance levels
c       nx, ny, nz = number of grid points in x, y, z directions
--------------------------------------------------------------------*/


    ipr = IPR_DEFAULT;
    inorm = INORM_DEFAULT;
    itmax = ITMAX_DEFAULT;
    dt = DT_DEFAULT;
    omega = OMEGA_DEFAULT;
    tolrsd[0] = TOLRSD1_DEF;
    tolrsd[1] = TOLRSD2_DEF;
    tolrsd[2] = TOLRSD3_DEF;
    tolrsd[3] = TOLRSD4_DEF;
    tolrsd[4] = TOLRSD5_DEF;
    nx0 = ISIZ1;
    ny0 = ISIZ2;
    nz0 = ISIZ3;

/*--------------------------------------------------------------------
c   check problem size
--------------------------------------------------------------------*/

  if ( nx0 < 4 || ny0 < 4 || nz0 < 4 ) {
#ifdef PFDEBUG
    printf("     PROBLEM SIZE IS TOO SMALL - \n"
	   "     SET EACH OF NX, NY AND NZ AT LEAST EQUAL TO 5\n");
#endif
    return -1;
  }

  if ( nx0 > ISIZ1 || ny0 > ISIZ2 || nz0 > ISIZ3 ) {
#ifdef PFDEBUG
    printf("     PROBLEM SIZE IS TOO LARGE - \n"
	   "     NX, NY AND NZ SHOULD BE EQUAL TO \n"
	   "     ISIZ1, ISIZ2 AND ISIZ3 RESPECTIVELY\n");
#endif
    return -1;
  }

#ifdef PFDEBUG
  printf(" Size: %3dx%3dx%3d\n", nx0, ny0, nz0);
  printf(" Iterations: %3d\n", itmax);
#endif

  return 0;
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void rhs(void) {

#pragma omp parallel 
{

/*--------------------------------------------------------------------
c   compute the right hand sides
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int L1, L2;
  int ist1, iend1;
  int jst1, jend1;
  element_t  q;
  element_t  u21, u31, u41;
  element_t  tmp;
  element_t  u21i, u31i, u41i, u51i;
  element_t  u21j, u31j, u41j, u51j;
  element_t  u21k, u31k, u41k, u51k;
  element_t  u21im1, u31im1, u41im1, u51im1;
  element_t  u21jm1, u31jm1, u41jm1, u51jm1;
  element_t  u21km1, u31km1, u41km1, u51km1;

#pragma omp for  
  for (i = 0; i <= nx-1; i++) {
    for (j = 0; j <= ny-1; j++) {
      for (k = 0; k <= nz-1; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = zero - frct[i][j][k][m];
	}
      }
    }
  }

/*--------------------------------------------------------------------
c   xi-direction flux differences
--------------------------------------------------------------------*/

  L1 = 0;
  L2 = nx-1;

#pragma omp for  
  for (i = L1; i <= L2; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 1; k <= nz - 2; k++) {
	flux[i][j][k][0] = u[i][j][k][1];
	u21 = u[i][j][k][1] / u[i][j][k][0];

	q = (one/two) * (  u[i][j][k][1] * u[i][j][k][1]
		      + u[i][j][k][2] * u[i][j][k][2]
		      + u[i][j][k][3] * u[i][j][k][3] )
	  / u[i][j][k][0];

	flux[i][j][k][1] = u[i][j][k][1] * u21 + C2 * 
	  ( u[i][j][k][4] - q );
	flux[i][j][k][2] = u[i][j][k][2] * u21;
	flux[i][j][k][3] = u[i][j][k][3] * u21;
	flux[i][j][k][4] = ( C1 * u[i][j][k][4] - C2 * q ) * u21;
      }
    } 
  } 

#pragma omp for  
  for (j = jst; j <= jend; j++) {
    for (k = 1; k <= nz - 2; k++) {
      for (i = ist; i <= iend; i++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  rsd[i][j][k][m]
	    - tx2 * ( flux[i+1][j][k][m] - flux[i-1][j][k][m] );
	}
      }

      L2 = nx-1;

      for (i = ist; i <= L2; i++) {
	tmp = one / u[i][j][k][0];

	u21i = tmp * u[i][j][k][1];
	u31i = tmp * u[i][j][k][2];
	u41i = tmp * u[i][j][k][3];
	u51i = tmp * u[i][j][k][4];

	tmp = one / u[i-1][j][k][0];

	u21im1 = tmp * u[i-1][j][k][1];
	u31im1 = tmp * u[i-1][j][k][2];
	u41im1 = tmp * u[i-1][j][k][3];
	u51im1 = tmp * u[i-1][j][k][4];

	flux[i][j][k][1] = (four/three) * tx3 * (u21i-u21im1);
	flux[i][j][k][2] = tx3 * ( u31i - u31im1 );
	flux[i][j][k][3] = tx3 * ( u41i - u41im1 );
	flux[i][j][k][4] = (one/two) * ( one - C1*C5 )
	  * tx3 * (   ( pow2(u21i)   + pow2(u31i)   + pow2(u41i) )
		    - ( pow2(u21im1) + pow2(u31im1) + pow2(u41im1) ) )
	  + (one/six)
	  * tx3 * ( pow2(u21i) - pow2(u21im1) )
	  + C1 * C5 * tx3 * ( u51i - u51im1 );
      }

      for (i = ist; i <= iend; i++) {
	rsd[i][j][k][0] = rsd[i][j][k][0]
	  + dx1 * tx1 * (            u[i-1][j][k][0]
				     - two * u[i][j][k][0]
				     +           u[i+1][j][k][0] );
	rsd[i][j][k][1] = rsd[i][j][k][1]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][1] - flux[i][j][k][1] )
	  + dx2 * tx1 * (            u[i-1][j][k][1]
				     - two * u[i][j][k][1]
				     +           u[i+1][j][k][1] );
	rsd[i][j][k][2] = rsd[i][j][k][2]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][2] - flux[i][j][k][2] )
	  + dx3 * tx1 * (            u[i-1][j][k][2]
				     - two * u[i][j][k][2]
				     +           u[i+1][j][k][2] );
	rsd[i][j][k][3] = rsd[i][j][k][3]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][3] - flux[i][j][k][3] )
	  + dx4 * tx1 * (            u[i-1][j][k][3]
				     - two * u[i][j][k][3]
				     +           u[i+1][j][k][3] );
	rsd[i][j][k][4] = rsd[i][j][k][4]
	  + tx3 * C3 * C4 * ( flux[i+1][j][k][4] - flux[i][j][k][4] )
	  + dx5 * tx1 * (            u[i-1][j][k][4]
				     - two * u[i][j][k][4]
				     +           u[i+1][j][k][4] );
      }

/*--------------------------------------------------------------------
c   Fourth-order dissipation
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	rsd[1][j][k][m] = rsd[1][j][k][m]
	  - dssp * ( + five * u[1][j][k][m]
		     - four * u[2][j][k][m]
		     +           u[3][j][k][m] );
	rsd[2][j][k][m] = rsd[2][j][k][m]
	  - dssp * ( - four * u[1][j][k][m]
		     + six * u[2][j][k][m]
		     - four * u[3][j][k][m]
		     +           u[4][j][k][m] );
      }

      ist1 = 3;
      iend1 = nx - 4;

      for (i = ist1; i <= iend1; i++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = rsd[i][j][k][m]
	    - dssp * (            u[i-2][j][k][m]
				  - four * u[i-1][j][k][m]
				  + six * u[i][j][k][m]
				  - four * u[i+1][j][k][m]
				  +           u[i+2][j][k][m] );
	}
      }


      for (m = 0; m < 5; m++) {
	rsd[nx-3][j][k][m] = rsd[nx-3][j][k][m]
	  - dssp * (             u[nx-5][j][k][m]
				 - four * u[nx-4][j][k][m]
				 + six * u[nx-3][j][k][m]
				 - four * u[nx-2][j][k][m]  );
	rsd[nx-2][j][k][m] = rsd[nx-2][j][k][m]
	  - dssp * (             u[nx-4][j][k][m]
				 - four * u[nx-3][j][k][m]
				 + five * u[nx-2][j][k][m] );
      }
    }
  }

/*--------------------------------------------------------------------
c   eta-direction flux differences
--------------------------------------------------------------------*/

  L1 = 0;
  L2 = ny-1;

#pragma omp for  
  for (i = ist; i <= iend; i++) {
    for (j = L1; j <= L2; j++) {
      for (k = 1; k <= nz - 2; k++) {
	flux[i][j][k][0] = u[i][j][k][2];
	u31 = u[i][j][k][2] / u[i][j][k][0];

	q = (one/two) * (  u[i][j][k][1] * u[i][j][k][1]
		      + u[i][j][k][2] * u[i][j][k][2]
		      + u[i][j][k][3] * u[i][j][k][3] )
	  / u[i][j][k][0];

	flux[i][j][k][1] = u[i][j][k][1] * u31;
	flux[i][j][k][2] = u[i][j][k][2] * u31 + C2 * (u[i][j][k][4]-q);
	flux[i][j][k][3] = u[i][j][k][3] * u31;
	flux[i][j][k][4] = ( C1 * u[i][j][k][4] - C2 * q ) * u31;
      }
    }
  }

#pragma omp for  
  for (i = ist; i <= iend; i++) {
    for (k = 1; k <= nz - 2; k++) {
      for (j = jst; j <= jend; j++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  rsd[i][j][k][m]
	    - ty2 * ( flux[i][j+1][k][m] - flux[i][j-1][k][m] );
	}
      }

      L2 = ny-1;
      for (j = jst; j <= L2; j++) {
	tmp = one / u[i][j][k][0];

	u21j = tmp * u[i][j][k][1];
	u31j = tmp * u[i][j][k][2];
	u41j = tmp * u[i][j][k][3];
	u51j = tmp * u[i][j][k][4];

	tmp = one / u[i][j-1][k][0];
	u21jm1 = tmp * u[i][j-1][k][1];
	u31jm1 = tmp * u[i][j-1][k][2];
	u41jm1 = tmp * u[i][j-1][k][3];
	u51jm1 = tmp * u[i][j-1][k][4];

	flux[i][j][k][1] = ty3 * ( u21j - u21jm1 );
	flux[i][j][k][2] = (four/three) * ty3 * (u31j-u31jm1);
	flux[i][j][k][3] = ty3 * ( u41j - u41jm1 );
	flux[i][j][k][4] = (one/two) * ( one - C1*C5 )
	  * ty3 * (   ( pow2(u21j)   + pow2(u31j)   + pow2(u41j) )
		    - ( pow2(u21jm1) + pow2(u31jm1) + pow2(u41jm1) ) )
	  + (one/six)
	  * ty3 * ( pow2(u31j) - pow2(u31jm1) )
	  + C1 * C5 * ty3 * ( u51j - u51jm1 );
      }

      for (j = jst; j <= jend; j++) {

	rsd[i][j][k][0] = rsd[i][j][k][0]
	  + dy1 * ty1 * (            u[i][j-1][k][0]
				     - two * u[i][j][k][0]
				     +           u[i][j+1][k][0] );

	rsd[i][j][k][1] = rsd[i][j][k][1]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][1] - flux[i][j][k][1] )
	  + dy2 * ty1 * (            u[i][j-1][k][1]
				     - two * u[i][j][k][1]
				     +           u[i][j+1][k][1] );

	rsd[i][j][k][2] = rsd[i][j][k][2]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][2] - flux[i][j][k][2] )
	  + dy3 * ty1 * (            u[i][j-1][k][2]
				     - two * u[i][j][k][2]
				     +           u[i][j+1][k][2] );

	rsd[i][j][k][3] = rsd[i][j][k][3]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][3] - flux[i][j][k][3] )
	  + dy4 * ty1 * (            u[i][j-1][k][3]
				     - two * u[i][j][k][3]
				     +           u[i][j+1][k][3] );

	rsd[i][j][k][4] = rsd[i][j][k][4]
	  + ty3 * C3 * C4 * ( flux[i][j+1][k][4] - flux[i][j][k][4] )
	  + dy5 * ty1 * (            u[i][j-1][k][4]
				     - two * u[i][j][k][4]
				     +           u[i][j+1][k][4] );

      }

/*--------------------------------------------------------------------
c   fourth-order dissipation
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	rsd[i][1][k][m] = rsd[i][1][k][m]
	  - dssp * ( + five * u[i][1][k][m]
		     - four * u[i][2][k][m]
		     +           u[i][3][k][m] );
	rsd[i][2][k][m] = rsd[i][2][k][m]
	  - dssp * ( - four * u[i][1][k][m]
		     + six * u[i][2][k][m]
		     - four * u[i][3][k][m]
		     +           u[i][4][k][m] );
      }

      jst1 = 3;
      jend1 = ny - 4;
      for (j = jst1; j <= jend1; j++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = rsd[i][j][k][m]
	    - dssp * (            u[i][j-2][k][m]
				  - four * u[i][j-1][k][m]
				  + six * u[i][j][k][m]
				  - four * u[i][j+1][k][m]
				  +           u[i][j+2][k][m] );
	}
      }

      for (m = 0; m < 5; m++) {
	rsd[i][ny-3][k][m] = rsd[i][ny-3][k][m]
	  - dssp * (             u[i][ny-5][k][m]
				 - four * u[i][ny-4][k][m]
				 + six * u[i][ny-3][k][m]
				 - four * u[i][ny-2][k][m]  );
	rsd[i][ny-2][k][m] = rsd[i][ny-2][k][m]
	  - dssp * (             u[i][ny-4][k][m]
				 - four * u[i][ny-3][k][m]
				 + five * u[i][ny-2][k][m] );
      }
    }
  }

/*--------------------------------------------------------------------
c   zeta-direction flux differences
--------------------------------------------------------------------*/
#pragma omp for  
  for (i = ist; i <= iend; i++) {
    for (j = jst; j <= jend; j++) {
      for (k = 0; k <= nz-1; k++) {
	flux[i][j][k][0] = u[i][j][k][3];
	u41 = u[i][j][k][3] / u[i][j][k][0];

	q = (one/two) * (  u[i][j][k][1] * u[i][j][k][1]
		      + u[i][j][k][2] * u[i][j][k][2]
		      + u[i][j][k][3] * u[i][j][k][3] )
	  / u[i][j][k][0];

	flux[i][j][k][1] = u[i][j][k][1] * u41;
	flux[i][j][k][2] = u[i][j][k][2] * u41; 
	flux[i][j][k][3] = u[i][j][k][3] * u41 + C2 * (u[i][j][k][4]-q);
	flux[i][j][k][4] = ( C1 * u[i][j][k][4] - C2 * q ) * u41;
      }

      for (k = 1; k <= nz - 2; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] =  rsd[i][j][k][m]
	    - tz2 * ( flux[i][j][k+1][m] - flux[i][j][k-1][m] );
	}
      }

      for (k = 1; k <= nz-1; k++) {
	tmp = one / u[i][j][k][0];

	u21k = tmp * u[i][j][k][1];
	u31k = tmp * u[i][j][k][2];
	u41k = tmp * u[i][j][k][3];
	u51k = tmp * u[i][j][k][4];

	tmp = one / u[i][j][k-1][0];

	u21km1 = tmp * u[i][j][k-1][1];
	u31km1 = tmp * u[i][j][k-1][2];
	u41km1 = tmp * u[i][j][k-1][3];
	u51km1 = tmp * u[i][j][k-1][4];

	flux[i][j][k][1] = tz3 * ( u21k - u21km1 );
	flux[i][j][k][2] = tz3 * ( u31k - u31km1 );
	flux[i][j][k][3] = (four/three) * tz3 * (u41k-u41km1);
	flux[i][j][k][4] = (one/two) * ( one - C1*C5 )
	  * tz3 * (   ( pow2(u21k)   + pow2(u31k)   + pow2(u41k) )
		    - ( pow2(u21km1) + pow2(u31km1) + pow2(u41km1) ) )
	  + (one/six)
	  * tz3 * ( pow2(u41k) - pow2(u41km1) )
	  + C1 * C5 * tz3 * ( u51k - u51km1 );
      }

      for (k = 1; k <= nz - 2; k++) {
	rsd[i][j][k][0] = rsd[i][j][k][0]
	  + dz1 * tz1 * (            u[i][j][k-1][0]
				     - two * u[i][j][k][0]
				     +           u[i][j][k+1][0] );
	rsd[i][j][k][1] = rsd[i][j][k][1]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][1] - flux[i][j][k][1] )
	  + dz2 * tz1 * (            u[i][j][k-1][1]
				     - two * u[i][j][k][1]
				     +           u[i][j][k+1][1] );
	rsd[i][j][k][2] = rsd[i][j][k][2]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][2] - flux[i][j][k][2] )
	  + dz3 * tz1 * (            u[i][j][k-1][2]
				     - two * u[i][j][k][2]
				     +           u[i][j][k+1][2] );
	rsd[i][j][k][3] = rsd[i][j][k][3]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][3] - flux[i][j][k][3] )
	  + dz4 * tz1 * (            u[i][j][k-1][3]
				     - two * u[i][j][k][3]
				     +           u[i][j][k+1][3] );
	rsd[i][j][k][4] = rsd[i][j][k][4]
	  + tz3 * C3 * C4 * ( flux[i][j][k+1][4] - flux[i][j][k][4] )
	  + dz5 * tz1 * (            u[i][j][k-1][4]
				     - two * u[i][j][k][4]
				     +           u[i][j][k+1][4] );
      }

/*--------------------------------------------------------------------
c   fourth-order dissipation
--------------------------------------------------------------------*/
      for (m = 0; m < 5; m++) {
	rsd[i][j][1][m] = rsd[i][j][1][m]
	  - dssp * ( + five * u[i][j][1][m]
		     - four * u[i][j][2][m]
		     +           u[i][j][3][m] );
	rsd[i][j][2][m] = rsd[i][j][2][m]
	  - dssp * ( - four * u[i][j][1][m]
		     + six * u[i][j][2][m]
		     - four * u[i][j][3][m]
		     +           u[i][j][4][m] );
      }

      for (k = 3; k <= nz - 4; k++) {
	for (m = 0; m < 5; m++) {
	  rsd[i][j][k][m] = rsd[i][j][k][m]
	    - dssp * (            u[i][j][k-2][m]
				  - four * u[i][j][k-1][m]
				  + six * u[i][j][k][m]
				  - four * u[i][j][k+1][m]
				  +           u[i][j][k+2][m] );
	}
      }

      for (m = 0; m < 5; m++) {
	rsd[i][j][nz-3][m] = rsd[i][j][nz-3][m]
	  - dssp * (             u[i][j][nz-5][m]
				 - four * u[i][j][nz-4][m]
				 + six * u[i][j][nz-3][m]
				 - four * u[i][j][nz-2][m]  );
	rsd[i][j][nz-2][m] = rsd[i][j][nz-2][m]
	  - dssp * (             u[i][j][nz-4][m]
				 - four * u[i][j][nz-3][m]
				 + five * u[i][j][nz-2][m] );
      }
    }
  }
}

}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void setbv(void) {

#pragma omp parallel
{

/*--------------------------------------------------------------------
c   set the boundary values of dependent variables
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c   local variables
--------------------------------------------------------------------*/
  int i, j, k;
  int iglob, jglob;

/*--------------------------------------------------------------------
c   set the dependent variable values along the top and bottom faces
--------------------------------------------------------------------*/
#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (j = 0; j < ny; j++) {
      jglob = j;
      exact( iglob, jglob, 0, &u[i][j][0][0] );
      exact( iglob, jglob, nz-1, &u[i][j][nz-1][0] );
    }
  }

/*--------------------------------------------------------------------
c   set the dependent variable values along north and south faces
--------------------------------------------------------------------*/
#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, 0, k, &u[i][0][k][0] );
    }
  }

#pragma omp for
  for (i = 0; i < nx; i++) {
    iglob = i;
    for (k = 0; k < nz; k++) {
      exact( iglob, ny0-1,  k, &u[i][ny-1][k][0] );
    }
  }

/*--------------------------------------------------------------------
c   set the dependent variable values along east and west faces
--------------------------------------------------------------------*/
#pragma omp for
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( 0, jglob, k, &u[0][j][k][0] );
    }
  }

#pragma omp for
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 0; k < nz; k++) {
      exact( nx0-1, jglob, k, &u[nx-1][j][k][0] );
    }
  }
}

}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void setcoeff(void) {

/*--------------------------------------------------------------------
c   set up coefficients
--------------------------------------------------------------------*/
  dxi = one / ( nx0 - 1 );
  deta = one / ( ny0 - 1 );
  dzeta = one / ( nz0 - 1 );

  tx1 = one / ( dxi * dxi );
  tx2 = one / ( two * dxi );
  tx3 = one / dxi;

  ty1 = one / ( deta * deta );
  ty2 = one / ( two * deta );
  ty3 = one / deta;

  tz1 = one / ( dzeta * dzeta );
  tz2 = one / ( two * dzeta );
  tz3 = one / dzeta;

  ii1 = 1;
  ii2 = nx0 - 2;
  ji1 = 1;
  ji2 = ny0 - 3;
  ki1 = 2;
  ki2 = nz0 - 2;

/*--------------------------------------------------------------------
c   diffusion coefficients
--------------------------------------------------------------------*/
  dx1 = three/four;
  dx2 = dx1;
  dx3 = dx1;
  dx4 = dx1;
  dx5 = dx1;

  dy1 = three/four;
  dy2 = dy1;
  dy3 = dy1;
  dy4 = dy1;
  dy5 = dy1;

  dz1 = one;
  dz2 = dz1;
  dz3 = dz1;
  dz4 = dz1;
  dz5 = dz1;

/*--------------------------------------------------------------------
c   fourth difference dissipation
--------------------------------------------------------------------*/
  dssp = ( max (dx1, max(dy1, dz1) ) ) / four;

/*--------------------------------------------------------------------
c   coefficients of the exact solution to the first pde
--------------------------------------------------------------------*/
  ce[0][0] = two;
  ce[0][1] = zero;
  ce[0][2] = zero;
  ce[0][3] = four;
  ce[0][4] = five;
  ce[0][5] = three;
  ce[0][6] = five/ten;
  ce[0][7] = two/hundred;
  ce[0][8] = one/hundred;
  ce[0][9] = three/hundred;
  ce[0][10] = five/ten;
  ce[0][11] = four/ten;
  ce[0][12] = three/ten;

/*--------------------------------------------------------------------
c   coefficients of the exact solution to the second pde
--------------------------------------------------------------------*/
  ce[1][0] = one;
  ce[1][1] = zero;
  ce[1][2] = zero;
  ce[1][3] = zero;
  ce[1][4] = one;
  ce[1][5] = two;
  ce[1][6] = three;
  ce[1][7] = one/hundred;
  ce[1][8] = three/hundred;
  ce[1][9] = two/hundred;
  ce[1][10] = four/ten;
  ce[1][11] = three/ten;
  ce[1][12] = five/ten;

/*--------------------------------------------------------------------
c   coefficients of the exact solution to the third pde
--------------------------------------------------------------------*/
  ce[2][0] = two;
  ce[2][1] = two;
  ce[2][2] = zero;
  ce[2][3] = zero;
  ce[2][4] = zero;
  ce[2][5] = two;
  ce[2][6] = three;
  ce[2][7] = four/hundred;
  ce[2][8] = three/hundred;
  ce[2][9] = five/hundred;
  ce[2][10] = three/ten;
  ce[2][11] = five/ten;
  ce[2][12] = four/ten;

/*--------------------------------------------------------------------
c   coefficients of the exact solution to the fourth pde
--------------------------------------------------------------------*/
  ce[3][0] = two;
  ce[3][1] = two;
  ce[3][2] = zero;
  ce[3][3] = zero;
  ce[3][4] = zero;
  ce[3][5] = two;
  ce[3][6] = three;
  ce[3][7] = three/hundred;
  ce[3][8] = five/hundred;
  ce[3][9] = four/hundred;
  ce[3][10] = two/ten;
  ce[3][11] = one/ten;
  ce[3][12] = three/ten;

/*--------------------------------------------------------------------
c   coefficients of the exact solution to the fifth pde
--------------------------------------------------------------------*/
  ce[4][0] = five;
  ce[4][1] = four;
  ce[4][2] = three;
  ce[4][3] = two;
  ce[4][4] = one/ten;
  ce[4][5] = four/ten;
  ce[4][6] = three/ten;
  ce[4][7] = five/hundred;
  ce[4][8] = four/hundred;
  ce[4][9] = three/hundred;
  ce[4][10] = one/ten;
  ce[4][11] = three/ten;
  ce[4][12] = two/ten;
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void setiv(void) {

#pragma omp parallel
{
/*--------------------------------------------------------------------
c
c   set the initial values of independent variables based on tri-linear
c   interpolation of boundary values in the computational space.
c
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int iglob, jglob;
  element_t  xi, eta, zeta;
  element_t  pxi, peta, pzeta;
  element_t  ue_1jk[5],ue_nx0jk[5],ue_i1k[5],
    ue_iny0k[5],ue_ij1[5],ue_ijnz[5];

#pragma omp for
  for (j = 0; j < ny; j++) {
    jglob = j;
    for (k = 1; k < nz - 1; k++) {
      zeta = ((element_t)k) / (nz-1);
      if (jglob != 0 && jglob != ny0-1) {
	eta = ( (element_t) (jglob) ) / (ny0-1);
	for (i = 0; i < nx; i++) {
	  iglob = i;
	  if(iglob != 0 && iglob != nx0-1) {
	    xi = ( (element_t) (iglob) ) / (nx0-1);
	    exact (0,jglob,k,ue_1jk);
	    exact (nx0-1,jglob,k,ue_nx0jk);
	    exact (iglob,0,k,ue_i1k);
	    exact (iglob,ny0-1,k,ue_iny0k);
	    exact (iglob,jglob,0,ue_ij1);
	    exact (iglob,jglob,nz-1,ue_ijnz);
	    for (m = 0; m < 5; m++) {
	      pxi =   ( one - xi ) * ue_1jk[m]
		+ xi   * ue_nx0jk[m];
	      peta =  ( one - eta ) * ue_i1k[m]
		+ eta   * ue_iny0k[m];
	      pzeta = ( one - zeta ) * ue_ij1[m]
		+ zeta   * ue_ijnz[m];

	      u[i][j][k][m] = pxi + peta + pzeta
		- pxi * peta - peta * pzeta - pzeta * pxi
		+ pxi * peta * pzeta;
	    }
	  }
	}
      }
    }
  }
}
}
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void ssor(void) {

/*--------------------------------------------------------------------
c   to perform pseudo-time stepping SSOR iterations
c   for five nonlinear pde s.
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c  local variables
--------------------------------------------------------------------*/
  int i, j, k, m;
  int istep;
  element_t  tmp;
  element_t  delunm[5], tv[ISIZ1][ISIZ2][5];

/*--------------------------------------------------------------------
c   begin pseudo-time stepping iterations
--------------------------------------------------------------------*/
  tmp = one / ( omega * ( two - omega ) ) ;

/*--------------------------------------------------------------------
c   initialize a,b,c,d to zero (guarantees that page tables have been
c   formed, if applicable on given architecture, before timestepping).
--------------------------------------------------------------------*/
#pragma omp parallel private(i,j,k,m)
{
#pragma omp for    
  for (i = 0; i < ISIZ1; i++) {
    for (j = 0; j < ISIZ2; j++) {
      for (k = 0; k < 5; k++) {
	for (m = 0; m < 5; m++) {
	  a[i][j][k][m] = zero;
	  b[i][j][k][m] = zero;
	  c[i][j][k][m] = zero;
	  d[i][j][k][m] = zero;
	}
      }
    }
  }
}
/*--------------------------------------------------------------------
c   compute the steady-state residuals
--------------------------------------------------------------------*/
  rhs();
 
/*--------------------------------------------------------------------
c   compute the L2 norms of newton iteration residuals
--------------------------------------------------------------------*/
  l2norm( nx0, ny0, nz0,
	  ist, iend, jst, jend,
	  rsd, rsdnm );
 
/*--------------------------------------------------------------------
c   the timestep loop
--------------------------------------------------------------------*/


  for (istep = 1; istep <= itmax; istep++) {

    if (istep%20  ==  0 || istep  ==  itmax || istep  ==  1) {
#pragma omp master
#ifdef PFDEBUG
      printf(" Time step %4d\n", istep);
#endif
    }

#pragma omp parallel private(istep,i,j,k,m)
{  
 
/*--------------------------------------------------------------------
c   perform SSOR iteration
--------------------------------------------------------------------*/
#pragma omp for    
    for (i = ist; i <= iend; i++) {
      for (j = jst; j <= jend; j++) {
	for (k = 1; k <= nz - 2; k++) {
	  for (m = 0; m < 5; m++) {
	    rsd[i][j][k][m] = dt * rsd[i][j][k][m];
	  }
	}
      }
    }

    for (k = 1; k <= nz - 2; k++) {
/*--------------------------------------------------------------------
c   form the lower triangular part of the jacobian matrix
--------------------------------------------------------------------*/
      jacld(k);
 
/*--------------------------------------------------------------------
c   perform the lower triangular solution
--------------------------------------------------------------------*/
      blts(nx, ny, nz, k,
	   omega,
	   rsd,
	   a, b, c, d,
	   ist, iend, jst, jend, 
	   nx0, ny0 );
    }
    
#pragma omp barrier

    for (k = nz - 2; k >= 1; k--) {
/*--------------------------------------------------------------------
c   form the strictly upper triangular part of the jacobian matrix
--------------------------------------------------------------------*/
      jacu(k);

/*--------------------------------------------------------------------
c   perform the upper triangular solution
--------------------------------------------------------------------*/
      buts(nx, ny, nz, k,
	   omega,
	   rsd, tv,
	   d, a, b, c,
	   ist, iend, jst, jend,
	   nx0, ny0 );
    }
#pragma omp barrier 
 
/*--------------------------------------------------------------------
c   update the variables
--------------------------------------------------------------------*/

#pragma omp for
    for (i = ist; i <= iend; i++) {
      for (j = jst; j <= jend; j++) {
	for (k = 1; k <= nz-2; k++) {
	  for (m = 0; m < 5; m++) {
	    u[i][j][k][m] = u[i][j][k][m]
	      + tmp * rsd[i][j][k][m];
	  }
	}
      }
    }
} /* end parallel */
/*--------------------------------------------------------------------
c   compute the max-norms of newton iteration corrections
--------------------------------------------------------------------*/
    if ( istep % inorm  ==  0 ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, delunm );
    }
 
/*--------------------------------------------------------------------
c   compute the steady-state residuals
--------------------------------------------------------------------*/
    rhs();
 
/*--------------------------------------------------------------------
c   compute the max-norms of newton iteration residuals
--------------------------------------------------------------------*/
    if ( ( istep % inorm  ==  0 ) ||
	 ( istep  ==  itmax ) ) {
      l2norm( nx0, ny0, nz0,
	      ist, iend, jst, jend,
	      rsd, rsdnm );
    }

/*--------------------------------------------------------------------
c   check the newton-iteration residuals against the tolerance levels
--------------------------------------------------------------------*/
    if ( ( rsdnm[0] < tolrsd[0] ) &&
	 ( rsdnm[1] < tolrsd[1] ) &&
	 ( rsdnm[2] < tolrsd[2] ) &&
	 ( rsdnm[3] < tolrsd[3] ) &&
	 ( rsdnm[4] < tolrsd[4] ) ) {
    }
  }

}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void verify(element_t xcr[5], element_t xce[5], element_t xci,
		   char *class, boolean *verified) {

/*--------------------------------------------------------------------
c  verification routine                         
--------------------------------------------------------------------*/

  element_t xcrref[5],xceref[5],xciref,
    xcrdif[5],xcedif[5],xcidif,
    epsilon, dtref;
  int m;

/*--------------------------------------------------------------------
c   tolerance level
--------------------------------------------------------------------*/
  epsilon = c_epsilon;

  *class = 'U';
  *verified = TRUE;

  for (m = 0; m < 5; m++) {
    xcrref[m] = one;
    xceref[m] = one;
  }
  xciref = one;

  if ( nx0 == 12 && ny0 == 12 && nz0 == 12 && itmax == 50)  {
    *class = 'S';
    dtref = five/ten;

/*--------------------------------------------------------------------
c   Reference values of RMS-norms of residual, for the (12X12X12) grid,
c   after 50 time steps, with  DT = 5.0d-01
--------------------------------------------------------------------*/
#ifdef WITH_POSIT_32
    *((uint32_t*)&xcrref[0]) = 0x2825718c;
    *((uint32_t*)&xcrref[1]) = 0x1e401b70;
    *((uint32_t*)&xcrref[2]) = 0x1d1bdd8a;
    *((uint32_t*)&xcrref[3]) = 0x1d13fbaa;
    *((uint32_t*)&xcrref[4]) = 0x2c62c3e0;
#elif (defined WITH_POSIT_16)
    *((uint32_t*)&xcrref[0]) = 0x1825;
    *((uint32_t*)&xcrref[1]) = 0xe40;
    *((uint32_t*)&xcrref[2]) = 0xd1b;
    *((uint32_t*)&xcrref[3]) = 0xd13;
    *((uint32_t*)&xcrref[4]) = 0x1c62;
#elif (defined WITH_POSIT_8)
    *((uint32_t*)&xcrref[0]) = 0x8;
    *((uint32_t*)&xcrref[1]) = 0x3;
    *((uint32_t*)&xcrref[2]) = 0x2;
    *((uint32_t*)&xcrref[3]) = 0x2;
    *((uint32_t*)&xcrref[4]) = 0xc;
#else
    xcrref[0] = 1.6196343210976702e-02;
    xcrref[1] = 2.1976745164821318e-03;
    xcrref[2] = 1.5179927653399185e-03;
    xcrref[3] = 1.5029584435994323e-03;
    xcrref[4] = 3.4264073155896461e-02;
#endif
/*--------------------------------------------------------------------
c   Reference values of RMS-norms of solution error, for the (12X12X12) grid,
c   after 50 time steps, with  DT = 5.0d-01
--------------------------------------------------------------------*/
#ifdef WITH_POSIT_32
    *((uint32_t*)&xceref[0]) = 0x1aa16e29;
    *((uint32_t*)&xceref[1]) = 0x14c1da99;
    *((uint32_t*)&xceref[2]) = 0x13d6f2b3;
    *((uint32_t*)&xceref[3]) = 0x13d508df;
    *((uint32_t*)&xceref[4]) = 0x1caefe28;
#elif (defined WITH_POSIT_16)
    *((uint32_t*)&xcrref[0]) = 0xaa1;
    *((uint32_t*)&xcrref[1]) = 0x660;
    *((uint32_t*)&xcrref[2]) = 0x5eb;
    *((uint32_t*)&xcrref[3]) = 0x5ea;
    *((uint32_t*)&xcrref[4]) = 0xcae;
#elif (defined WITH_POSIT_8)
    *((uint32_t*)&xcrref[0]) = 0x1;
    *((uint32_t*)&xcrref[1]) = 0x1;
    *((uint32_t*)&xcrref[2]) = 0x1;
    *((uint32_t*)&xcrref[3]) = 0x1;
    *((uint32_t*)&xcrref[4]) = 0x2;
#else
    xceref[0] = 6.4223319957960924e-04;
    xceref[1] = 8.4144342047347926e-05;
    xceref[2] = 5.8588269616485186e-05;
    xceref[3] = 5.8474222595157350e-05;
    xceref[4] = 1.3103347914111294e-03;
#endif
/*--------------------------------------------------------------------
c   Reference value of surface integral, for the (12X12X12) grid,
c   after 50 time steps, with DT = 5.0d-01
--------------------------------------------------------------------*/
#ifdef WITH_POSIT_32
    *((uint32_t*)&xciref) = 0x4bd7864a;
#elif (defined WITH_POSIT_16)
    *((uint32_t*)&xciref) = 0x57af;
#elif (defined WITH_POSIT_8)
    *((uint32_t*)&xciref) = 0x67;
#else
    xciref = 7.8418928865937083;
#endif
  }

/*
  } else if ( nx0 == 33 && ny0 == 33 && nz0 == 33 && itmax == 300) {

    *class = 'W';   // SPEC95fp size
    dtref = 1.5e-3;
//--------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (33x33x33) grid,
//   after 300 time steps, with  DT = 1.5d-3
//--------------------------------------------------------------------
    xcrref[0] =   0.1236511638192e+02;
    xcrref[1] =   0.1317228477799e+01;
    xcrref[2] =   0.2550120713095e+01;
    xcrref[3] =   0.2326187750252e+01;
    xcrref[4] =   0.2826799444189e+02;

//--------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (33X33X33) grid,
//--------------------------------------------------------------------
    xceref[0] =   0.4867877144216;
    xceref[1] =   0.5064652880982e-01;
    xceref[2] =   0.9281818101960e-01;
    xceref[3] =   0.8570126542733e-01;
    xceref[4] =   0.1084277417792e+01;

//--------------------------------------------------------------------
//   Reference value of surface integral, for the (33X33X33) grid,
//   after 300 time steps, with  DT = 1.5d-3
//--------------------------------------------------------------------
    xciref    =   0.1161399311023e+02;

  } else if ( nx0 == 64 && ny0 == 64 && nz0 == 64 && itmax == 250) {

    *class = 'A';
    dtref = two;
//--------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (64X64X64) grid,
//   after 250 time steps, with  DT = 2.0d+0.0
//--------------------------------------------------------------------
    xcrref[0] = 7.7902107606689367e+02;
    xcrref[1] = 6.3402765259692870e+01;
    xcrref[2] = 1.9499249727292479e+02;
    xcrref[3] = 1.7845301160418537e+02;
    xcrref[4] = 1.8384760349464247e+03;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of solution error, for the (64X64X64) grid,
//  after 250 time steps, with  DT = 2.0d+0.0
//--------------------------------------------------------------------
    xceref[0] = 2.9964085685471943e+01;
    xceref[1] = 2.8194576365003349;
    xceref[2] = 7.3473412698774742;
    xceref[3] = 6.7139225687777051;
    xceref[4] = 7.0715315688392578e+01;

//--------------------------------------------------------------------
//  Reference value of surface integral, for the (64X64X64) grid,
//   after 250 time steps, with DT = 2.0d+0.0
//--------------------------------------------------------------------
    xciref = 2.6030925604886277e+01;

    } else if ( nx0 == 102 && ny0 == 102 && nz0 == 102 && itmax == 250) {

      *class = 'B';
      dtref = two;

//--------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (102X102X102) grid,
//   after 250 time steps, with  DT = 2.0d+0.0
//--------------------------------------------------------------------
      xcrref[0] = 3.5532672969982736e+03;
      xcrref[1] = 2.6214750795310692e+02;
      xcrref[2] = 8.8333721850952190e+02;
      xcrref[3] = 7.7812774739425265e+02;
      xcrref[4] = 7.3087969592545314e+03;

//--------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (102X102X102)
//   grid, after 250 time steps, with  DT = 2.0d+0.0
//--------------------------------------------------------------------
      xceref[0] = 1.1401176380212709e+02;
      xceref[1] = 8.1098963655421574;
      xceref[2] = 2.8480597317698308e+01;
      xceref[3] = 2.5905394567832939e+01;
      xceref[4] = 2.6054907504857413e+02;

//--------------------------------------------------------------------
//   Reference value of surface integral, for the (102X102X102) grid,
//   after 250 time steps, with DT = 2.0d+0.0
//--------------------------------------------------------------------
      xciref = 4.7887162703308227e+01;

      } else if ( nx0 == 162 && ny0 == 162 && nz0 == 162 && itmax == 250) {

	*class = 'C';
	dtref = two;

//--------------------------------------------------------------------
//   Reference values of RMS-norms of residual, for the (162X162X162) grid,
//   after 250 time steps, with  DT = 2.0d+0.0
//--------------------------------------------------------------------
	xcrref[0] = 1.03766980323537846e+04;
	xcrref[1] = 8.92212458801008552e+02;
	xcrref[2] = 2.56238814582660871e+03;
	xcrref[3] = 2.19194343857831427e+03;
	xcrref[4] = 1.78078057261061185e+04;

//--------------------------------------------------------------------
//   Reference values of RMS-norms of solution error, for the (162X162X162)
//   grid, after 250 time steps, with  DT = 2.0d+0.0
//--------------------------------------------------------------------
	xceref[0] = 2.15986399716949279e+02;
	xceref[1] = 1.55789559239863600e+01;
	xceref[2] = 5.41318863077207766e+01;
	xceref[3] = 4.82262643154045421e+01;
	xceref[4] = 4.55902910043250358e+02;

//--------------------------------------------------------------------
//   Reference value of surface integral, for the (162X162X162) grid,
//   after 250 time steps, with DT = 2.0d+0.0
//--------------------------------------------------------------------
	xciref = 6.66404553572181300e+01;
      } else {
	*verified = FALSE;
      }
*/

/*--------------------------------------------------------------------
c    verification test for residuals if gridsize is either 12X12X12 or 
c    64X64X64 or 102X102X102 or 162X162X162
--------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c    Compute the difference of solution values and the known reference values.
--------------------------------------------------------------------*/
  for (m = 0; m < 5; m++) {
           
    xcrdif[m] = my_fabs((xcr[m]-xcrref[m])/xcrref[m]);
    xcedif[m] = my_fabs((xce[m]-xceref[m])/xceref[m]);
           
  }
  xcidif = my_fabs((xci - xciref)/xciref);

/*--------------------------------------------------------------------
c    Output the comparison of computed results to known cases.
--------------------------------------------------------------------*/

  if (*class != 'U') {
#ifdef PFDEBUG
    printf("\n Verification being performed for class %1c\n", *class);
    printf(" Accuracy setting for epsilon = %20.13e\n", epsilon);
#endif
    if (my_fabs(dt-dtref) > epsilon) {
      *verified = FALSE;
      *class = 'U';
#ifdef PFDEBUG
      printf(" DT does not match the reference value of %15.8e\n", dtref);
#endif
    }
  }
#ifdef PFDEBUG
  else {
    printf(" Unknown class\n");
  }
#endif
#ifdef PFDEBUG
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of residual\n");
  } else {
    printf(" RMS-norms of residual\n");
  }
#endif

  for (m = 0; m < 5; m++) {
    if (*class  ==  'U') {
#ifdef PFDEBUG
      printf("          %2d  %20.13e\n", m, xcr[m]);
#endif
    } else if (xcrdif[m] > epsilon) {
      *verified = FALSE;
#ifdef PFDEBUG
      printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n",
	     m,xcr[m],xcrref[m],xcrdif[m]);
#endif
    } else {
#ifdef PFDEBUG
      printf("          %2d  %20.13e%20.13e%20.13e\n",
	     m,xcr[m],xcrref[m],xcrdif[m]);
#endif
    }
  }
#ifdef PFDEBUG
  if (*class != 'U') {
    printf(" Comparison of RMS-norms of solution error\n");
  } else {
    printf(" RMS-norms of solution error\n");
  }
#endif
  for (m = 0; m < 5; m++) {
    if (*class  ==  'U') {
#ifdef PFDEBUG
      printf("          %2d  %20.13e\n", m, xce[m]);
#endif
    } else if (xcedif[m] > epsilon) {
      *verified = FALSE;
#ifdef PFDEBUG
      printf(" FAILURE: %2d  %20.13e%20.13e%20.13e\n",
	     m,xce[m],xceref[m],xcedif[m]);
#endif
    }
#ifdef PFDEBUG
    else {
      printf("          %2d  %20.13e%20.13e%20.13e\n",
	     m,xce[m],xceref[m],xcedif[m]);
    }
#endif
  }
#ifdef PFDEBUG
  if (*class != 'U') {
    printf(" Comparison of surface integral\n");
  } else {
    printf(" Surface integral\n");
  }
#endif
  if (*class  ==  'U') {
#ifdef PFDEBUG
    printf("              %20.13e\n", xci);
#endif
  } else if (xcidif > epsilon) {
    *verified = FALSE;
#ifdef PFDEBUG
    printf(" FAILURE:     %20.13e%20.13e%20.13e\n", 
	   xci, xciref, xcidif);
#endif
  } else {
#ifdef PFDEBUG
    printf("              %20.13e%20.13e%20.13e\n",
	   xci, xciref, xcidif);
#endif
  }
#ifdef PFDEBUG
  if (*class  ==  'U') {
    printf(" No reference values provided\n");
    printf(" No verification performed\n");
  } else if (*verified) {
    printf(" Verification Successful\n");
  } else {
    printf(" Verification failed\n");
  }
#endif
}
