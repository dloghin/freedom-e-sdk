/*--------------------------------------------------------------------
  
  NAS Parallel Benchmarks 3.0 structured OpenMP C versions - BT

  This benchmark is an OpenMP C version of the NPB BT code.
  
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

  Authors: R. Van der Wijngaart
           T. Harris
           M. Yarrow

  OpenMP C version: S. Satoh

  3.0 structure translation: M. Popov
  
--------------------------------------------------------------------*/

// https://github.com/benchmark-subsetting/NPB3-0-omp-C

/* global variables */
#include "header.h"

#include <stdint.h>

#include "../common/perf.h"

// #define WITH_POSIT_32
// #define PFDEBUG

#ifdef PFDEBUG
#include <stdio.h>
#endif

// constants
element_t zero, one, two, three, four, five, six, ten, hundred, c_epsilon, c_dtref;

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
uint32_t posit_1 = 0x32666666;
uint32_t posit_01 = 0x251eb851;
uint32_t posit_001 = 0x1c0c49ba;
uint32_t posit_0001 = 0x1546dc5d;
uint32_t posit_00001 = 0xf4f8b58;
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
	*((uint32_t*)&c_dtref) = posit_01;
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
	*((uint32_t*)&c_epsilon) = fp32_00001;
	*((uint32_t*)&c_dtref) = fp32_01;
#endif /* WITH_POSIT */
}

#define	DT_DEFAULT	(one/hundred)

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
static void add(void);
static void adi(void);
static void error_norm(element_t rms[5]);
static void rhs_norm(element_t rms[5]);
static void exact_rhs(void);
static void exact_solution(element_t xi, element_t eta, element_t zeta,
		element_t dtemp[5]);
static void initialize(void);
static void lhsinit(void);
static void lhsx(void);
static void lhsy(void);
static void lhsz(void);
static void compute_rhs(void);
static void set_constants(void);
static void verify(int no_time_steps, char *class, boolean *verified);
static void x_solve(void);
static void x_backsubstitute(void);
static void x_solve_cell(void);
static void matvec_sub(element_t ablock[5][5], element_t avec[5], element_t bvec[5]);
static void matmul_sub(element_t ablock[5][5], element_t bblock[5][5],
		element_t cblock[5][5]);
static void binvcrhs(element_t lhs[5][5], element_t c[5][5], element_t r[5]);
static void binvrhs(element_t lhs[5][5], element_t r[5]);
static void y_solve(void);
static void y_backsubstitute(void);
static void y_solve_cell(void);
static void z_solve(void);
static void z_backsubstitute(void);
static void z_solve_cell(void);

/*--------------------------------------------------------------------
      program BT
c-------------------------------------------------------------------*/
int main(int argc, char **argv) {

	int niter, step, n3;
	int nthreads = 1;
	element_t navg, mflops;

	element_t tmax;
	boolean verified;
	char class;

	init_constants();

	/*--------------------------------------------------------------------
c      Root node reads input file (if it exists) else takes
c      defaults from parameters
c-------------------------------------------------------------------*/

	niter = NITER_DEFAULT;
	dt    = DT_DEFAULT;
	grid_points[0] = PROBLEM_SIZE;
	grid_points[1] = PROBLEM_SIZE;
	grid_points[2] = PROBLEM_SIZE;

#ifdef PFDEBUG
	printf("\n\n NAS Parallel Benchmarks 3_0 structured OpenMP C version"
			" - BT Benchmark\n\n");
	printf(" Size: %3dx%3dx%3d\n",
			grid_points[0], grid_points[1], grid_points[2]);
	printf(" Iterations: %3d   dt: %10.6f\n", niter, dt);
#endif

	if (grid_points[0] > IMAX ||
			grid_points[1] > JMAX ||
			grid_points[2] > KMAX) {
#ifdef PFDEBUG
		printf(" %dx%dx%d\n", grid_points[0], grid_points[1], grid_points[2]);
		printf(" Problem size too big for compiled array sizes\n");
#endif
		return -1;
	}

	unsigned long long startc = read_cycles();

	set_constants();

	initialize();
	lhsinit();
	exact_rhs();

	/*--------------------------------------------------------------------
c      do one time step to touch all code, and reinitialize
c-------------------------------------------------------------------*/

	adi();

	initialize();

	//  timer_clear(1);
	//  timer_start(1);

	for (step = 1; step <= niter; step++) {
#ifdef PFDEBUG
		if (step%20 == 0 || step == 1) {
			printf(" Time step %4d\n", step);
		}
#endif
		adi();
	}

	verify(niter, &class, &verified);

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
c-------------------------------------------------------------------*/

static void add(void) {

	/*--------------------------------------------------------------------
c     addition of update to the vector u
c-------------------------------------------------------------------*/

	int i, j, k, m;

#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = u[i][j][k][m] + rhs[i][j][k][m];
				}
			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void adi(void) {
#pragma omp parallel
	compute_rhs();

#pragma omp parallel
	x_solve();

#pragma omp parallel
	y_solve();

#pragma omp parallel
	z_solve();

#pragma omp parallel
	add();
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void error_norm(element_t rms[5]) {

	/*--------------------------------------------------------------------
c     this function computes the norm of the difference between the
c     computed solution and the exact solution
c-------------------------------------------------------------------*/

	int i, j, k, m, d;
	element_t xi, eta, zeta, u_exact[5], add;

	for (m = 0; m < 5; m++) {
		rms[m] = zero;
	}

	for (i = 0; i < grid_points[0]; i++) {
		xi = (element_t)i * dnxm1;
		for (j = 0; j < grid_points[1]; j++) {
			eta = (element_t)j * dnym1;
			for (k = 0; k < grid_points[2]; k++) {
				zeta = (element_t)k * dnzm1;
				exact_solution(xi, eta, zeta, u_exact);

				for (m = 0; m < 5; m++) {
					add = u[i][j][k][m] - u_exact[m];
					rms[m] = rms[m] + add*add;
				}
			}
		}
	}

  for (m = 0; m < 5; m++) {
    for (d = 0; d <= 2; d++) {
      rms[m] = rms[m] / (element_t)(grid_points[d]-2);
    }
    rms[m] = sqrt_asm(rms[m]);
  }
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void rhs_norm(element_t rms[5]) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	int i, j, k, d, m;
	element_t add;

	for (m = 0; m < 5; m++) {
		rms[m] = zero;
	}

	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				for (m = 0; m < 5; m++) {
					add = rhs[i][j][k][m];
					rms[m] = rms[m] + add*add;
				}
			}
		}
	}

	for (m = 0; m < 5; m++) {
		for (d = 0; d <= 2; d++) {
			rms[m] = rms[m] / (element_t)(grid_points[d]-2);
		}
		rms[m] = sqrt_asm(rms[m]);
	}
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void exact_rhs(void) {

#pragma omp parallel
	{
		/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

		/*--------------------------------------------------------------------
c     compute the right hand side based on exact solution
c-------------------------------------------------------------------*/

		element_t dtemp[5], xi, eta, zeta, dtpp;
		int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

		/*--------------------------------------------------------------------
c     initialize                                  
c-------------------------------------------------------------------*/
#pragma omp for  
		for (i = 0; i < grid_points[0]; i++) {
			for (j = 0; j < grid_points[1]; j++) {
				for (k = 0; k < grid_points[2]; k++) {
					for (m = 0; m < 5; m++) {
						forcing[i][j][k][m] = zero;
					}
				}
			}
		}

		/*--------------------------------------------------------------------
c     xi-direction flux differences                      
c-------------------------------------------------------------------*/
#pragma omp for
		for (j = 1; j < grid_points[1]-1; j++) {
			eta = (element_t)j * dnym1;

			for (k = 1; k < grid_points[2]-1; k++) {
				zeta = (element_t)k * dnzm1;

				for (i = 0; i < grid_points[0]; i++) {
					xi = (element_t)i * dnxm1;

					exact_solution(xi, eta, zeta, dtemp);
					for (m = 0; m < 5; m++) {
						ue[i][m] = dtemp[m];
					}

					dtpp = one / dtemp[0];

					for (m = 1; m <= 4; m++) {
						buf[i][m] = dtpp * dtemp[m];
					}

					cuf[i]   = buf[i][1] * buf[i][1];
					buf[i][0] = cuf[i] + buf[i][2] * buf[i][2] +
							buf[i][3] * buf[i][3];
					q[i] = one/two*(buf[i][1]*ue[i][1] + buf[i][2]*ue[i][2] +
							buf[i][3]*ue[i][3]);
				}

				for (i = 1; i < grid_points[0]-1; i++) {
					im1 = i-1;
					ip1 = i+1;

					forcing[i][j][k][0] = forcing[i][j][k][0] -
							tx2*(ue[ip1][1]-ue[im1][1])+
							dx1tx1*(ue[ip1][0]-two*ue[i][0]+ue[im1][0]);

					forcing[i][j][k][1] = forcing[i][j][k][1] -
							tx2 * ((ue[ip1][1]*buf[ip1][1]+c2*(ue[ip1][4]-q[ip1]))-
									(ue[im1][1]*buf[im1][1]+c2*(ue[im1][4]-q[im1])))+
									xxcon1*(buf[ip1][1]-two*buf[i][1]+buf[im1][1])+
									dx2tx1*( ue[ip1][1]-two* ue[i][1]+ ue[im1][1]);

					forcing[i][j][k][2] = forcing[i][j][k][2] -
							tx2 * (ue[ip1][2]*buf[ip1][1]-ue[im1][2]*buf[im1][1])+
							xxcon2*(buf[ip1][2]-two*buf[i][2]+buf[im1][2])+
							dx3tx1*( ue[ip1][2]-two* ue[i][2]+ ue[im1][2]);

					forcing[i][j][k][3] = forcing[i][j][k][3] -
							tx2*(ue[ip1][3]*buf[ip1][1]-ue[im1][3]*buf[im1][1])+
							xxcon2*(buf[ip1][3]-two*buf[i][3]+buf[im1][3])+
							dx4tx1*( ue[ip1][3]-two* ue[i][3]+ ue[im1][3]);

					forcing[i][j][k][4] = forcing[i][j][k][4] -
							tx2*(buf[ip1][1]*(c1*ue[ip1][4]-c2*q[ip1])-
									buf[im1][1]*(c1*ue[im1][4]-c2*q[im1]))+
									one/two*xxcon3*(buf[ip1][0]-two*buf[i][0]+buf[im1][0])+
									xxcon4*(cuf[ip1]-two*cuf[i]+cuf[im1])+
									xxcon5*(buf[ip1][4]-two*buf[i][4]+buf[im1][4])+
									dx5tx1*( ue[ip1][4]-two* ue[i][4]+ ue[im1][4]);
				}

				/*--------------------------------------------------------------------
c     Fourth-order dissipation                         
c-------------------------------------------------------------------*/

				for (m = 0; m < 5; m++) {
					i = 1;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(five*ue[i][m] - four*ue[i+1][m] +ue[i+2][m]);
					i = 2;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(-four*ue[i-1][m] + six*ue[i][m] -
									four*ue[i+1][m] +     ue[i+2][m]);
				}

				for (m = 0; m < 5; m++) {
					for (i = 1*3; i <= grid_points[0]-3*1-1; i++) {
						forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
								(ue[i-2][m] - four*ue[i-1][m] +
										six*ue[i][m] - four*ue[i+1][m] + ue[i+2][m]);
					}
				}

				for (m = 0; m < 5; m++) {
					i = grid_points[0]-3;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(ue[i-2][m] - four*ue[i-1][m] +
									six*ue[i][m] - four*ue[i+1][m]);
					i = grid_points[0]-2;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(ue[i-2][m] - four*ue[i-1][m] + five*ue[i][m]);
				}

			}
		}

		/*--------------------------------------------------------------------
c     eta-direction flux differences             
c-------------------------------------------------------------------*/
#pragma omp for
		for (i = 1; i < grid_points[0]-1; i++) {
			xi = (element_t)i * dnxm1;

			for (k = 1; k < grid_points[2]-1; k++) {
				zeta = (element_t)k * dnzm1;

				for (j = 0; j < grid_points[1]; j++) {
					eta = (element_t)j * dnym1;

					exact_solution(xi, eta, zeta, dtemp);
					for (m = 0; m < 5; m++) {
						ue[j][m] = dtemp[m];
					}

					dtpp = one/dtemp[0];

					for (m = 1; m <= 4; m++) {
						buf[j][m] = dtpp * dtemp[m];
					}

					cuf[j]   = buf[j][2] * buf[j][2];
					buf[j][0] = cuf[j] + buf[j][1] * buf[j][1] +
							buf[j][3] * buf[j][3];
					q[j] = one/two*(buf[j][1]*ue[j][1] + buf[j][2]*ue[j][2] +
							buf[j][3]*ue[j][3]);
				}

				for (j = 1; j < grid_points[1]-1; j++) {
					jm1 = j-1;
					jp1 = j+1;

					forcing[i][j][k][0] = forcing[i][j][k][0] -
							ty2*( ue[jp1][2]-ue[jm1][2] )+
							dy1ty1*(ue[jp1][0]-two*ue[j][0]+ue[jm1][0]);

					forcing[i][j][k][1] = forcing[i][j][k][1] -
							ty2*(ue[jp1][1]*buf[jp1][2]-ue[jm1][1]*buf[jm1][2])+
							yycon2*(buf[jp1][1]-two*buf[j][1]+buf[jm1][1])+
							dy2ty1*( ue[jp1][1]-two* ue[j][1]+ ue[jm1][1]);

					forcing[i][j][k][2] = forcing[i][j][k][2] -
							ty2*((ue[jp1][2]*buf[jp1][2]+c2*(ue[jp1][4]-q[jp1]))-
									(ue[jm1][2]*buf[jm1][2]+c2*(ue[jm1][4]-q[jm1])))+
									yycon1*(buf[jp1][2]-two*buf[j][2]+buf[jm1][2])+
									dy3ty1*( ue[jp1][2]-two*ue[j][2] +ue[jm1][2]);

					forcing[i][j][k][3] = forcing[i][j][k][3] -
							ty2*(ue[jp1][3]*buf[jp1][2]-ue[jm1][3]*buf[jm1][2])+
							yycon2*(buf[jp1][3]-two*buf[j][3]+buf[jm1][3])+
							dy4ty1*( ue[jp1][3]-two*ue[j][3]+ ue[jm1][3]);

					forcing[i][j][k][4] = forcing[i][j][k][4] -
							ty2*(buf[jp1][2]*(c1*ue[jp1][4]-c2*q[jp1])-
									buf[jm1][2]*(c1*ue[jm1][4]-c2*q[jm1]))+
									one/two*yycon3*(buf[jp1][0]-two*buf[j][0]+
											buf[jm1][0])+
											yycon4*(cuf[jp1]-two*cuf[j]+cuf[jm1])+
											yycon5*(buf[jp1][4]-two*buf[j][4]+buf[jm1][4])+
											dy5ty1*(ue[jp1][4]-two*ue[j][4]+ue[jm1][4]);
				}

				/*--------------------------------------------------------------------
c     Fourth-order dissipation                      
c-------------------------------------------------------------------*/
				for (m = 0; m < 5; m++) {
					j = 1;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(five*ue[j][m] - four*ue[j+1][m] +ue[j+2][m]);
					j = 2;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(-four*ue[j-1][m] + six*ue[j][m] -
									four*ue[j+1][m] +       ue[j+2][m]);
				}

				for (m = 0; m < 5; m++) {
					for (j = 1*3; j <= grid_points[1]-3*1-1; j++) {
						forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
								(ue[j-2][m] - four*ue[j-1][m] +
										six*ue[j][m] - four*ue[j+1][m] + ue[j+2][m]);
					}
				}

				for (m = 0; m < 5; m++) {
					j = grid_points[1]-3;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(ue[j-2][m] - four*ue[j-1][m] +
									six*ue[j][m] - four*ue[j+1][m]);
					j = grid_points[1]-2;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(ue[j-2][m] - four*ue[j-1][m] + five*ue[j][m]);
				}

			}
		}


		/*--------------------------------------------------------------------
c     zeta-direction flux differences                      
c-------------------------------------------------------------------*/
#pragma omp for
		for (i = 1; i < grid_points[0]-1; i++) {
			xi = (element_t)i * dnxm1;

			for (j = 1; j < grid_points[1]-1; j++) {
				eta = (element_t)j * dnym1;

				for (k = 0; k < grid_points[2]; k++) {
					zeta = (element_t)k * dnzm1;

					exact_solution(xi, eta, zeta, dtemp);
					for (m = 0; m < 5; m++) {
						ue[k][m] = dtemp[m];
					}

					dtpp = one/dtemp[0];

					for (m = 1; m <= 4; m++) {
						buf[k][m] = dtpp * dtemp[m];
					}

					cuf[k]   = buf[k][3] * buf[k][3];
					buf[k][0] = cuf[k] + buf[k][1] * buf[k][1] +
							buf[k][2] * buf[k][2];
					q[k] = one/two*(buf[k][1]*ue[k][1] + buf[k][2]*ue[k][2] +
							buf[k][3]*ue[k][3]);
				}

				for (k = 1; k < grid_points[2]-1; k++) {
					km1 = k-1;
					kp1 = k+1;

					forcing[i][j][k][0] = forcing[i][j][k][0] -
							tz2*( ue[kp1][3]-ue[km1][3] )+
							dz1tz1*(ue[kp1][0]-two*ue[k][0]+ue[km1][0]);

					forcing[i][j][k][1] = forcing[i][j][k][1] -
							tz2 * (ue[kp1][1]*buf[kp1][3]-ue[km1][1]*buf[km1][3])+
							zzcon2*(buf[kp1][1]-two*buf[k][1]+buf[km1][1])+
							dz2tz1*( ue[kp1][1]-two* ue[k][1]+ ue[km1][1]);

					forcing[i][j][k][2] = forcing[i][j][k][2] -
							tz2 * (ue[kp1][2]*buf[kp1][3]-ue[km1][2]*buf[km1][3])+
							zzcon2*(buf[kp1][2]-two*buf[k][2]+buf[km1][2])+
							dz3tz1*(ue[kp1][2]-two*ue[k][2]+ue[km1][2]);

					forcing[i][j][k][3] = forcing[i][j][k][3] -
							tz2 * ((ue[kp1][3]*buf[kp1][3]+c2*(ue[kp1][4]-q[kp1]))-
									(ue[km1][3]*buf[km1][3]+c2*(ue[km1][4]-q[km1])))+
									zzcon1*(buf[kp1][3]-two*buf[k][3]+buf[km1][3])+
									dz4tz1*( ue[kp1][3]-two*ue[k][3] +ue[km1][3]);

					forcing[i][j][k][4] = forcing[i][j][k][4] -
							tz2 * (buf[kp1][3]*(c1*ue[kp1][4]-c2*q[kp1])-
									buf[km1][3]*(c1*ue[km1][4]-c2*q[km1]))+
									one/two*zzcon3*(buf[kp1][0]-two*buf[k][0]
																		   +buf[km1][0])+
																		   zzcon4*(cuf[kp1]-two*cuf[k]+cuf[km1])+
																		   zzcon5*(buf[kp1][4]-two*buf[k][4]+buf[km1][4])+
																		   dz5tz1*( ue[kp1][4]-two*ue[k][4]+ ue[km1][4]);
				}

				/*--------------------------------------------------------------------
c     Fourth-order dissipation                        
c-------------------------------------------------------------------*/
				for (m = 0; m < 5; m++) {
					k = 1;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(five*ue[k][m] - four*ue[k+1][m] +ue[k+2][m]);
					k = 2;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(-four*ue[k-1][m] + six*ue[k][m] -
									four*ue[k+1][m] +       ue[k+2][m]);
				}

				for (m = 0; m < 5; m++) {
					for (k = 1*3; k <= grid_points[2]-3*1-1; k++) {
						forcing[i][j][k][m] = forcing[i][j][k][m] - dssp*
								(ue[k-2][m] - four*ue[k-1][m] +
										six*ue[k][m] - four*ue[k+1][m] + ue[k+2][m]);
					}
				}

				for (m = 0; m < 5; m++) {
					k = grid_points[2]-3;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(ue[k-2][m] - four*ue[k-1][m] +
									six*ue[k][m] - four*ue[k+1][m]);
					k = grid_points[2]-2;
					forcing[i][j][k][m] = forcing[i][j][k][m] - dssp *
							(ue[k-2][m] - four*ue[k-1][m] + five*ue[k][m]);
				}

			}
		}

		/*--------------------------------------------------------------------
c     now change the sign of the forcing function, 
c-------------------------------------------------------------------*/
#pragma omp for  
		for (i = 1; i < grid_points[0]-1; i++) {
			for (j = 1; j < grid_points[1]-1; j++) {
				for (k = 1; k < grid_points[2]-1; k++) {
					for (m = 0; m < 5; m++) {
						forcing[i][j][k][m] = -one * forcing[i][j][k][m];
					}
				}
			}
		}
	}

}
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void exact_solution(element_t xi, element_t eta, element_t zeta,
		element_t dtemp[5]) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     this function returns the exact solution at point xi, eta, zeta  
c-------------------------------------------------------------------*/

	int m;

	for (m = 0; m < 5; m++) {
		dtemp[m] =  ce[m][0] +
				xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7]
														+ xi*ce[m][10]))) +
														eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8]
																								   + eta*ce[m][11])))+
																								   zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] +
																										   zeta*ce[m][12])));
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void initialize(void) {

#pragma omp parallel
	{
		/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

		/*--------------------------------------------------------------------
c     This subroutine initializes the field variable u using 
c     tri-linear transfinite interpolation of the boundary values     
c-------------------------------------------------------------------*/

		int i, j, k, m, ix, iy, iz;
		element_t xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];

		/*--------------------------------------------------------------------
c  Later (in compute_rhs) we compute 1/u for every element. A few of 
c  the corner elements are not used, but it convenient (and faster) 
c  to compute the whole thing with a simple loop. Make sure those 
c  values are nonzero by initializing the whole thing here. 
c-------------------------------------------------------------------*/
#pragma omp for
		for (i = 0; i < IMAX; i++) {
			for (j = 0; j < IMAX; j++) {
				for (k = 0; k < IMAX; k++) {
					for (m = 0; m < 5; m++) {
						u[i][j][k][m] = one;
					}
				}
			}
		}

		/*--------------------------------------------------------------------
c     first store the "interpolated" values everywhere on the grid    
c-------------------------------------------------------------------*/

#pragma omp for
		for (i = 0; i < grid_points[0]; i++) {
			xi = (element_t)i * dnxm1;

			for (j = 0; j < grid_points[1]; j++) {
				eta = (element_t)j * dnym1;

				for (k = 0; k < grid_points[2]; k++) {
					zeta = (element_t)k * dnzm1;

					for (ix = 0; ix < 2; ix++) {
						exact_solution((element_t)ix, eta, zeta,
								&(Pface[ix][0][0]));
					}

					for (iy = 0; iy < 2; iy++) {
						exact_solution(xi, (element_t)iy , zeta,
								&Pface[iy][1][0]);
					}

					for (iz = 0; iz < 2; iz++) {
						exact_solution(xi, eta, (element_t)iz,
								&Pface[iz][2][0]);
					}

					for (m = 0; m < 5; m++) {
						Pxi   = xi   * Pface[1][0][m] +
								(one-xi)   * Pface[0][0][m];
						Peta  = eta  * Pface[1][1][m] +
								(one-eta)  * Pface[0][1][m];
						Pzeta = zeta * Pface[1][2][m] +
								(one-zeta) * Pface[0][2][m];

						u[i][j][k][m] = Pxi + Peta + Pzeta -
								Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
								Pxi*Peta*Pzeta;
					}
				}
			}
		}

		/*--------------------------------------------------------------------
c     now store the exact values on the boundaries        
c-------------------------------------------------------------------*/

		/*--------------------------------------------------------------------
c     west face                                                  
c-------------------------------------------------------------------*/
		i = 0;
		xi = zero;
#pragma omp for nowait
		for (j = 0; j < grid_points[1]; j++) {
			eta = (element_t)j * dnym1;
			for (k = 0; k < grid_points[2]; k++) {
				zeta = (element_t)k * dnzm1;
				exact_solution(xi, eta, zeta, temp);
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = temp[m];
				}
			}
		}

		/*--------------------------------------------------------------------
c     east face                                                      
c-------------------------------------------------------------------*/

		i = grid_points[0]-1;
		xi = one;
#pragma omp for
		for (j = 0; j < grid_points[1]; j++) {
			eta = (element_t)j * dnym1;
			for (k = 0; k < grid_points[2]; k++) {
				zeta = (element_t)k * dnzm1;
				exact_solution(xi, eta, zeta, temp);
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = temp[m];
				}
			}
		}

		/*--------------------------------------------------------------------
c     south face                                                 
c-------------------------------------------------------------------*/
		j = 0;
		eta = zero;
#pragma omp for nowait
		for (i = 0; i < grid_points[0]; i++) {
			xi = (element_t)i * dnxm1;
			for (k = 0; k < grid_points[2]; k++) {
				zeta = (element_t)k * dnzm1;
				exact_solution(xi, eta, zeta, temp);
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = temp[m];
				}
			}
		}

		/*--------------------------------------------------------------------
c     north face                                    
c-------------------------------------------------------------------*/
		j = grid_points[1]-1;
		eta = one;
#pragma omp for
		for (i = 0; i < grid_points[0]; i++) {
			xi = (element_t)i * dnxm1;
			for (k = 0; k < grid_points[2]; k++) {
				zeta = (element_t)k * dnzm1;
				exact_solution(xi, eta, zeta, temp);
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = temp[m];
				}
			}
		}

		/*--------------------------------------------------------------------
c     bottom face                                       
c-------------------------------------------------------------------*/
		k = 0;
		zeta = zero;
#pragma omp for nowait
		for (i = 0; i < grid_points[0]; i++) {
			xi = (element_t)i *dnxm1;
			for (j = 0; j < grid_points[1]; j++) {
				eta = (element_t)j * dnym1;
				exact_solution(xi, eta, zeta, temp);
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = temp[m];
				}
			}
		}

/*--------------------------------------------------------------------
c     top face     
c-------------------------------------------------------------------*/
		k = grid_points[2]-1;
		zeta = one;
#pragma omp for
		for (i = 0; i < grid_points[0]; i++) {
			xi = (element_t)i * dnxm1;
			for (j = 0; j < grid_points[1]; j++) {
				eta = (element_t)j * dnym1;
				exact_solution(xi, eta, zeta, temp);
				for (m = 0; m < 5; m++) {
					u[i][j][k][m] = temp[m];
				}
			}
		}
	}

}
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void lhsinit(void) {

#pragma omp parallel
	{
		int i, j, k, m, n;

		/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

		/*--------------------------------------------------------------------
c     zero the whole left hand side for starters
c-------------------------------------------------------------------*/
#pragma omp for
		for (i = 0; i < grid_points[0]; i++) {
			for (j = 0; j < grid_points[1]; j++) {
				for (k = 0; k < grid_points[2]; k++) {
					for (m = 0; m < 5; m++) {
						for (n = 0; n < 5; n++) {
							lhs[i][j][k][0][m][n] = zero;
							lhs[i][j][k][1][m][n] = zero;
							lhs[i][j][k][2][m][n] = zero;
						}
					}
				}
			}
		}

		/*--------------------------------------------------------------------
c     next, set all diagonal values to 1. This is overkill, but convenient
c-------------------------------------------------------------------*/
#pragma omp for  
		for (i = 0; i < grid_points[0]; i++) {
			for (j = 0; j < grid_points[1]; j++) {
				for (k = 0; k < grid_points[2]; k++) {
					for (m = 0; m < 5; m++) {
						lhs[i][j][k][1][m][m] = one;
					}
				}
			}
		}
	}

}
/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void lhsx(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     This function computes the left hand side in the xi-direction
c-------------------------------------------------------------------*/

	int i, j, k;

	/*--------------------------------------------------------------------
c     determine a (labeled f) and n jacobians
c-------------------------------------------------------------------*/
#pragma omp for  
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (i = 0; i < grid_points[0]; i++) {

				tmp1 = one / u[i][j][k][0];
				tmp2 = tmp1 * tmp1;
				tmp3 = tmp1 * tmp2;
				/*--------------------------------------------------------------------
c     
c-------------------------------------------------------------------*/
				fjac[ i][ j][ k][0][0] = zero;
				fjac[ i][ j][ k][0][1] = one;
				fjac[ i][ j][ k][0][2] = zero;
				fjac[ i][ j][ k][0][3] = zero;
				fjac[ i][ j][ k][0][4] = zero;

				fjac[ i][ j][ k][1][0] = -(u[i][j][k][1] * tmp2 *
						u[i][j][k][1])
						+ c2 * one/two * (u[i][j][k][1] * u[i][j][k][1]
																	 + u[i][j][k][2] * u[i][j][k][2]
																								  + u[i][j][k][3] * u[i][j][k][3] ) * tmp2;
				fjac[i][j][k][1][1] = ( two - c2 )
					  * ( u[i][j][k][1] / u[i][j][k][0] );
				fjac[i][j][k][1][2] = - c2 * ( u[i][j][k][2] * tmp1 );
				fjac[i][j][k][1][3] = - c2 * ( u[i][j][k][3] * tmp1 );
				fjac[i][j][k][1][4] = c2;

				fjac[i][j][k][2][0] = - ( u[i][j][k][1]*u[i][j][k][2] ) * tmp2;
				fjac[i][j][k][2][1] = u[i][j][k][2] * tmp1;
				fjac[i][j][k][2][2] = u[i][j][k][1] * tmp1;
				fjac[i][j][k][2][3] = zero;
				fjac[i][j][k][2][4] = zero;

				fjac[i][j][k][3][0] = - ( u[i][j][k][1]*u[i][j][k][3] ) * tmp2;
				fjac[i][j][k][3][1] = u[i][j][k][3] * tmp1;
				fjac[i][j][k][3][2] = zero;
				fjac[i][j][k][3][3] = u[i][j][k][1] * tmp1;
				fjac[i][j][k][3][4] = zero;

				fjac[i][j][k][4][0] = ( c2 * ( u[i][j][k][1] * u[i][j][k][1]
																		  + u[i][j][k][2] * u[i][j][k][2]
																									   + u[i][j][k][3] * u[i][j][k][3] ) * tmp2
						- c1 * ( u[i][j][k][4] * tmp1 ) )
						* ( u[i][j][k][1] * tmp1 );
				fjac[i][j][k][4][1] = c1 *  u[i][j][k][4] * tmp1
						- one/two * c2
						* (  three*u[i][j][k][1]*u[i][j][k][1]
															+ u[i][j][k][2]*u[i][j][k][2]
																					   + u[i][j][k][3]*u[i][j][k][3] ) * tmp2;
				fjac[i][j][k][4][2] = - c2 * ( u[i][j][k][2]*u[i][j][k][1] )
					  * tmp2;
				fjac[i][j][k][4][3] = - c2 * ( u[i][j][k][3]*u[i][j][k][1] )
					  * tmp2;
				fjac[i][j][k][4][4] = c1 * ( u[i][j][k][1] * tmp1 );

				njac[i][j][k][0][0] = zero;
				njac[i][j][k][0][1] = zero;
				njac[i][j][k][0][2] = zero;
				njac[i][j][k][0][3] = zero;
				njac[i][j][k][0][4] = zero;

				njac[i][j][k][1][0] = - con43 * c3c4 * tmp2 * u[i][j][k][1];
				njac[i][j][k][1][1] =   con43 * c3c4 * tmp1;
				njac[i][j][k][1][2] =   zero;
				njac[i][j][k][1][3] =   zero;
				njac[i][j][k][1][4] =   zero;

				njac[i][j][k][2][0] = - c3c4 * tmp2 * u[i][j][k][2];
				njac[i][j][k][2][1] =   zero;
				njac[i][j][k][2][2] =   c3c4 * tmp1;
				njac[i][j][k][2][3] =   zero;
				njac[i][j][k][2][4] =   zero;

				njac[i][j][k][3][0] = - c3c4 * tmp2 * u[i][j][k][3];
				njac[i][j][k][3][1] =   zero;
				njac[i][j][k][3][2] =   zero;
				njac[i][j][k][3][3] =   c3c4 * tmp1;
				njac[i][j][k][3][4] =   zero;

				njac[i][j][k][4][0] = - ( con43 * c3c4
						- c1345 ) * tmp3 * (pow2(u[i][j][k][1]))
						- ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][2]))
						- ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][3]))
						- c1345 * tmp2 * u[i][j][k][4];

				njac[i][j][k][4][1] = ( con43 * c3c4
						- c1345 ) * tmp2 * u[i][j][k][1];
				njac[i][j][k][4][2] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][2];
				njac[i][j][k][4][3] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][3];
				njac[i][j][k][4][4] = ( c1345 ) * tmp1;

			}
			/*--------------------------------------------------------------------
c     now jacobians set, so form left hand side in x direction
c-------------------------------------------------------------------*/
			for (i = 1; i < grid_points[0]-1; i++) {

				tmp1 = dt * tx1;
				tmp2 = dt * tx2;

				lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i-1][j][k][0][0]
																	 - tmp1 * njac[i-1][j][k][0][0]
																								 - tmp1 * dx1;
				lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i-1][j][k][0][1]
																	 - tmp1 * njac[i-1][j][k][0][1];
				lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i-1][j][k][0][2]
																	 - tmp1 * njac[i-1][j][k][0][2];
				lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i-1][j][k][0][3]
																	 - tmp1 * njac[i-1][j][k][0][3];
				lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i-1][j][k][0][4]
																	 - tmp1 * njac[i-1][j][k][0][4];

				lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i-1][j][k][1][0]
																	 - tmp1 * njac[i-1][j][k][1][0];
				lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i-1][j][k][1][1]
																	 - tmp1 * njac[i-1][j][k][1][1]
																								 - tmp1 * dx2;
				lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i-1][j][k][1][2]
																	 - tmp1 * njac[i-1][j][k][1][2];
				lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i-1][j][k][1][3]
																	 - tmp1 * njac[i-1][j][k][1][3];
				lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i-1][j][k][1][4]
																	 - tmp1 * njac[i-1][j][k][1][4];

				lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i-1][j][k][2][0]
																	 - tmp1 * njac[i-1][j][k][2][0];
				lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i-1][j][k][2][1]
																	 - tmp1 * njac[i-1][j][k][2][1];
				lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i-1][j][k][2][2]
																	 - tmp1 * njac[i-1][j][k][2][2]
																								 - tmp1 * dx3;
				lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i-1][j][k][2][3]
																	 - tmp1 * njac[i-1][j][k][2][3];
				lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i-1][j][k][2][4]
																	 - tmp1 * njac[i-1][j][k][2][4];

				lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i-1][j][k][3][0]
																	 - tmp1 * njac[i-1][j][k][3][0];
				lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i-1][j][k][3][1]
																	 - tmp1 * njac[i-1][j][k][3][1];
				lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i-1][j][k][3][2]
																	 - tmp1 * njac[i-1][j][k][3][2];
				lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i-1][j][k][3][3]
																	 - tmp1 * njac[i-1][j][k][3][3]
																								 - tmp1 * dx4;
				lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i-1][j][k][3][4]
																	 - tmp1 * njac[i-1][j][k][3][4];

				lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i-1][j][k][4][0]
																	 - tmp1 * njac[i-1][j][k][4][0];
				lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i-1][j][k][4][1]
																	 - tmp1 * njac[i-1][j][k][4][1];
				lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i-1][j][k][4][2]
																	 - tmp1 * njac[i-1][j][k][4][2];
				lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i-1][j][k][4][3]
																	 - tmp1 * njac[i-1][j][k][4][3];
				lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i-1][j][k][4][4]
																	 - tmp1 * njac[i-1][j][k][4][4]
																								 - tmp1 * dx5;

				lhs[i][j][k][BB][0][0] = one
						+ tmp1 * two * njac[i][j][k][0][0]
														+ tmp1 * two * dx1;
				lhs[i][j][k][BB][0][1] = tmp1 * two * njac[i][j][k][0][1];
				lhs[i][j][k][BB][0][2] = tmp1 * two * njac[i][j][k][0][2];
				lhs[i][j][k][BB][0][3] = tmp1 * two * njac[i][j][k][0][3];
				lhs[i][j][k][BB][0][4] = tmp1 * two * njac[i][j][k][0][4];

				lhs[i][j][k][BB][1][0] = tmp1 * two * njac[i][j][k][1][0];
				lhs[i][j][k][BB][1][1] = one
						+ tmp1 * two * njac[i][j][k][1][1]
														+ tmp1 * two * dx2;
				lhs[i][j][k][BB][1][2] = tmp1 * two * njac[i][j][k][1][2];
				lhs[i][j][k][BB][1][3] = tmp1 * two * njac[i][j][k][1][3];
				lhs[i][j][k][BB][1][4] = tmp1 * two * njac[i][j][k][1][4];

				lhs[i][j][k][BB][2][0] = tmp1 * two * njac[i][j][k][2][0];
				lhs[i][j][k][BB][2][1] = tmp1 * two * njac[i][j][k][2][1];
				lhs[i][j][k][BB][2][2] = one
						+ tmp1 * two * njac[i][j][k][2][2]
														+ tmp1 * two * dx3;
				lhs[i][j][k][BB][2][3] = tmp1 * two * njac[i][j][k][2][3];
				lhs[i][j][k][BB][2][4] = tmp1 * two * njac[i][j][k][2][4];

				lhs[i][j][k][BB][3][0] = tmp1 * two * njac[i][j][k][3][0];
				lhs[i][j][k][BB][3][1] = tmp1 * two * njac[i][j][k][3][1];
				lhs[i][j][k][BB][3][2] = tmp1 * two * njac[i][j][k][3][2];
				lhs[i][j][k][BB][3][3] = one
						+ tmp1 * two * njac[i][j][k][3][3]
														+ tmp1 * two * dx4;
				lhs[i][j][k][BB][3][4] = tmp1 * two * njac[i][j][k][3][4];

				lhs[i][j][k][BB][4][0] = tmp1 * two * njac[i][j][k][4][0];
				lhs[i][j][k][BB][4][1] = tmp1 * two * njac[i][j][k][4][1];
				lhs[i][j][k][BB][4][2] = tmp1 * two * njac[i][j][k][4][2];
				lhs[i][j][k][BB][4][3] = tmp1 * two * njac[i][j][k][4][3];
				lhs[i][j][k][BB][4][4] = one
						+ tmp1 * two * njac[i][j][k][4][4]
														+ tmp1 * two * dx5;

				lhs[i][j][k][CC][0][0] =  tmp2 * fjac[i+1][j][k][0][0]
																	- tmp1 * njac[i+1][j][k][0][0]
																								- tmp1 * dx1;
				lhs[i][j][k][CC][0][1] =  tmp2 * fjac[i+1][j][k][0][1]
																	- tmp1 * njac[i+1][j][k][0][1];
				lhs[i][j][k][CC][0][2] =  tmp2 * fjac[i+1][j][k][0][2]
																	- tmp1 * njac[i+1][j][k][0][2];
				lhs[i][j][k][CC][0][3] =  tmp2 * fjac[i+1][j][k][0][3]
																	- tmp1 * njac[i+1][j][k][0][3];
				lhs[i][j][k][CC][0][4] =  tmp2 * fjac[i+1][j][k][0][4]
																	- tmp1 * njac[i+1][j][k][0][4];

				lhs[i][j][k][CC][1][0] =  tmp2 * fjac[i+1][j][k][1][0]
																	- tmp1 * njac[i+1][j][k][1][0];
				lhs[i][j][k][CC][1][1] =  tmp2 * fjac[i+1][j][k][1][1]
																	- tmp1 * njac[i+1][j][k][1][1]
																								- tmp1 * dx2;
				lhs[i][j][k][CC][1][2] =  tmp2 * fjac[i+1][j][k][1][2]
																	- tmp1 * njac[i+1][j][k][1][2];
				lhs[i][j][k][CC][1][3] =  tmp2 * fjac[i+1][j][k][1][3]
																	- tmp1 * njac[i+1][j][k][1][3];
				lhs[i][j][k][CC][1][4] =  tmp2 * fjac[i+1][j][k][1][4]
																	- tmp1 * njac[i+1][j][k][1][4];

				lhs[i][j][k][CC][2][0] =  tmp2 * fjac[i+1][j][k][2][0]
																	- tmp1 * njac[i+1][j][k][2][0];
				lhs[i][j][k][CC][2][1] =  tmp2 * fjac[i+1][j][k][2][1]
																	- tmp1 * njac[i+1][j][k][2][1];
				lhs[i][j][k][CC][2][2] =  tmp2 * fjac[i+1][j][k][2][2]
																	- tmp1 * njac[i+1][j][k][2][2]
																								- tmp1 * dx3;
				lhs[i][j][k][CC][2][3] =  tmp2 * fjac[i+1][j][k][2][3]
																	- tmp1 * njac[i+1][j][k][2][3];
				lhs[i][j][k][CC][2][4] =  tmp2 * fjac[i+1][j][k][2][4]
																	- tmp1 * njac[i+1][j][k][2][4];

				lhs[i][j][k][CC][3][0] =  tmp2 * fjac[i+1][j][k][3][0]
																	- tmp1 * njac[i+1][j][k][3][0];
				lhs[i][j][k][CC][3][1] =  tmp2 * fjac[i+1][j][k][3][1]
																	- tmp1 * njac[i+1][j][k][3][1];
				lhs[i][j][k][CC][3][2] =  tmp2 * fjac[i+1][j][k][3][2]
																	- tmp1 * njac[i+1][j][k][3][2];
				lhs[i][j][k][CC][3][3] =  tmp2 * fjac[i+1][j][k][3][3]
																	- tmp1 * njac[i+1][j][k][3][3]
																								- tmp1 * dx4;
				lhs[i][j][k][CC][3][4] =  tmp2 * fjac[i+1][j][k][3][4]
																	- tmp1 * njac[i+1][j][k][3][4];

				lhs[i][j][k][CC][4][0] =  tmp2 * fjac[i+1][j][k][4][0]
																	- tmp1 * njac[i+1][j][k][4][0];
				lhs[i][j][k][CC][4][1] =  tmp2 * fjac[i+1][j][k][4][1]
																	- tmp1 * njac[i+1][j][k][4][1];
				lhs[i][j][k][CC][4][2] =  tmp2 * fjac[i+1][j][k][4][2]
																	- tmp1 * njac[i+1][j][k][4][2];
				lhs[i][j][k][CC][4][3] =  tmp2 * fjac[i+1][j][k][4][3]
																	- tmp1 * njac[i+1][j][k][4][3];
				lhs[i][j][k][CC][4][4] =  tmp2 * fjac[i+1][j][k][4][4]
																	- tmp1 * njac[i+1][j][k][4][4]
																								- tmp1 * dx5;

			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void lhsy(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     This function computes the left hand side for the three y-factors   
c-------------------------------------------------------------------*/

	int i, j, k;

	/*--------------------------------------------------------------------
c     Compute the indices for storing the tri-diagonal matrix;
c     determine a (labeled f) and n jacobians for cell c
c-------------------------------------------------------------------*/
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 0; j < grid_points[1]; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {

				tmp1 = one / u[i][j][k][0];
				tmp2 = tmp1 * tmp1;
				tmp3 = tmp1 * tmp2;

				fjac[ i][ j][ k][0][0] = zero;
				fjac[ i][ j][ k][0][1] = zero;
				fjac[ i][ j][ k][0][2] = one;
				fjac[ i][ j][ k][0][3] = zero;
				fjac[ i][ j][ k][0][4] = zero;

				fjac[i][j][k][1][0] = - ( u[i][j][k][1]*u[i][j][k][2] )
					  * tmp2;
				fjac[i][j][k][1][1] = u[i][j][k][2] * tmp1;
				fjac[i][j][k][1][2] = u[i][j][k][1] * tmp1;
				fjac[i][j][k][1][3] = zero;
				fjac[i][j][k][1][4] = zero;

				fjac[i][j][k][2][0] = - ( u[i][j][k][2]*u[i][j][k][2]*tmp2)
					  + one/two * c2 * ( (  u[i][j][k][1] * u[i][j][k][1]
																	   + u[i][j][k][2] * u[i][j][k][2]
																									+ u[i][j][k][3] * u[i][j][k][3] )
							  * tmp2 );
				fjac[i][j][k][2][1] = - c2 *  u[i][j][k][1] * tmp1;
				fjac[i][j][k][2][2] = ( two - c2 )
					  *  u[i][j][k][2] * tmp1;
				fjac[i][j][k][2][3] = - c2 * u[i][j][k][3] * tmp1;
				fjac[i][j][k][2][4] = c2;

				fjac[i][j][k][3][0] = - ( u[i][j][k][2]*u[i][j][k][3] )
					  * tmp2;
				fjac[i][j][k][3][1] = zero;
				fjac[i][j][k][3][2] = u[i][j][k][3] * tmp1;
				fjac[i][j][k][3][3] = u[i][j][k][2] * tmp1;
				fjac[i][j][k][3][4] = zero;

				fjac[i][j][k][4][0] = ( c2 * (  u[i][j][k][1] * u[i][j][k][1]
																		   + u[i][j][k][2] * u[i][j][k][2]
																										+ u[i][j][k][3] * u[i][j][k][3] )
						* tmp2
						- c1 * u[i][j][k][4] * tmp1 )
						* u[i][j][k][2] * tmp1;
				fjac[i][j][k][4][1] = - c2 * u[i][j][k][1]*u[i][j][k][2]
																	  * tmp2;
				fjac[i][j][k][4][2] = c1 * u[i][j][k][4] * tmp1
						- one/two * c2
						* ( (  u[i][j][k][1]*u[i][j][k][1]
														+ three * u[i][j][k][2]*u[i][j][k][2]
																						   + u[i][j][k][3]*u[i][j][k][3] )
								* tmp2 );
				fjac[i][j][k][4][3] = - c2 * ( u[i][j][k][2]*u[i][j][k][3] )
					  * tmp2;
				fjac[i][j][k][4][4] = c1 * u[i][j][k][2] * tmp1;

				njac[i][j][k][0][0] = zero;
				njac[i][j][k][0][1] = zero;
				njac[i][j][k][0][2] = zero;
				njac[i][j][k][0][3] = zero;
				njac[i][j][k][0][4] = zero;

				njac[i][j][k][1][0] = - c3c4 * tmp2 * u[i][j][k][1];
				njac[i][j][k][1][1] =   c3c4 * tmp1;
				njac[i][j][k][1][2] =   zero;
				njac[i][j][k][1][3] =   zero;
				njac[i][j][k][1][4] =   zero;

				njac[i][j][k][2][0] = - con43 * c3c4 * tmp2 * u[i][j][k][2];
				njac[i][j][k][2][1] =   zero;
				njac[i][j][k][2][2] =   con43 * c3c4 * tmp1;
				njac[i][j][k][2][3] =   zero;
				njac[i][j][k][2][4] =   zero;

				njac[i][j][k][3][0] = - c3c4 * tmp2 * u[i][j][k][3];
				njac[i][j][k][3][1] =   zero;
				njac[i][j][k][3][2] =   zero;
				njac[i][j][k][3][3] =   c3c4 * tmp1;
				njac[i][j][k][3][4] =   zero;

				njac[i][j][k][4][0] = - (  c3c4
						- c1345 ) * tmp3 * (pow2(u[i][j][k][1]))
						- ( con43 * c3c4
								- c1345 ) * tmp3 * (pow2(u[i][j][k][2]))
								- ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][3]))
								- c1345 * tmp2 * u[i][j][k][4];

				njac[i][j][k][4][1] = (  c3c4 - c1345 ) * tmp2 * u[i][j][k][1];
				njac[i][j][k][4][2] = ( con43 * c3c4
						- c1345 ) * tmp2 * u[i][j][k][2];
				njac[i][j][k][4][3] = ( c3c4 - c1345 ) * tmp2 * u[i][j][k][3];
				njac[i][j][k][4][4] = ( c1345 ) * tmp1;

			}
		}
	}

	/*--------------------------------------------------------------------
c     now joacobians set, so form left hand side in y direction
c-------------------------------------------------------------------*/
#pragma omp for  
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {

				tmp1 = dt * ty1;
				tmp2 = dt * ty2;

				lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i][j-1][k][0][0]
																	 - tmp1 * njac[i][j-1][k][0][0]
																								 - tmp1 * dy1;
				lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i][j-1][k][0][1]
																	 - tmp1 * njac[i][j-1][k][0][1];
				lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i][j-1][k][0][2]
																	 - tmp1 * njac[i][j-1][k][0][2];
				lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i][j-1][k][0][3]
																	 - tmp1 * njac[i][j-1][k][0][3];
				lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i][j-1][k][0][4]
																	 - tmp1 * njac[i][j-1][k][0][4];

				lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i][j-1][k][1][0]
																	 - tmp1 * njac[i][j-1][k][1][0];
				lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i][j-1][k][1][1]
																	 - tmp1 * njac[i][j-1][k][1][1]
																								 - tmp1 * dy2;
				lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i][j-1][k][1][2]
																	 - tmp1 * njac[i][j-1][k][1][2];
				lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i][j-1][k][1][3]
																	 - tmp1 * njac[i][j-1][k][1][3];
				lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i][j-1][k][1][4]
																	 - tmp1 * njac[i][j-1][k][1][4];

				lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i][j-1][k][2][0]
																	 - tmp1 * njac[i][j-1][k][2][0];
				lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i][j-1][k][2][1]
																	 - tmp1 * njac[i][j-1][k][2][1];
				lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i][j-1][k][2][2]
																	 - tmp1 * njac[i][j-1][k][2][2]
																								 - tmp1 * dy3;
				lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i][j-1][k][2][3]
																	 - tmp1 * njac[i][j-1][k][2][3];
				lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i][j-1][k][2][4]
																	 - tmp1 * njac[i][j-1][k][2][4];

				lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i][j-1][k][3][0]
																	 - tmp1 * njac[i][j-1][k][3][0];
				lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i][j-1][k][3][1]
																	 - tmp1 * njac[i][j-1][k][3][1];
				lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i][j-1][k][3][2]
																	 - tmp1 * njac[i][j-1][k][3][2];
				lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i][j-1][k][3][3]
																	 - tmp1 * njac[i][j-1][k][3][3]
																								 - tmp1 * dy4;
				lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i][j-1][k][3][4]
																	 - tmp1 * njac[i][j-1][k][3][4];

				lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i][j-1][k][4][0]
																	 - tmp1 * njac[i][j-1][k][4][0];
				lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i][j-1][k][4][1]
																	 - tmp1 * njac[i][j-1][k][4][1];
				lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i][j-1][k][4][2]
																	 - tmp1 * njac[i][j-1][k][4][2];
				lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i][j-1][k][4][3]
																	 - tmp1 * njac[i][j-1][k][4][3];
				lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i][j-1][k][4][4]
																	 - tmp1 * njac[i][j-1][k][4][4]
																								 - tmp1 * dy5;

				lhs[i][j][k][BB][0][0] = one
						+ tmp1 * two * njac[i][j][k][0][0]
														+ tmp1 * two * dy1;
				lhs[i][j][k][BB][0][1] = tmp1 * two * njac[i][j][k][0][1];
				lhs[i][j][k][BB][0][2] = tmp1 * two * njac[i][j][k][0][2];
				lhs[i][j][k][BB][0][3] = tmp1 * two * njac[i][j][k][0][3];
				lhs[i][j][k][BB][0][4] = tmp1 * two * njac[i][j][k][0][4];

				lhs[i][j][k][BB][1][0] = tmp1 * two * njac[i][j][k][1][0];
				lhs[i][j][k][BB][1][1] = one
						+ tmp1 * two * njac[i][j][k][1][1]
														+ tmp1 * two * dy2;
				lhs[i][j][k][BB][1][2] = tmp1 * two * njac[i][j][k][1][2];
				lhs[i][j][k][BB][1][3] = tmp1 * two * njac[i][j][k][1][3];
				lhs[i][j][k][BB][1][4] = tmp1 * two * njac[i][j][k][1][4];

				lhs[i][j][k][BB][2][0] = tmp1 * two * njac[i][j][k][2][0];
				lhs[i][j][k][BB][2][1] = tmp1 * two * njac[i][j][k][2][1];
				lhs[i][j][k][BB][2][2] = one
						+ tmp1 * two * njac[i][j][k][2][2]
														+ tmp1 * two * dy3;
				lhs[i][j][k][BB][2][3] = tmp1 * two * njac[i][j][k][2][3];
				lhs[i][j][k][BB][2][4] = tmp1 * two * njac[i][j][k][2][4];

				lhs[i][j][k][BB][3][0] = tmp1 * two * njac[i][j][k][3][0];
				lhs[i][j][k][BB][3][1] = tmp1 * two * njac[i][j][k][3][1];
				lhs[i][j][k][BB][3][2] = tmp1 * two * njac[i][j][k][3][2];
				lhs[i][j][k][BB][3][3] = one
						+ tmp1 * two * njac[i][j][k][3][3]
														+ tmp1 * two * dy4;
				lhs[i][j][k][BB][3][4] = tmp1 * two * njac[i][j][k][3][4];

				lhs[i][j][k][BB][4][0] = tmp1 * two * njac[i][j][k][4][0];
				lhs[i][j][k][BB][4][1] = tmp1 * two * njac[i][j][k][4][1];
				lhs[i][j][k][BB][4][2] = tmp1 * two * njac[i][j][k][4][2];
				lhs[i][j][k][BB][4][3] = tmp1 * two * njac[i][j][k][4][3];
				lhs[i][j][k][BB][4][4] = one
						+ tmp1 * two * njac[i][j][k][4][4]
														+ tmp1 * two * dy5;

				lhs[i][j][k][CC][0][0] =  tmp2 * fjac[i][j+1][k][0][0]
																	- tmp1 * njac[i][j+1][k][0][0]
																								- tmp1 * dy1;
				lhs[i][j][k][CC][0][1] =  tmp2 * fjac[i][j+1][k][0][1]
																	- tmp1 * njac[i][j+1][k][0][1];
				lhs[i][j][k][CC][0][2] =  tmp2 * fjac[i][j+1][k][0][2]
																	- tmp1 * njac[i][j+1][k][0][2];
				lhs[i][j][k][CC][0][3] =  tmp2 * fjac[i][j+1][k][0][3]
																	- tmp1 * njac[i][j+1][k][0][3];
				lhs[i][j][k][CC][0][4] =  tmp2 * fjac[i][j+1][k][0][4]
																	- tmp1 * njac[i][j+1][k][0][4];

				lhs[i][j][k][CC][1][0] =  tmp2 * fjac[i][j+1][k][1][0]
																	- tmp1 * njac[i][j+1][k][1][0];
				lhs[i][j][k][CC][1][1] =  tmp2 * fjac[i][j+1][k][1][1]
																	- tmp1 * njac[i][j+1][k][1][1]
																								- tmp1 * dy2;
				lhs[i][j][k][CC][1][2] =  tmp2 * fjac[i][j+1][k][1][2]
																	- tmp1 * njac[i][j+1][k][1][2];
				lhs[i][j][k][CC][1][3] =  tmp2 * fjac[i][j+1][k][1][3]
																	- tmp1 * njac[i][j+1][k][1][3];
				lhs[i][j][k][CC][1][4] =  tmp2 * fjac[i][j+1][k][1][4]
																	- tmp1 * njac[i][j+1][k][1][4];

				lhs[i][j][k][CC][2][0] =  tmp2 * fjac[i][j+1][k][2][0]
																	- tmp1 * njac[i][j+1][k][2][0];
				lhs[i][j][k][CC][2][1] =  tmp2 * fjac[i][j+1][k][2][1]
																	- tmp1 * njac[i][j+1][k][2][1];
				lhs[i][j][k][CC][2][2] =  tmp2 * fjac[i][j+1][k][2][2]
																	- tmp1 * njac[i][j+1][k][2][2]
																								- tmp1 * dy3;
				lhs[i][j][k][CC][2][3] =  tmp2 * fjac[i][j+1][k][2][3]
																	- tmp1 * njac[i][j+1][k][2][3];
				lhs[i][j][k][CC][2][4] =  tmp2 * fjac[i][j+1][k][2][4]
																	- tmp1 * njac[i][j+1][k][2][4];

				lhs[i][j][k][CC][3][0] =  tmp2 * fjac[i][j+1][k][3][0]
																	- tmp1 * njac[i][j+1][k][3][0];
				lhs[i][j][k][CC][3][1] =  tmp2 * fjac[i][j+1][k][3][1]
																	- tmp1 * njac[i][j+1][k][3][1];
				lhs[i][j][k][CC][3][2] =  tmp2 * fjac[i][j+1][k][3][2]
																	- tmp1 * njac[i][j+1][k][3][2];
				lhs[i][j][k][CC][3][3] =  tmp2 * fjac[i][j+1][k][3][3]
																	- tmp1 * njac[i][j+1][k][3][3]
																								- tmp1 * dy4;
				lhs[i][j][k][CC][3][4] =  tmp2 * fjac[i][j+1][k][3][4]
																	- tmp1 * njac[i][j+1][k][3][4];

				lhs[i][j][k][CC][4][0] =  tmp2 * fjac[i][j+1][k][4][0]
																	- tmp1 * njac[i][j+1][k][4][0];
				lhs[i][j][k][CC][4][1] =  tmp2 * fjac[i][j+1][k][4][1]
																	- tmp1 * njac[i][j+1][k][4][1];
				lhs[i][j][k][CC][4][2] =  tmp2 * fjac[i][j+1][k][4][2]
																	- tmp1 * njac[i][j+1][k][4][2];
				lhs[i][j][k][CC][4][3] =  tmp2 * fjac[i][j+1][k][4][3]
																	- tmp1 * njac[i][j+1][k][4][3];
				lhs[i][j][k][CC][4][4] =  tmp2 * fjac[i][j+1][k][4][4]
																	- tmp1 * njac[i][j+1][k][4][4]
																								- tmp1 * dy5;

			}
		}
	}
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void lhsz(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     This function computes the left hand side for the three z-factors   
c-------------------------------------------------------------------*/

	int i, j, k;

	/*--------------------------------------------------------------------
c     Compute the indices for storing the block-diagonal matrix;
c     determine c (labeled f) and s jacobians
c---------------------------------------------------------------------*/
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 0; k < grid_points[2]; k++) {

				tmp1 = one / u[i][j][k][0];
				tmp2 = tmp1 * tmp1;
				tmp3 = tmp1 * tmp2;

				fjac[i][j][k][0][0] = zero;
				fjac[i][j][k][0][1] = zero;
				fjac[i][j][k][0][2] = zero;
				fjac[i][j][k][0][3] = one;
				fjac[i][j][k][0][4] = zero;

				fjac[i][j][k][1][0] = - ( u[i][j][k][1]*u[i][j][k][3] )
					  * tmp2;
				fjac[i][j][k][1][1] = u[i][j][k][3] * tmp1;
				fjac[i][j][k][1][2] = zero;
				fjac[i][j][k][1][3] = u[i][j][k][1] * tmp1;
				fjac[i][j][k][1][4] = zero;

				fjac[i][j][k][2][0] = - ( u[i][j][k][2]*u[i][j][k][3] )
					  * tmp2;
				fjac[i][j][k][2][1] = zero;
				fjac[i][j][k][2][2] = u[i][j][k][3] * tmp1;
				fjac[i][j][k][2][3] = u[i][j][k][2] * tmp1;
				fjac[i][j][k][2][4] = zero;

				fjac[i][j][k][3][0] = - (u[i][j][k][3]*u[i][j][k][3] * tmp2 )
					  + one/two * c2 * ( (  u[i][j][k][1] * u[i][j][k][1]
																	   + u[i][j][k][2] * u[i][j][k][2]
																									+ u[i][j][k][3] * u[i][j][k][3] ) * tmp2 );
				fjac[i][j][k][3][1] = - c2 *  u[i][j][k][1] * tmp1;
				fjac[i][j][k][3][2] = - c2 *  u[i][j][k][2] * tmp1;
				fjac[i][j][k][3][3] = ( two - c2 )
					  *  u[i][j][k][3] * tmp1;
				fjac[i][j][k][3][4] = c2;

				fjac[i][j][k][4][0] = ( c2 * (  u[i][j][k][1] * u[i][j][k][1]
																		   + u[i][j][k][2] * u[i][j][k][2]
																										+ u[i][j][k][3] * u[i][j][k][3] )
						* tmp2
						- c1 * ( u[i][j][k][4] * tmp1 ) )
						* ( u[i][j][k][3] * tmp1 );
				fjac[i][j][k][4][1] = - c2 * ( u[i][j][k][1]*u[i][j][k][3] )
					  * tmp2;
				fjac[i][j][k][4][2] = - c2 * ( u[i][j][k][2]*u[i][j][k][3] )
					  * tmp2;
				fjac[i][j][k][4][3] = c1 * ( u[i][j][k][4] * tmp1 )
					  - one/two * c2
					  * ( (  u[i][j][k][1]*u[i][j][k][1]
													  + u[i][j][k][2]*u[i][j][k][2]
																				 + three*u[i][j][k][3]*u[i][j][k][3] )
							  * tmp2 );
				fjac[i][j][k][4][4] = c1 * u[i][j][k][3] * tmp1;

				njac[i][j][k][0][0] = zero;
				njac[i][j][k][0][1] = zero;
				njac[i][j][k][0][2] = zero;
				njac[i][j][k][0][3] = zero;
				njac[i][j][k][0][4] = zero;

				njac[i][j][k][1][0] = - c3c4 * tmp2 * u[i][j][k][1];
				njac[i][j][k][1][1] =   c3c4 * tmp1;
				njac[i][j][k][1][2] =   zero;
				njac[i][j][k][1][3] =   zero;
				njac[i][j][k][1][4] =   zero;

				njac[i][j][k][2][0] = - c3c4 * tmp2 * u[i][j][k][2];
				njac[i][j][k][2][1] =   zero;
				njac[i][j][k][2][2] =   c3c4 * tmp1;
				njac[i][j][k][2][3] =   zero;
				njac[i][j][k][2][4] =   zero;

				njac[i][j][k][3][0] = - con43 * c3c4 * tmp2 * u[i][j][k][3];
				njac[i][j][k][3][1] =   zero;
				njac[i][j][k][3][2] =   zero;
				njac[i][j][k][3][3] =   con43 * c3 * c4 * tmp1;
				njac[i][j][k][3][4] =   zero;

				njac[i][j][k][4][0] = - (  c3c4
						- c1345 ) * tmp3 * (pow2(u[i][j][k][1]))
						- ( c3c4 - c1345 ) * tmp3 * (pow2(u[i][j][k][2]))
						- ( con43 * c3c4
								- c1345 ) * tmp3 * (pow2(u[i][j][k][3]))
								- c1345 * tmp2 * u[i][j][k][4];

				njac[i][j][k][4][1] = (  c3c4 - c1345 ) * tmp2 * u[i][j][k][1];
				njac[i][j][k][4][2] = (  c3c4 - c1345 ) * tmp2 * u[i][j][k][2];
				njac[i][j][k][4][3] = ( con43 * c3c4
						- c1345 ) * tmp2 * u[i][j][k][3];
				njac[i][j][k][4][4] = ( c1345 )* tmp1;

			}
		}
	}

	/*--------------------------------------------------------------------
c     now jacobians set, so form left hand side in z direction
c-------------------------------------------------------------------*/
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {

				tmp1 = dt * tz1;
				tmp2 = dt * tz2;

				lhs[i][j][k][AA][0][0] = - tmp2 * fjac[i][j][k-1][0][0]
																	 - tmp1 * njac[i][j][k-1][0][0]
																								 - tmp1 * dz1;
				lhs[i][j][k][AA][0][1] = - tmp2 * fjac[i][j][k-1][0][1]
																	 - tmp1 * njac[i][j][k-1][0][1];
				lhs[i][j][k][AA][0][2] = - tmp2 * fjac[i][j][k-1][0][2]
																	 - tmp1 * njac[i][j][k-1][0][2];
				lhs[i][j][k][AA][0][3] = - tmp2 * fjac[i][j][k-1][0][3]
																	 - tmp1 * njac[i][j][k-1][0][3];
				lhs[i][j][k][AA][0][4] = - tmp2 * fjac[i][j][k-1][0][4]
																	 - tmp1 * njac[i][j][k-1][0][4];

				lhs[i][j][k][AA][1][0] = - tmp2 * fjac[i][j][k-1][1][0]
																	 - tmp1 * njac[i][j][k-1][1][0];
				lhs[i][j][k][AA][1][1] = - tmp2 * fjac[i][j][k-1][1][1]
																	 - tmp1 * njac[i][j][k-1][1][1]
																								 - tmp1 * dz2;
				lhs[i][j][k][AA][1][2] = - tmp2 * fjac[i][j][k-1][1][2]
																	 - tmp1 * njac[i][j][k-1][1][2];
				lhs[i][j][k][AA][1][3] = - tmp2 * fjac[i][j][k-1][1][3]
																	 - tmp1 * njac[i][j][k-1][1][3];
				lhs[i][j][k][AA][1][4] = - tmp2 * fjac[i][j][k-1][1][4]
																	 - tmp1 * njac[i][j][k-1][1][4];

				lhs[i][j][k][AA][2][0] = - tmp2 * fjac[i][j][k-1][2][0]
																	 - tmp1 * njac[i][j][k-1][2][0];
				lhs[i][j][k][AA][2][1] = - tmp2 * fjac[i][j][k-1][2][1]
																	 - tmp1 * njac[i][j][k-1][2][1];
				lhs[i][j][k][AA][2][2] = - tmp2 * fjac[i][j][k-1][2][2]
																	 - tmp1 * njac[i][j][k-1][2][2]
																								 - tmp1 * dz3;
				lhs[i][j][k][AA][2][3] = - tmp2 * fjac[i][j][k-1][2][3]
																	 - tmp1 * njac[i][j][k-1][2][3];
				lhs[i][j][k][AA][2][4] = - tmp2 * fjac[i][j][k-1][2][4]
																	 - tmp1 * njac[i][j][k-1][2][4];

				lhs[i][j][k][AA][3][0] = - tmp2 * fjac[i][j][k-1][3][0]
																	 - tmp1 * njac[i][j][k-1][3][0];
				lhs[i][j][k][AA][3][1] = - tmp2 * fjac[i][j][k-1][3][1]
																	 - tmp1 * njac[i][j][k-1][3][1];
				lhs[i][j][k][AA][3][2] = - tmp2 * fjac[i][j][k-1][3][2]
																	 - tmp1 * njac[i][j][k-1][3][2];
				lhs[i][j][k][AA][3][3] = - tmp2 * fjac[i][j][k-1][3][3]
																	 - tmp1 * njac[i][j][k-1][3][3]
																								 - tmp1 * dz4;
				lhs[i][j][k][AA][3][4] = - tmp2 * fjac[i][j][k-1][3][4]
																	 - tmp1 * njac[i][j][k-1][3][4];

				lhs[i][j][k][AA][4][0] = - tmp2 * fjac[i][j][k-1][4][0]
																	 - tmp1 * njac[i][j][k-1][4][0];
				lhs[i][j][k][AA][4][1] = - tmp2 * fjac[i][j][k-1][4][1]
																	 - tmp1 * njac[i][j][k-1][4][1];
				lhs[i][j][k][AA][4][2] = - tmp2 * fjac[i][j][k-1][4][2]
																	 - tmp1 * njac[i][j][k-1][4][2];
				lhs[i][j][k][AA][4][3] = - tmp2 * fjac[i][j][k-1][4][3]
																	 - tmp1 * njac[i][j][k-1][4][3];
				lhs[i][j][k][AA][4][4] = - tmp2 * fjac[i][j][k-1][4][4]
																	 - tmp1 * njac[i][j][k-1][4][4]
																								 - tmp1 * dz5;

				lhs[i][j][k][BB][0][0] = one
						+ tmp1 * two * njac[i][j][k][0][0]
														+ tmp1 * two * dz1;
				lhs[i][j][k][BB][0][1] = tmp1 * two * njac[i][j][k][0][1];
				lhs[i][j][k][BB][0][2] = tmp1 * two * njac[i][j][k][0][2];
				lhs[i][j][k][BB][0][3] = tmp1 * two * njac[i][j][k][0][3];
				lhs[i][j][k][BB][0][4] = tmp1 * two * njac[i][j][k][0][4];

				lhs[i][j][k][BB][1][0] = tmp1 * two * njac[i][j][k][1][0];
				lhs[i][j][k][BB][1][1] = one
						+ tmp1 * two * njac[i][j][k][1][1]
														+ tmp1 * two * dz2;
				lhs[i][j][k][BB][1][2] = tmp1 * two * njac[i][j][k][1][2];
				lhs[i][j][k][BB][1][3] = tmp1 * two * njac[i][j][k][1][3];
				lhs[i][j][k][BB][1][4] = tmp1 * two * njac[i][j][k][1][4];

				lhs[i][j][k][BB][2][0] = tmp1 * two * njac[i][j][k][2][0];
				lhs[i][j][k][BB][2][1] = tmp1 * two * njac[i][j][k][2][1];
				lhs[i][j][k][BB][2][2] = one
						+ tmp1 * two * njac[i][j][k][2][2]
														+ tmp1 * two * dz3;
				lhs[i][j][k][BB][2][3] = tmp1 * two * njac[i][j][k][2][3];
				lhs[i][j][k][BB][2][4] = tmp1 * two * njac[i][j][k][2][4];

				lhs[i][j][k][BB][3][0] = tmp1 * two * njac[i][j][k][3][0];
				lhs[i][j][k][BB][3][1] = tmp1 * two * njac[i][j][k][3][1];
				lhs[i][j][k][BB][3][2] = tmp1 * two * njac[i][j][k][3][2];
				lhs[i][j][k][BB][3][3] = one
						+ tmp1 * two * njac[i][j][k][3][3]
														+ tmp1 * two * dz4;
				lhs[i][j][k][BB][3][4] = tmp1 * two * njac[i][j][k][3][4];

				lhs[i][j][k][BB][4][0] = tmp1 * two * njac[i][j][k][4][0];
				lhs[i][j][k][BB][4][1] = tmp1 * two * njac[i][j][k][4][1];
				lhs[i][j][k][BB][4][2] = tmp1 * two * njac[i][j][k][4][2];
				lhs[i][j][k][BB][4][3] = tmp1 * two * njac[i][j][k][4][3];
				lhs[i][j][k][BB][4][4] = one
						+ tmp1 * two * njac[i][j][k][4][4]
														+ tmp1 * two * dz5;

				lhs[i][j][k][CC][0][0] =  tmp2 * fjac[i][j][k+1][0][0]
																	- tmp1 * njac[i][j][k+1][0][0]
																								- tmp1 * dz1;
				lhs[i][j][k][CC][0][1] =  tmp2 * fjac[i][j][k+1][0][1]
																	- tmp1 * njac[i][j][k+1][0][1];
				lhs[i][j][k][CC][0][2] =  tmp2 * fjac[i][j][k+1][0][2]
																	- tmp1 * njac[i][j][k+1][0][2];
				lhs[i][j][k][CC][0][3] =  tmp2 * fjac[i][j][k+1][0][3]
																	- tmp1 * njac[i][j][k+1][0][3];
				lhs[i][j][k][CC][0][4] =  tmp2 * fjac[i][j][k+1][0][4]
																	- tmp1 * njac[i][j][k+1][0][4];

				lhs[i][j][k][CC][1][0] =  tmp2 * fjac[i][j][k+1][1][0]
																	- tmp1 * njac[i][j][k+1][1][0];
				lhs[i][j][k][CC][1][1] =  tmp2 * fjac[i][j][k+1][1][1]
																	- tmp1 * njac[i][j][k+1][1][1]
																								- tmp1 * dz2;
				lhs[i][j][k][CC][1][2] =  tmp2 * fjac[i][j][k+1][1][2]
																	- tmp1 * njac[i][j][k+1][1][2];
				lhs[i][j][k][CC][1][3] =  tmp2 * fjac[i][j][k+1][1][3]
																	- tmp1 * njac[i][j][k+1][1][3];
				lhs[i][j][k][CC][1][4] =  tmp2 * fjac[i][j][k+1][1][4]
																	- tmp1 * njac[i][j][k+1][1][4];

				lhs[i][j][k][CC][2][0] =  tmp2 * fjac[i][j][k+1][2][0]
																	- tmp1 * njac[i][j][k+1][2][0];
				lhs[i][j][k][CC][2][1] =  tmp2 * fjac[i][j][k+1][2][1]
																	- tmp1 * njac[i][j][k+1][2][1];
				lhs[i][j][k][CC][2][2] =  tmp2 * fjac[i][j][k+1][2][2]
																	- tmp1 * njac[i][j][k+1][2][2]
																								- tmp1 * dz3;
				lhs[i][j][k][CC][2][3] =  tmp2 * fjac[i][j][k+1][2][3]
																	- tmp1 * njac[i][j][k+1][2][3];
				lhs[i][j][k][CC][2][4] =  tmp2 * fjac[i][j][k+1][2][4]
																	- tmp1 * njac[i][j][k+1][2][4];

				lhs[i][j][k][CC][3][0] =  tmp2 * fjac[i][j][k+1][3][0]
																	- tmp1 * njac[i][j][k+1][3][0];
				lhs[i][j][k][CC][3][1] =  tmp2 * fjac[i][j][k+1][3][1]
																	- tmp1 * njac[i][j][k+1][3][1];
				lhs[i][j][k][CC][3][2] =  tmp2 * fjac[i][j][k+1][3][2]
																	- tmp1 * njac[i][j][k+1][3][2];
				lhs[i][j][k][CC][3][3] =  tmp2 * fjac[i][j][k+1][3][3]
																	- tmp1 * njac[i][j][k+1][3][3]
																								- tmp1 * dz4;
				lhs[i][j][k][CC][3][4] =  tmp2 * fjac[i][j][k+1][3][4]
																	- tmp1 * njac[i][j][k+1][3][4];

				lhs[i][j][k][CC][4][0] =  tmp2 * fjac[i][j][k+1][4][0]
																	- tmp1 * njac[i][j][k+1][4][0];
				lhs[i][j][k][CC][4][1] =  tmp2 * fjac[i][j][k+1][4][1]
																	- tmp1 * njac[i][j][k+1][4][1];
				lhs[i][j][k][CC][4][2] =  tmp2 * fjac[i][j][k+1][4][2]
																	- tmp1 * njac[i][j][k+1][4][2];
				lhs[i][j][k][CC][4][3] =  tmp2 * fjac[i][j][k+1][4][3]
																	- tmp1 * njac[i][j][k+1][4][3];
				lhs[i][j][k][CC][4][4] =  tmp2 * fjac[i][j][k+1][4][4]
																	- tmp1 * njac[i][j][k+1][4][4]
																								- tmp1 * dz5;

			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void compute_rhs(void) {

	int i, j, k, m;
	element_t rho_inv, uijk, up1, um1, vijk, vp1, vm1, wijk, wp1, wm1;

	/*--------------------------------------------------------------------
c     compute the reciprocal of density, and the kinetic energy, 
c     and the speed of sound.
c-------------------------------------------------------------------*/
#pragma omp for nowait
	for (i = 0; i < grid_points[0]; i++) {
		for (j = 0; j < grid_points[1]; j++) {
			for (k = 0; k < grid_points[2]; k++) {
				rho_inv = one/u[i][j][k][0];
				rho_i[i][j][k] = rho_inv;
				us[i][j][k] = u[i][j][k][1] * rho_inv;
				vs[i][j][k] = u[i][j][k][2] * rho_inv;
				ws[i][j][k] = u[i][j][k][3] * rho_inv;
				square[i][j][k] = one/two * (u[i][j][k][1]*u[i][j][k][1] +
						u[i][j][k][2]*u[i][j][k][2] +
						u[i][j][k][3]*u[i][j][k][3] ) * rho_inv;
				qs[i][j][k] = square[i][j][k] * rho_inv;
			}
		}
	}

	/*--------------------------------------------------------------------
c copy the exact forcing term to the right hand side;  because 
c this forcing term is known, we can store it on the whole grid
c including the boundary                   
c-------------------------------------------------------------------*/

#pragma omp for
	for (i = 0; i < grid_points[0]; i++) {
		for (j = 0; j < grid_points[1]; j++) {
			for (k = 0; k < grid_points[2]; k++) {
				for (m = 0; m < 5; m++) {
					rhs[i][j][k][m] = forcing[i][j][k][m];
				}
			}
		}
	}

	/*--------------------------------------------------------------------
c     compute xi-direction fluxes 
c-------------------------------------------------------------------*/
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				uijk = us[i][j][k];
				up1  = us[i+1][j][k];
				um1  = us[i-1][j][k];

				rhs[i][j][k][0] = rhs[i][j][k][0] + dx1tx1 *
						(u[i+1][j][k][0] - two*u[i][j][k][0] +
								u[i-1][j][k][0]) -
								tx2 * (u[i+1][j][k][1] - u[i-1][j][k][1]);

				rhs[i][j][k][1] = rhs[i][j][k][1] + dx2tx1 *
						(u[i+1][j][k][1] - two*u[i][j][k][1] +
								u[i-1][j][k][1]) +
								xxcon2*con43 * (up1 - two*uijk + um1) -
								tx2 * (u[i+1][j][k][1]*up1 -
										u[i-1][j][k][1]*um1 +
										(u[i+1][j][k][4]- square[i+1][j][k]-
												u[i-1][j][k][4]+ square[i-1][j][k])*
												c2);

				rhs[i][j][k][2] = rhs[i][j][k][2] + dx3tx1 *
						(u[i+1][j][k][2] - two*u[i][j][k][2] +
								u[i-1][j][k][2]) +
								xxcon2 * (vs[i+1][j][k] - two*vs[i][j][k] +
										vs[i-1][j][k]) -
										tx2 * (u[i+1][j][k][2]*up1 -
												u[i-1][j][k][2]*um1);

				rhs[i][j][k][3] = rhs[i][j][k][3] + dx4tx1 *
						(u[i+1][j][k][3] - two*u[i][j][k][3] +
								u[i-1][j][k][3]) +
								xxcon2 * (ws[i+1][j][k] - two*ws[i][j][k] +
										ws[i-1][j][k]) -
										tx2 * (u[i+1][j][k][3]*up1 -
												u[i-1][j][k][3]*um1);

				rhs[i][j][k][4] = rhs[i][j][k][4] + dx5tx1 *
						(u[i+1][j][k][4] - two*u[i][j][k][4] +
								u[i-1][j][k][4]) +
								xxcon3 * (qs[i+1][j][k] - two*qs[i][j][k] +
										qs[i-1][j][k]) +
										xxcon4 * (up1*up1 -       two*uijk*uijk +
												um1*um1) +
												xxcon5 * (u[i+1][j][k][4]*rho_i[i+1][j][k] -
														two*u[i][j][k][4]*rho_i[i][j][k] +
														u[i-1][j][k][4]*rho_i[i-1][j][k]) -
														tx2 * ( (c1*u[i+1][j][k][4] -
																c2*square[i+1][j][k])*up1 -
																(c1*u[i-1][j][k][4] -
																		c2*square[i-1][j][k])*um1 );
			}
		}
	}

	/*--------------------------------------------------------------------
c     add fourth order xi-direction dissipation               
c-------------------------------------------------------------------*/
	i = 1;
#pragma omp for nowait
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m]- dssp *
						( five*u[i][j][k][m] - four*u[i+1][j][k][m] +
								u[i+2][j][k][m]);
			}
		}
	}

	i = 2;
#pragma omp for nowait
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						(-four*u[i-1][j][k][m] + six*u[i][j][k][m] -
								four*u[i+1][j][k][m] + u[i+2][j][k][m]);
			}
		}
	}

#pragma omp for nowait
	for (i = 3; i < grid_points[0]-3; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				for (m = 0; m < 5; m++) {
					rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
							(  u[i-2][j][k][m] - four*u[i-1][j][k][m] +
									six*u[i][j][k][m] - four*u[i+1][j][k][m] +
									u[i+2][j][k][m] );
				}
			}
		}
	}

	i = grid_points[0]-3;
#pragma omp for nowait
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						( u[i-2][j][k][m] - four*u[i-1][j][k][m] +
								six*u[i][j][k][m] - four*u[i+1][j][k][m] );
			}
		}
	}

	i = grid_points[0]-2;
#pragma omp for
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						( u[i-2][j][k][m] - four*u[i-1][j][k][m] +
								five*u[i][j][k][m] );
			}
		}
	}

	/*--------------------------------------------------------------------
c     compute eta-direction fluxes 
c-------------------------------------------------------------------*/
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				vijk = vs[i][j][k];
				vp1  = vs[i][j+1][k];
				vm1  = vs[i][j-1][k];
				rhs[i][j][k][0] = rhs[i][j][k][0] + dy1ty1 *
						(u[i][j+1][k][0] - two*u[i][j][k][0] +
								u[i][j-1][k][0]) -
								ty2 * (u[i][j+1][k][2] - u[i][j-1][k][2]);
				rhs[i][j][k][1] = rhs[i][j][k][1] + dy2ty1 *
						(u[i][j+1][k][1] - two*u[i][j][k][1] +
								u[i][j-1][k][1]) +
								yycon2 * (us[i][j+1][k] - two*us[i][j][k] +
										us[i][j-1][k]) -
										ty2 * (u[i][j+1][k][1]*vp1 -
												u[i][j-1][k][1]*vm1);
				rhs[i][j][k][2] = rhs[i][j][k][2] + dy3ty1 *
						(u[i][j+1][k][2] - two*u[i][j][k][2] +
								u[i][j-1][k][2]) +
								yycon2*con43 * (vp1 - two*vijk + vm1) -
								ty2 * (u[i][j+1][k][2]*vp1 -
										u[i][j-1][k][2]*vm1 +
										(u[i][j+1][k][4] - square[i][j+1][k] -
												u[i][j-1][k][4] + square[i][j-1][k])
												*c2);
				rhs[i][j][k][3] = rhs[i][j][k][3] + dy4ty1 *
						(u[i][j+1][k][3] - two*u[i][j][k][3] +
								u[i][j-1][k][3]) +
								yycon2 * (ws[i][j+1][k] - two*ws[i][j][k] +
										ws[i][j-1][k]) -
										ty2 * (u[i][j+1][k][3]*vp1 -
												u[i][j-1][k][3]*vm1);
				rhs[i][j][k][4] = rhs[i][j][k][4] + dy5ty1 *
						(u[i][j+1][k][4] - two*u[i][j][k][4] +
								u[i][j-1][k][4]) +
								yycon3 * (qs[i][j+1][k] - two*qs[i][j][k] +
										qs[i][j-1][k]) +
										yycon4 * (vp1*vp1       - two*vijk*vijk +
												vm1*vm1) +
												yycon5 * (u[i][j+1][k][4]*rho_i[i][j+1][k] -
														two*u[i][j][k][4]*rho_i[i][j][k] +
														u[i][j-1][k][4]*rho_i[i][j-1][k]) -
														ty2 * ((c1*u[i][j+1][k][4] -
																c2*square[i][j+1][k]) * vp1 -
																(c1*u[i][j-1][k][4] -
																		c2*square[i][j-1][k]) * vm1);
			}
		}
	}

	/*--------------------------------------------------------------------
c     add fourth order eta-direction dissipation         
c-------------------------------------------------------------------*/
	j = 1;
#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m]- dssp *
						( five*u[i][j][k][m] - four*u[i][j+1][k][m] +
								u[i][j+2][k][m]);
			}
		}
	}

	j = 2;
#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						(-four*u[i][j-1][k][m] + six*u[i][j][k][m] -
								four*u[i][j+1][k][m] + u[i][j+2][k][m]);
			}
		}
	}

#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 3; j < grid_points[1]-3; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				for (m = 0; m < 5; m++) {
					rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
							(  u[i][j-2][k][m] - four*u[i][j-1][k][m] +
									six*u[i][j][k][m] - four*u[i][j+1][k][m] +
									u[i][j+2][k][m] );
				}
			}
		}
	}

	j = grid_points[1]-3;
#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						( u[i][j-2][k][m] - four*u[i][j-1][k][m] +
								six*u[i][j][k][m] - four*u[i][j+1][k][m] );
			}
		}
	}

	j = grid_points[1]-2;
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						( u[i][j-2][k][m] - four*u[i][j-1][k][m] +
								five*u[i][j][k][m] );
			}
		}
	}

	/*--------------------------------------------------------------------
c     compute zeta-direction fluxes 
c-------------------------------------------------------------------*/
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				wijk = ws[i][j][k];
				wp1  = ws[i][j][k+1];
				wm1  = ws[i][j][k-1];

				rhs[i][j][k][0] = rhs[i][j][k][0] + dz1tz1 *
						(u[i][j][k+1][0] - two*u[i][j][k][0] +
								u[i][j][k-1][0]) -
								tz2 * (u[i][j][k+1][3] - u[i][j][k-1][3]);
				rhs[i][j][k][1] = rhs[i][j][k][1] + dz2tz1 *
						(u[i][j][k+1][1] - two*u[i][j][k][1] +
								u[i][j][k-1][1]) +
								zzcon2 * (us[i][j][k+1] - two*us[i][j][k] +
										us[i][j][k-1]) -
										tz2 * (u[i][j][k+1][1]*wp1 -
												u[i][j][k-1][1]*wm1);
				rhs[i][j][k][2] = rhs[i][j][k][2] + dz3tz1 *
						(u[i][j][k+1][2] - two*u[i][j][k][2] +
								u[i][j][k-1][2]) +
								zzcon2 * (vs[i][j][k+1] - two*vs[i][j][k] +
										vs[i][j][k-1]) -
										tz2 * (u[i][j][k+1][2]*wp1 -
												u[i][j][k-1][2]*wm1);
				rhs[i][j][k][3] = rhs[i][j][k][3] + dz4tz1 *
						(u[i][j][k+1][3] - two*u[i][j][k][3] +
								u[i][j][k-1][3]) +
								zzcon2*con43 * (wp1 - two*wijk + wm1) -
								tz2 * (u[i][j][k+1][3]*wp1 -
										u[i][j][k-1][3]*wm1 +
										(u[i][j][k+1][4] - square[i][j][k+1] -
												u[i][j][k-1][4] + square[i][j][k-1])
												*c2);
				rhs[i][j][k][4] = rhs[i][j][k][4] + dz5tz1 *
						(u[i][j][k+1][4] - two*u[i][j][k][4] +
								u[i][j][k-1][4]) +
								zzcon3 * (qs[i][j][k+1] - two*qs[i][j][k] +
										qs[i][j][k-1]) +
										zzcon4 * (wp1*wp1 - two*wijk*wijk +
												wm1*wm1) +
												zzcon5 * (u[i][j][k+1][4]*rho_i[i][j][k+1] -
														two*u[i][j][k][4]*rho_i[i][j][k] +
														u[i][j][k-1][4]*rho_i[i][j][k-1]) -
														tz2 * ( (c1*u[i][j][k+1][4] -
																c2*square[i][j][k+1])*wp1 -
																(c1*u[i][j][k-1][4] -
																		c2*square[i][j][k-1])*wm1);
			}
		}
	}

	/*--------------------------------------------------------------------
c     add fourth order zeta-direction dissipation                
c-------------------------------------------------------------------*/
	k = 1;
#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m]- dssp *
						( five*u[i][j][k][m] - four*u[i][j][k+1][m] +
								u[i][j][k+2][m]);
			}
		}
	}

	k = 2;
#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						(-four*u[i][j][k-1][m] + six*u[i][j][k][m] -
								four*u[i][j][k+1][m] + u[i][j][k+2][m]);
			}
		}
	}

#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 3; k < grid_points[2]-3; k++) {
				for (m = 0; m < 5; m++) {
					rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
							(  u[i][j][k-2][m] - four*u[i][j][k-1][m] +
									six*u[i][j][k][m] - four*u[i][j][k+1][m] +
									u[i][j][k+2][m] );
				}
			}
		}
	}

	k = grid_points[2]-3;
#pragma omp for nowait
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						( u[i][j][k-2][m] - four*u[i][j][k-1][m] +
								six*u[i][j][k][m] - four*u[i][j][k+1][m] );
			}
		}
	}

	k = grid_points[2]-2;
#pragma omp for
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (m = 0; m < 5; m++) {
				rhs[i][j][k][m] = rhs[i][j][k][m] - dssp *
						( u[i][j][k-2][m] - four*u[i][j][k-1][m] +
								five*u[i][j][k][m] );
			}
		}
	}

#pragma omp for
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {
			for (m = 0; m < 5; m++) {
				for (i = 1; i < grid_points[0]-1; i++) {
					rhs[i][j][k][m] = rhs[i][j][k][m] * dt;
				}
			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void set_constants(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	ce[0][0]  = two;
	ce[0][1]  = zero;
	ce[0][2]  = zero;
	ce[0][3]  = four;
	ce[0][4]  = five;
	ce[0][5]  = three;
	ce[0][6]  = one/two;
	ce[0][7]  = two/hundred;
	ce[0][8]  = one/hundred;
	ce[0][9]  = three/hundred;
	ce[0][10] = one/two;
	ce[0][11] = four/ten;
	ce[0][12] = three/ten;

	ce[1][0]  = one;
	ce[1][1]  = zero;
	ce[1][2]  = zero;
	ce[1][3]  = zero;
	ce[1][4]  = one;
	ce[1][5]  = two;
	ce[1][6]  = three;
	ce[1][7]  = one/hundred;
	ce[1][8]  = three/hundred;
	ce[1][9]  = two/hundred;
	ce[1][10] = four/ten;
	ce[1][11] = three/ten;
	ce[1][12] = one/two;

	ce[2][0]  = two;
	ce[2][1]  = two;
	ce[2][2]  = zero;
	ce[2][3]  = zero;
	ce[2][4]  = zero;
	ce[2][5]  = two;
	ce[2][6]  = three;
	ce[2][7]  = four/hundred;
	ce[2][8]  = three/hundred;
	ce[2][9]  = five/hundred;
	ce[2][10] = three/ten;
	ce[2][11] = one/two;
	ce[2][12] = four/ten;

	ce[3][0]  = two;
	ce[3][1]  = two;
	ce[3][2]  = zero;
	ce[3][3]  = zero;
	ce[3][4]  = zero;
	ce[3][5]  = two;
	ce[3][6]  = three;
	ce[3][7]  = three/hundred;
	ce[3][8]  = five/hundred;
	ce[3][9]  = four/hundred;
	ce[3][10] = two/ten;
	ce[3][11] = one/ten;
	ce[3][12] = three/ten;

	ce[4][0]  = five;
	ce[4][1]  = four;
	ce[4][2]  = three;
	ce[4][3]  = two;
	ce[4][4]  = one/ten;
	ce[4][5]  = four/ten;
	ce[4][6]  = three/ten;
	ce[4][7]  = five/hundred;
	ce[4][8]  = four/hundred;
	ce[4][9]  = three/hundred;
	ce[4][10] = one/ten;
	ce[4][11] = three/ten;
	ce[4][12] = two/ten;

	c1 = one + four/ten;
	c2 = four/ten;
	c3 = one/ten;
	c4 = one;
	c5 = one + four/ten;

	dnxm1 = one / (element_t)(grid_points[0]-1);
	dnym1 = one / (element_t)(grid_points[1]-1);
	dnzm1 = one / (element_t)(grid_points[2]-1);

	c1c2 = c1 * c2;
	c1c5 = c1 * c5;
	c3c4 = c3 * c4;
	c1345 = c1c5 * c3c4;

	conz1 = (one-c1c5);

	tx1 = one / (dnxm1 * dnxm1);
	tx2 = one / (two * dnxm1);
	tx3 = one / dnxm1;

	ty1 = one / (dnym1 * dnym1);
	ty2 = one / (two * dnym1);
	ty3 = one / dnym1;

	tz1 = one / (dnzm1 * dnzm1);
	tz2 = one / (two * dnzm1);
	tz3 = one / dnzm1;

	dx1 = three/four;
	dx2 = three/four;
	dx3 = three/four;
	dx4 = three/four;
	dx5 = three/four;

	dy1 = three/four;
	dy2 = three/four;
	dy3 = three/four;
	dy4 = three/four;
	dy5 = three/four;

	dz1 = one;
	dz2 = one;
	dz3 = one;
	dz4 = one;
	dz5 = one;

	dxmax = max(dx3, dx4);
	dymax = max(dy2, dy4);
	dzmax = max(dz2, dz3);

	dssp = one/four * max(dx1, max(dy1, dz1) );

	c4dssp = four * dssp;
	c5dssp = five * dssp;

	dttx1 = dt*tx1;
	dttx2 = dt*tx2;
	dtty1 = dt*ty1;
	dtty2 = dt*ty2;
	dttz1 = dt*tz1;
	dttz2 = dt*tz2;

	c2dttx1 = two*dttx1;
	c2dtty1 = two*dtty1;
	c2dttz1 = two*dttz1;

	dtdssp = dt*dssp;

	comz1  = dtdssp;
	comz4  = four*dtdssp;
	comz5  = five*dtdssp;
	comz6  = six*dtdssp;

	c3c4tx3 = c3c4*tx3;
	c3c4ty3 = c3c4*ty3;
	c3c4tz3 = c3c4*tz3;

	dx1tx1 = dx1*tx1;
	dx2tx1 = dx2*tx1;
	dx3tx1 = dx3*tx1;
	dx4tx1 = dx4*tx1;
	dx5tx1 = dx5*tx1;

	dy1ty1 = dy1*ty1;
	dy2ty1 = dy2*ty1;
	dy3ty1 = dy3*ty1;
	dy4ty1 = dy4*ty1;
	dy5ty1 = dy5*ty1;

	dz1tz1 = dz1*tz1;
	dz2tz1 = dz2*tz1;
	dz3tz1 = dz3*tz1;
	dz4tz1 = dz4*tz1;
	dz5tz1 = dz5*tz1;

	c2iv  = five/two;
	con43 = four/three;
	con16 = one/six;

	xxcon1 = c3c4tx3*con43*tx3;
	xxcon2 = c3c4tx3*tx3;
	xxcon3 = c3c4tx3*conz1*tx3;
	xxcon4 = c3c4tx3*con16*tx3;
	xxcon5 = c3c4tx3*c1c5*tx3;

	yycon1 = c3c4ty3*con43*ty3;
	yycon2 = c3c4ty3*ty3;
	yycon3 = c3c4ty3*conz1*ty3;
	yycon4 = c3c4ty3*con16*ty3;
	yycon5 = c3c4ty3*c1c5*ty3;

	zzcon1 = c3c4tz3*con43*tz3;
	zzcon2 = c3c4tz3*tz3;
	zzcon3 = c3c4tz3*conz1*tz3;
	zzcon4 = c3c4tz3*con16*tz3;
	zzcon5 = c3c4tz3*c1c5*tz3;
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void verify(int no_time_steps, char *class, boolean *verified) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c  verification routine                         
c-------------------------------------------------------------------*/

	element_t xcrref[5],xceref[5],xcrdif[5],xcedif[5],
	epsilon, xce[5], xcr[5], dtref;
	int m;

	/*--------------------------------------------------------------------
c   tolerance level
c-------------------------------------------------------------------*/
	epsilon = c_epsilon;


	/*--------------------------------------------------------------------
c   compute the error norm and the residual norm, and exit if not printing
c-------------------------------------------------------------------*/
	error_norm(xce);
	compute_rhs();

	rhs_norm(xcr);

	for (m = 0; m < 5; m++) {
		xcr[m] = xcr[m] / dt;
	}

	*class = 'U';
	*verified = TRUE;

	for (m = 0; m < 5; m++) {
		xcrref[m] = one;
		xceref[m] = one;
	}

	/*--------------------------------------------------------------------
c    reference data for 12X12X12 grids after 100 time steps, with DT = 1.0d-02
c-------------------------------------------------------------------*/
	if (grid_points[0] == 12 &&
			grid_points[1] == 12 &&
			grid_points[2] == 12 &&
			no_time_steps == 60) {

		*class = 'S';
		dtref = c_dtref;

		/*--------------------------------------------------------------------
c  Reference values of RMS-norms of residual.
c-------------------------------------------------------------------*/
#if defined WITH_POSIT_32
		*((uint32_t*)&xcrref[0]) = 0x357372d2;
		*((uint32_t*)&xcrref[1]) = 0x26a4b136;
		*((uint32_t*)&xcrref[2]) = 0x2c29e007;
		*((uint32_t*)&xcrref[3]) = 0x2ac4898c;
		*((uint32_t*)&xcrref[4]) = 0x3625d450;
#elif (defined WITH_POSIT_16)
		*((uint32_t*)&xcrref[0]) = 0x2ae6;
		*((uint32_t*)&xcrref[1]) = 0x16a4;
		*((uint32_t*)&xcrref[2]) = 0x1c29;
		*((uint32_t*)&xcrref[3]) = 0x1ac4;
		*((uint32_t*)&xcrref[4]) = 0x2c4b;
#elif (defined WITH_POSIT_8)
		*((uint32_t*)&xcrref[0]) = 0x1a;
		*((uint32_t*)&xcrref[1]) = 0x7;
		*((uint32_t*)&xcrref[2]) = 0xc;
		*((uint32_t*)&xcrref[3]) = 0xa;
		*((uint32_t*)&xcrref[4]) = 0x1c;
#else
		xcrref[0] = 1.7034283709541311e-01;
		xcrref[1] = 1.2975252070034097e-02;
		xcrref[2] = 3.2527926989486055e-02;
		xcrref[3] = 2.6436421275166801e-02;
		xcrref[4] = 1.9211784131744430e-01;
#endif

		/*--------------------------------------------------------------------
c  Reference values of RMS-norms of solution error.
c-------------------------------------------------------------------*/
#if defined WITH_POSIT_32
		*((uint32_t*)&xceref[0]) = 0x1a0c0bc1;
		*((uint32_t*)&xceref[1]) = 0x12f641e9;
		*((uint32_t*)&xceref[2]) = 0x146c8973;
		*((uint32_t*)&xceref[3]) = 0x146b41e7;
		*((uint32_t*)&xceref[4]) = 0x1ba80f57;
#elif (defined WITH_POSIT_16)
		*((uint32_t*)&xceref[0]) = 0xa0c;
		*((uint32_t*)&xceref[1]) = 0x57b;
		*((uint32_t*)&xceref[2]) = 0x636;
		*((uint32_t*)&xceref[3]) = 0x635;
		*((uint32_t*)&xceref[4]) = 0xba8;
#elif (defined WITH_POSIT_8)
		*((uint32_t*)&xceref[0]) = 0x1;
		*((uint32_t*)&xceref[1]) = 0x1;
		*((uint32_t*)&xceref[2]) = 0x1;
		*((uint32_t*)&xceref[3]) = 0x1;
		*((uint32_t*)&xceref[4]) = 0x1;
#else
		xceref[0] = 4.9976913345811579e-04;
		xceref[1] = 4.5195666782961927e-05;
		xceref[2] = 7.3973765172921357e-05;
		xceref[3] = 7.3821238632439731e-05;
		xceref[4] = 8.9269630987491446e-04;
#endif
	}
	/*
     TODO - for now we use only class S
//--------------------------------------------------------------------
//    reference data for 24X24X24 grids after 200 time steps, with DT = 0.8d-3
//-------------------------------------------------------------------//
    } else if (grid_points[0] == 24 &&
	       grid_points[1] == 24 &&
	       grid_points[2] == 24 &&
	       no_time_steps == 200) {

	 *class = 'W';
      dtref = 0.8e-3;
//--------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//-------------------------------------------------------------------//
      xcrref[0] = 0.1125590409344e+03;
      xcrref[1] = 0.1180007595731e+02;
      xcrref[2] = 0.2710329767846e+02;
      xcrref[3] = 0.2469174937669e+02;
      xcrref[4] = 0.2638427874317e+03;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//-------------------------------------------------------------------//
      xceref[0] = 0.4419655736008e+01;
      xceref[1] = 0.4638531260002e+00;
      xceref[2] = 0.1011551749967e+01;
      xceref[3] = 0.9235878729944e+00;
      xceref[4] = 0.1018045837718e+02;


//--------------------------------------------------------------------
//    reference data for 64X64X64 grids after 200 time steps, with DT = 0.8d-3
//-------------------------------------------------------------------//
    } else if (grid_points[0] == 64 &&
	       grid_points[1] == 64 &&
	       grid_points[2] == 64 &&
	       no_time_steps == 200) {

	 *class = 'A';
      dtref = 0.8e-3;
//--------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//-------------------------------------------------------------------//
      xcrref[0] = 1.0806346714637264e+02;
      xcrref[1] = 1.1319730901220813e+01;
      xcrref[2] = 2.5974354511582465e+01;
      xcrref[3] = 2.3665622544678910e+01;
      xcrref[4] = 2.5278963211748344e+02;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//-------------------------------------------------------------------//
      xceref[0] = 4.2348416040525025e+00;
      xceref[1] = 4.4390282496995698e-01;
      xceref[2] = 9.6692480136345650e-01;
      xceref[3] = 8.8302063039765474e-01;
      xceref[4] = 9.7379901770829278e+00;

//--------------------------------------------------------------------
//    reference data for 102X102X102 grids after 200 time steps,
//    with DT = 3.0d-04
//-------------------------------------------------------------------//
    } else if (grid_points[0] == 102 &&
	       grid_points[1] == 102 &&
	       grid_points[2] == 102 &&
	       no_time_steps == 200) {

	 *class = 'B';
      dtref = 3.0e-4;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//-------------------------------------------------------------------//
      xcrref[0] = 1.4233597229287254e+03;
      xcrref[1] = 9.9330522590150238e+01;
      xcrref[2] = 3.5646025644535285e+02;
      xcrref[3] = 3.2485447959084092e+02;
      xcrref[4] = 3.2707541254659363e+03;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//-------------------------------------------------------------------//
      xceref[0] = 5.2969847140936856e+01;
      xceref[1] = 4.4632896115670668e+00;
      xceref[2] = 1.3122573342210174e+01;
      xceref[3] = 1.2006925323559144e+01;
      xceref[4] = 1.2459576151035986e+02;

//--------------------------------------------------------------------
//    reference data for 162X162X162 grids after 200 time steps,
//    with DT = 1.0d-04
//-------------------------------------------------------------------//
    } else if (grid_points[0] == 162 &&
	       grid_points[1] == 162 &&
	       grid_points[2] == 162 &&
	       no_time_steps == 200) {

	 *class = 'C';
      dtref = 1.0e-4;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of residual.
//-------------------------------------------------------------------//
      xcrref[0] = 0.62398116551764615e+04;
      xcrref[1] = 0.50793239190423964e+03;
      xcrref[2] = 0.15423530093013596e+04;
      xcrref[3] = 0.13302387929291190e+04;
      xcrref[4] = 0.11604087428436455e+05;

//--------------------------------------------------------------------
//  Reference values of RMS-norms of solution error.
//-------------------------------------------------------------------//
      xceref[0] = 0.16462008369091265e+03;
      xceref[1] = 0.11497107903824313e+02;
      xceref[2] = 0.41207446207461508e+02;
      xceref[3] = 0.37087651059694167e+02;
      xceref[4] = 0.36211053051841265e+03;

    } else {
	 *verified = FALSE;
    }
	 */
	//--------------------------------------------------------------------
	//    verification test for residuals if gridsize is either 12X12X12 or
	//    64X64X64 or 102X102X102 or 162X162X162
	//-------------------------------------------------------------------//

	//--------------------------------------------------------------------
	//    Compute the difference of solution values and the known reference values.
	//-------------------------------------------------------------------//
	for (m = 0; m < 5; m++) {

		xcrdif[m] = my_fabs((xcr[m]-xcrref[m])/xcrref[m]);
		xcedif[m] = my_fabs((xce[m]-xceref[m])/xceref[m]);

	}

	//--------------------------------------------------------------------
	//    Output the comparison of computed results to known cases.
	//-------------------------------------------------------------------//

	if (*class != 'U') {
#ifdef PFDEBUG
		printf(" Verification being performed for class %1c\n", *class);
		printf(" accuracy setting for epsilon = %20.13e\n", epsilon);
#endif
		if (my_fabs(dt-dtref) > epsilon) {
			*verified = FALSE;
			*class = 'U';
#ifdef PFDEBUG
			printf(" DT does not match the reference value of %15.8e\n", dtref);
#endif
		}
	} else {
#ifdef PFDEBUG
		printf(" Unknown class\n");
#endif
	}

#ifdef PFDEBUG
	if (*class != 'U') {
		printf(" Comparison of RMS-norms of residual\n");
	} else {
		printf(" RMS-norms of residual\n");
	}
#endif

	for (m = 0; m < 5; m++) {
		if (*class == 'U') {
#ifdef PFDEBUG
			printf("          %2d%20.13e\n", m, xcr[m]);
#endif
		} else if (xcrdif[m] > epsilon) {
			*verified = FALSE;
#ifdef PFDEBUG
			printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n",m, xcr[m], xcrref[m], xcrdif[m]);
#endif
		} else {
#ifdef PFDEBUG
			printf("          %2d%20.13e%20.13e%20.13e\n",m, xcr[m], xcrref[m], xcrdif[m]);
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
		if (*class == 'U') {
#ifdef PFDEBUG
			printf("          %2d%20.13e\n", m, xce[m]);
#endif
		} else if (xcedif[m] > epsilon) {
			*verified = FALSE;
#ifdef PFDEBUG
			printf(" FAILURE: %2d%20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
#endif
		} else {
#ifdef PFDEBUG
			printf("          %2d%20.13e%20.13e%20.13e\n", m, xce[m], xceref[m], xcedif[m]);
#endif
		}
	}

#ifdef PFDEBUG
	if (*class == 'U') {
		printf(" No reference values provided\n");
		printf(" No verification performed\n");
	} else if (*verified == TRUE) {
		printf(" Verification Successful\n");
	} else {
		printf(" Verification failed\n");
	}
#endif
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void x_solve(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     
c     Performs line solves in X direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c     
c-------------------------------------------------------------------*/

	lhsx();
	x_solve_cell();
	x_backsubstitute();
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void x_backsubstitute(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     back solve: if last cell, then generate U(isize)=rhs[isize)
c     else assume U(isize) is loaded in un pack backsub_info
c     so just use it
c     after call u(istart) will be sent to next cell
c-------------------------------------------------------------------*/

	int i, j, k, m, n;

	for (i = grid_points[0]-2; i >= 0; i--) {
#pragma omp for
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				for (m = 0; m < BLOCK_SIZE; m++) {
					for (n = 0; n < BLOCK_SIZE; n++) {
						rhs[i][j][k][m] = rhs[i][j][k][m]
													   - lhs[i][j][k][CC][m][n]*rhs[i+1][j][k][n];
					}
				}
			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void x_solve_cell(void) {

	/*--------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(IMAX) and rhs'(IMAX) will be sent to next cell
c-------------------------------------------------------------------*/

	int i,j,k,isize;

	isize = grid_points[0]-1;

	/*--------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c-------------------------------------------------------------------*/
#pragma omp for
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {

			/*--------------------------------------------------------------------
c     multiply c(0,j,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c-------------------------------------------------------------------*/
			binvcrhs( lhs[0][j][k][BB],
					lhs[0][j][k][CC],
					rhs[0][j][k] );
		}
	}

	/*--------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c-------------------------------------------------------------------*/
	for (i = 1; i < isize; i++) {
#pragma omp for
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = 1; k < grid_points[2]-1; k++) {

				/*--------------------------------------------------------------------
c     rhs(i) = rhs(i) - A*rhs(i-1)
c-------------------------------------------------------------------*/
				matvec_sub(lhs[i][j][k][AA],
						rhs[i-1][j][k], rhs[i][j][k]);

				/*--------------------------------------------------------------------
c     B(i) = B(i) - C(i-1)*A(i)
c-------------------------------------------------------------------*/
				matmul_sub(lhs[i][j][k][AA],
						lhs[i-1][j][k][CC],
						lhs[i][j][k][BB]);


				/*--------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(1,j,k) by b_inverse(1,j,k) and copy to rhs
c-------------------------------------------------------------------*/
				binvcrhs( lhs[i][j][k][BB],
						lhs[i][j][k][CC],
						rhs[i][j][k] );

			}
		}
	}

#pragma omp for
	for (j = 1; j < grid_points[1]-1; j++) {
		for (k = 1; k < grid_points[2]-1; k++) {

			/*--------------------------------------------------------------------
c     rhs(isize) = rhs(isize) - A*rhs(isize-1)
c-------------------------------------------------------------------*/
			matvec_sub(lhs[isize][j][k][AA],
					rhs[isize-1][j][k], rhs[isize][j][k]);

			/*--------------------------------------------------------------------
c     B(isize) = B(isize) - C(isize-1)*A(isize)
c-------------------------------------------------------------------*/
			matmul_sub(lhs[isize][j][k][AA],
					lhs[isize-1][j][k][CC],
					lhs[isize][j][k][BB]);

			/*--------------------------------------------------------------------
c     multiply rhs() by b_inverse() and copy to rhs
c-------------------------------------------------------------------*/
			binvrhs( lhs[i][j][k][BB],
					rhs[i][j][k] );

		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void matvec_sub(element_t ablock[5][5], element_t avec[5], element_t bvec[5]) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     subtracts bvec=bvec - ablock*avec
c-------------------------------------------------------------------*/

	int i;

	for (i = 0; i < 5; i++) {
		/*--------------------------------------------------------------------
c            rhs(i,ic,jc,kc,ccell) = rhs(i,ic,jc,kc,ccell) 
c     $           - lhs[i,1,ablock,ia,ja,ka,acell)*
c-------------------------------------------------------------------*/
		bvec[i] = bvec[i] - ablock[i][0]*avec[0]
											  - ablock[i][1]*avec[1]
																  - ablock[i][2]*avec[2]
																					  - ablock[i][3]*avec[3]
																										  - ablock[i][4]*avec[4];
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void matmul_sub(element_t ablock[5][5], element_t bblock[5][5],
		element_t cblock[5][5]) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
c-------------------------------------------------------------------*/

	int j;

	for (j = 0; j < 5; j++) {
		cblock[0][j] = cblock[0][j] - ablock[0][0]*bblock[0][j]
															 - ablock[0][1]*bblock[1][j]
																					  - ablock[0][2]*bblock[2][j]
																											   - ablock[0][3]*bblock[3][j]
																																		- ablock[0][4]*bblock[4][j];
		cblock[1][j] = cblock[1][j] - ablock[1][0]*bblock[0][j]
															 - ablock[1][1]*bblock[1][j]
																					  - ablock[1][2]*bblock[2][j]
																											   - ablock[1][3]*bblock[3][j]
																																		- ablock[1][4]*bblock[4][j];
		cblock[2][j] = cblock[2][j] - ablock[2][0]*bblock[0][j]
															 - ablock[2][1]*bblock[1][j]
																					  - ablock[2][2]*bblock[2][j]
																											   - ablock[2][3]*bblock[3][j]
																																		- ablock[2][4]*bblock[4][j];
		cblock[3][j] = cblock[3][j] - ablock[3][0]*bblock[0][j]
															 - ablock[3][1]*bblock[1][j]
																					  - ablock[3][2]*bblock[2][j]
																											   - ablock[3][3]*bblock[3][j]
																																		- ablock[3][4]*bblock[4][j];
		cblock[4][j] = cblock[4][j] - ablock[4][0]*bblock[0][j]
															 - ablock[4][1]*bblock[1][j]
																					  - ablock[4][2]*bblock[2][j]
																											   - ablock[4][3]*bblock[3][j]
																																		- ablock[4][4]*bblock[4][j];
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void binvcrhs(element_t lhs[5][5], element_t c[5][5], element_t r[5]) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	element_t pivot, coeff;

	/*--------------------------------------------------------------------
c     
c-------------------------------------------------------------------*/

	pivot = one/lhs[0][0];
	lhs[0][1] = lhs[0][1]*pivot;
	lhs[0][2] = lhs[0][2]*pivot;
	lhs[0][3] = lhs[0][3]*pivot;
	lhs[0][4] = lhs[0][4]*pivot;
	c[0][0] = c[0][0]*pivot;
	c[0][1] = c[0][1]*pivot;
	c[0][2] = c[0][2]*pivot;
	c[0][3] = c[0][3]*pivot;
	c[0][4] = c[0][4]*pivot;
	r[0]   = r[0]  *pivot;

	coeff = lhs[1][0];
	lhs[1][1]= lhs[1][1] - coeff*lhs[0][1];
	lhs[1][2]= lhs[1][2] - coeff*lhs[0][2];
	lhs[1][3]= lhs[1][3] - coeff*lhs[0][3];
	lhs[1][4]= lhs[1][4] - coeff*lhs[0][4];
	c[1][0] = c[1][0] - coeff*c[0][0];
	c[1][1] = c[1][1] - coeff*c[0][1];
	c[1][2] = c[1][2] - coeff*c[0][2];
	c[1][3] = c[1][3] - coeff*c[0][3];
	c[1][4] = c[1][4] - coeff*c[0][4];
	r[1]   = r[1]   - coeff*r[0];

	coeff = lhs[2][0];
	lhs[2][1]= lhs[2][1] - coeff*lhs[0][1];
	lhs[2][2]= lhs[2][2] - coeff*lhs[0][2];
	lhs[2][3]= lhs[2][3] - coeff*lhs[0][3];
	lhs[2][4]= lhs[2][4] - coeff*lhs[0][4];
	c[2][0] = c[2][0] - coeff*c[0][0];
	c[2][1] = c[2][1] - coeff*c[0][1];
	c[2][2] = c[2][2] - coeff*c[0][2];
	c[2][3] = c[2][3] - coeff*c[0][3];
	c[2][4] = c[2][4] - coeff*c[0][4];
	r[2]   = r[2]   - coeff*r[0];

	coeff = lhs[3][0];
	lhs[3][1]= lhs[3][1] - coeff*lhs[0][1];
	lhs[3][2]= lhs[3][2] - coeff*lhs[0][2];
	lhs[3][3]= lhs[3][3] - coeff*lhs[0][3];
	lhs[3][4]= lhs[3][4] - coeff*lhs[0][4];
	c[3][0] = c[3][0] - coeff*c[0][0];
	c[3][1] = c[3][1] - coeff*c[0][1];
	c[3][2] = c[3][2] - coeff*c[0][2];
	c[3][3] = c[3][3] - coeff*c[0][3];
	c[3][4] = c[3][4] - coeff*c[0][4];
	r[3]   = r[3]   - coeff*r[0];

	coeff = lhs[4][0];
	lhs[4][1]= lhs[4][1] - coeff*lhs[0][1];
	lhs[4][2]= lhs[4][2] - coeff*lhs[0][2];
	lhs[4][3]= lhs[4][3] - coeff*lhs[0][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[0][4];
	c[4][0] = c[4][0] - coeff*c[0][0];
	c[4][1] = c[4][1] - coeff*c[0][1];
	c[4][2] = c[4][2] - coeff*c[0][2];
	c[4][3] = c[4][3] - coeff*c[0][3];
	c[4][4] = c[4][4] - coeff*c[0][4];
	r[4]   = r[4]   - coeff*r[0];


	pivot = one/lhs[1][1];
	lhs[1][2] = lhs[1][2]*pivot;
	lhs[1][3] = lhs[1][3]*pivot;
	lhs[1][4] = lhs[1][4]*pivot;
	c[1][0] = c[1][0]*pivot;
	c[1][1] = c[1][1]*pivot;
	c[1][2] = c[1][2]*pivot;
	c[1][3] = c[1][3]*pivot;
	c[1][4] = c[1][4]*pivot;
	r[1]   = r[1]  *pivot;

	coeff = lhs[0][1];
	lhs[0][2]= lhs[0][2] - coeff*lhs[1][2];
	lhs[0][3]= lhs[0][3] - coeff*lhs[1][3];
	lhs[0][4]= lhs[0][4] - coeff*lhs[1][4];
	c[0][0] = c[0][0] - coeff*c[1][0];
	c[0][1] = c[0][1] - coeff*c[1][1];
	c[0][2] = c[0][2] - coeff*c[1][2];
	c[0][3] = c[0][3] - coeff*c[1][3];
	c[0][4] = c[0][4] - coeff*c[1][4];
	r[0]   = r[0]   - coeff*r[1];

	coeff = lhs[2][1];
	lhs[2][2]= lhs[2][2] - coeff*lhs[1][2];
	lhs[2][3]= lhs[2][3] - coeff*lhs[1][3];
	lhs[2][4]= lhs[2][4] - coeff*lhs[1][4];
	c[2][0] = c[2][0] - coeff*c[1][0];
	c[2][1] = c[2][1] - coeff*c[1][1];
	c[2][2] = c[2][2] - coeff*c[1][2];
	c[2][3] = c[2][3] - coeff*c[1][3];
	c[2][4] = c[2][4] - coeff*c[1][4];
	r[2]   = r[2]   - coeff*r[1];

	coeff = lhs[3][1];
	lhs[3][2]= lhs[3][2] - coeff*lhs[1][2];
	lhs[3][3]= lhs[3][3] - coeff*lhs[1][3];
	lhs[3][4]= lhs[3][4] - coeff*lhs[1][4];
	c[3][0] = c[3][0] - coeff*c[1][0];
	c[3][1] = c[3][1] - coeff*c[1][1];
	c[3][2] = c[3][2] - coeff*c[1][2];
	c[3][3] = c[3][3] - coeff*c[1][3];
	c[3][4] = c[3][4] - coeff*c[1][4];
	r[3]   = r[3]   - coeff*r[1];

	coeff = lhs[4][1];
	lhs[4][2]= lhs[4][2] - coeff*lhs[1][2];
	lhs[4][3]= lhs[4][3] - coeff*lhs[1][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[1][4];
	c[4][0] = c[4][0] - coeff*c[1][0];
	c[4][1] = c[4][1] - coeff*c[1][1];
	c[4][2] = c[4][2] - coeff*c[1][2];
	c[4][3] = c[4][3] - coeff*c[1][3];
	c[4][4] = c[4][4] - coeff*c[1][4];
	r[4]   = r[4]   - coeff*r[1];


	pivot = one/lhs[2][2];
	lhs[2][3] = lhs[2][3]*pivot;
	lhs[2][4] = lhs[2][4]*pivot;
	c[2][0] = c[2][0]*pivot;
	c[2][1] = c[2][1]*pivot;
	c[2][2] = c[2][2]*pivot;
	c[2][3] = c[2][3]*pivot;
	c[2][4] = c[2][4]*pivot;
	r[2]   = r[2]  *pivot;

	coeff = lhs[0][2];
	lhs[0][3]= lhs[0][3] - coeff*lhs[2][3];
	lhs[0][4]= lhs[0][4] - coeff*lhs[2][4];
	c[0][0] = c[0][0] - coeff*c[2][0];
	c[0][1] = c[0][1] - coeff*c[2][1];
	c[0][2] = c[0][2] - coeff*c[2][2];
	c[0][3] = c[0][3] - coeff*c[2][3];
	c[0][4] = c[0][4] - coeff*c[2][4];
	r[0]   = r[0]   - coeff*r[2];

	coeff = lhs[1][2];
	lhs[1][3]= lhs[1][3] - coeff*lhs[2][3];
	lhs[1][4]= lhs[1][4] - coeff*lhs[2][4];
	c[1][0] = c[1][0] - coeff*c[2][0];
	c[1][1] = c[1][1] - coeff*c[2][1];
	c[1][2] = c[1][2] - coeff*c[2][2];
	c[1][3] = c[1][3] - coeff*c[2][3];
	c[1][4] = c[1][4] - coeff*c[2][4];
	r[1]   = r[1]   - coeff*r[2];

	coeff = lhs[3][2];
	lhs[3][3]= lhs[3][3] - coeff*lhs[2][3];
	lhs[3][4]= lhs[3][4] - coeff*lhs[2][4];
	c[3][0] = c[3][0] - coeff*c[2][0];
	c[3][1] = c[3][1] - coeff*c[2][1];
	c[3][2] = c[3][2] - coeff*c[2][2];
	c[3][3] = c[3][3] - coeff*c[2][3];
	c[3][4] = c[3][4] - coeff*c[2][4];
	r[3]   = r[3]   - coeff*r[2];

	coeff = lhs[4][2];
	lhs[4][3]= lhs[4][3] - coeff*lhs[2][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[2][4];
	c[4][0] = c[4][0] - coeff*c[2][0];
	c[4][1] = c[4][1] - coeff*c[2][1];
	c[4][2] = c[4][2] - coeff*c[2][2];
	c[4][3] = c[4][3] - coeff*c[2][3];
	c[4][4] = c[4][4] - coeff*c[2][4];
	r[4]   = r[4]   - coeff*r[2];


	pivot = one/lhs[3][3];
	lhs[3][4] = lhs[3][4]*pivot;
	c[3][0] = c[3][0]*pivot;
	c[3][1] = c[3][1]*pivot;
	c[3][2] = c[3][2]*pivot;
	c[3][3] = c[3][3]*pivot;
	c[3][4] = c[3][4]*pivot;
	r[3]   = r[3]  *pivot;

	coeff = lhs[0][3];
	lhs[0][4]= lhs[0][4] - coeff*lhs[3][4];
	c[0][0] = c[0][0] - coeff*c[3][0];
	c[0][1] = c[0][1] - coeff*c[3][1];
	c[0][2] = c[0][2] - coeff*c[3][2];
	c[0][3] = c[0][3] - coeff*c[3][3];
	c[0][4] = c[0][4] - coeff*c[3][4];
	r[0]   = r[0]   - coeff*r[3];

	coeff = lhs[1][3];
	lhs[1][4]= lhs[1][4] - coeff*lhs[3][4];
	c[1][0] = c[1][0] - coeff*c[3][0];
	c[1][1] = c[1][1] - coeff*c[3][1];
	c[1][2] = c[1][2] - coeff*c[3][2];
	c[1][3] = c[1][3] - coeff*c[3][3];
	c[1][4] = c[1][4] - coeff*c[3][4];
	r[1]   = r[1]   - coeff*r[3];

	coeff = lhs[2][3];
	lhs[2][4]= lhs[2][4] - coeff*lhs[3][4];
	c[2][0] = c[2][0] - coeff*c[3][0];
	c[2][1] = c[2][1] - coeff*c[3][1];
	c[2][2] = c[2][2] - coeff*c[3][2];
	c[2][3] = c[2][3] - coeff*c[3][3];
	c[2][4] = c[2][4] - coeff*c[3][4];
	r[2]   = r[2]   - coeff*r[3];

	coeff = lhs[4][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[3][4];
	c[4][0] = c[4][0] - coeff*c[3][0];
	c[4][1] = c[4][1] - coeff*c[3][1];
	c[4][2] = c[4][2] - coeff*c[3][2];
	c[4][3] = c[4][3] - coeff*c[3][3];
	c[4][4] = c[4][4] - coeff*c[3][4];
	r[4]   = r[4]   - coeff*r[3];


	pivot = one/lhs[4][4];
	c[4][0] = c[4][0]*pivot;
	c[4][1] = c[4][1]*pivot;
	c[4][2] = c[4][2]*pivot;
	c[4][3] = c[4][3]*pivot;
	c[4][4] = c[4][4]*pivot;
	r[4]   = r[4]  *pivot;

	coeff = lhs[0][4];
	c[0][0] = c[0][0] - coeff*c[4][0];
	c[0][1] = c[0][1] - coeff*c[4][1];
	c[0][2] = c[0][2] - coeff*c[4][2];
	c[0][3] = c[0][3] - coeff*c[4][3];
	c[0][4] = c[0][4] - coeff*c[4][4];
	r[0]   = r[0]   - coeff*r[4];

	coeff = lhs[1][4];
	c[1][0] = c[1][0] - coeff*c[4][0];
	c[1][1] = c[1][1] - coeff*c[4][1];
	c[1][2] = c[1][2] - coeff*c[4][2];
	c[1][3] = c[1][3] - coeff*c[4][3];
	c[1][4] = c[1][4] - coeff*c[4][4];
	r[1]   = r[1]   - coeff*r[4];

	coeff = lhs[2][4];
	c[2][0] = c[2][0] - coeff*c[4][0];
	c[2][1] = c[2][1] - coeff*c[4][1];
	c[2][2] = c[2][2] - coeff*c[4][2];
	c[2][3] = c[2][3] - coeff*c[4][3];
	c[2][4] = c[2][4] - coeff*c[4][4];
	r[2]   = r[2]   - coeff*r[4];

	coeff = lhs[3][4];
	c[3][0] = c[3][0] - coeff*c[4][0];
	c[3][1] = c[3][1] - coeff*c[4][1];
	c[3][2] = c[3][2] - coeff*c[4][2];
	c[3][3] = c[3][3] - coeff*c[4][3];
	c[3][4] = c[3][4] - coeff*c[4][4];
	r[3]   = r[3]   - coeff*r[4];
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void binvrhs( element_t lhs[5][5], element_t r[5] ) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	element_t pivot, coeff;

	/*--------------------------------------------------------------------
c     
c-------------------------------------------------------------------*/

	pivot = one/lhs[0][0];
	lhs[0][1] = lhs[0][1]*pivot;
	lhs[0][2] = lhs[0][2]*pivot;
	lhs[0][3] = lhs[0][3]*pivot;
	lhs[0][4] = lhs[0][4]*pivot;
	r[0]   = r[0]  *pivot;

	coeff = lhs[1][0];
	lhs[1][1]= lhs[1][1] - coeff*lhs[0][1];
	lhs[1][2]= lhs[1][2] - coeff*lhs[0][2];
	lhs[1][3]= lhs[1][3] - coeff*lhs[0][3];
	lhs[1][4]= lhs[1][4] - coeff*lhs[0][4];
	r[1]   = r[1]   - coeff*r[0];

	coeff = lhs[2][0];
	lhs[2][1]= lhs[2][1] - coeff*lhs[0][1];
	lhs[2][2]= lhs[2][2] - coeff*lhs[0][2];
	lhs[2][3]= lhs[2][3] - coeff*lhs[0][3];
	lhs[2][4]= lhs[2][4] - coeff*lhs[0][4];
	r[2]   = r[2]   - coeff*r[0];

	coeff = lhs[3][0];
	lhs[3][1]= lhs[3][1] - coeff*lhs[0][1];
	lhs[3][2]= lhs[3][2] - coeff*lhs[0][2];
	lhs[3][3]= lhs[3][3] - coeff*lhs[0][3];
	lhs[3][4]= lhs[3][4] - coeff*lhs[0][4];
	r[3]   = r[3]   - coeff*r[0];

	coeff = lhs[4][0];
	lhs[4][1]= lhs[4][1] - coeff*lhs[0][1];
	lhs[4][2]= lhs[4][2] - coeff*lhs[0][2];
	lhs[4][3]= lhs[4][3] - coeff*lhs[0][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[0][4];
	r[4]   = r[4]   - coeff*r[0];


	pivot = one/lhs[1][1];
	lhs[1][2] = lhs[1][2]*pivot;
	lhs[1][3] = lhs[1][3]*pivot;
	lhs[1][4] = lhs[1][4]*pivot;
	r[1]   = r[1]  *pivot;

	coeff = lhs[0][1];
	lhs[0][2]= lhs[0][2] - coeff*lhs[1][2];
	lhs[0][3]= lhs[0][3] - coeff*lhs[1][3];
	lhs[0][4]= lhs[0][4] - coeff*lhs[1][4];
	r[0]   = r[0]   - coeff*r[1];

	coeff = lhs[2][1];
	lhs[2][2]= lhs[2][2] - coeff*lhs[1][2];
	lhs[2][3]= lhs[2][3] - coeff*lhs[1][3];
	lhs[2][4]= lhs[2][4] - coeff*lhs[1][4];
	r[2]   = r[2]   - coeff*r[1];

	coeff = lhs[3][1];
	lhs[3][2]= lhs[3][2] - coeff*lhs[1][2];
	lhs[3][3]= lhs[3][3] - coeff*lhs[1][3];
	lhs[3][4]= lhs[3][4] - coeff*lhs[1][4];
	r[3]   = r[3]   - coeff*r[1];

	coeff = lhs[4][1];
	lhs[4][2]= lhs[4][2] - coeff*lhs[1][2];
	lhs[4][3]= lhs[4][3] - coeff*lhs[1][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[1][4];
	r[4]   = r[4]   - coeff*r[1];


	pivot = one/lhs[2][2];
	lhs[2][3] = lhs[2][3]*pivot;
	lhs[2][4] = lhs[2][4]*pivot;
	r[2]   = r[2]  *pivot;

	coeff = lhs[0][2];
	lhs[0][3]= lhs[0][3] - coeff*lhs[2][3];
	lhs[0][4]= lhs[0][4] - coeff*lhs[2][4];
	r[0]   = r[0]   - coeff*r[2];

	coeff = lhs[1][2];
	lhs[1][3]= lhs[1][3] - coeff*lhs[2][3];
	lhs[1][4]= lhs[1][4] - coeff*lhs[2][4];
	r[1]   = r[1]   - coeff*r[2];

	coeff = lhs[3][2];
	lhs[3][3]= lhs[3][3] - coeff*lhs[2][3];
	lhs[3][4]= lhs[3][4] - coeff*lhs[2][4];
	r[3]   = r[3]   - coeff*r[2];

	coeff = lhs[4][2];
	lhs[4][3]= lhs[4][3] - coeff*lhs[2][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[2][4];
	r[4]   = r[4]   - coeff*r[2];


	pivot = one/lhs[3][3];
	lhs[3][4] = lhs[3][4]*pivot;
	r[3]   = r[3]  *pivot;

	coeff = lhs[0][3];
	lhs[0][4]= lhs[0][4] - coeff*lhs[3][4];
	r[0]   = r[0]   - coeff*r[3];

	coeff = lhs[1][3];
	lhs[1][4]= lhs[1][4] - coeff*lhs[3][4];
	r[1]   = r[1]   - coeff*r[3];

	coeff = lhs[2][3];
	lhs[2][4]= lhs[2][4] - coeff*lhs[3][4];
	r[2]   = r[2]   - coeff*r[3];

	coeff = lhs[4][3];
	lhs[4][4]= lhs[4][4] - coeff*lhs[3][4];
	r[4]   = r[4]   - coeff*r[3];


	pivot = one/lhs[4][4];
	r[4]   = r[4]  *pivot;

	coeff = lhs[0][4];
	r[0]   = r[0]   - coeff*r[4];

	coeff = lhs[1][4];
	r[1]   = r[1]   - coeff*r[4];

	coeff = lhs[2][4];
	r[2]   = r[2]   - coeff*r[4];

	coeff = lhs[3][4];
	r[3]   = r[3]   - coeff*r[4];

}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void y_solve(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     Performs line solves in Y direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix][ 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c-------------------------------------------------------------------*/

	lhsy();
	y_solve_cell();
	y_backsubstitute();
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void y_backsubstitute(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     back solve: if last cell][ then generate U(jsize)=rhs(jsize)
c     else assume U(jsize) is loaded in un pack backsub_info
c     so just use it
c     after call u(jstart) will be sent to next cell
c-------------------------------------------------------------------*/

	int i, j, k, m, n;

	for (j = grid_points[1]-2; j >= 0; j--) {
#pragma omp for
		for (i = 1; i < grid_points[0]-1; i++) {
			for (k = 1; k < grid_points[2]-1; k++) {
				for (m = 0; m < BLOCK_SIZE; m++) {
					for (n = 0; n < BLOCK_SIZE; n++) {
						rhs[i][j][k][m] = rhs[i][j][k][m]
													   - lhs[i][j][k][CC][m][n]*rhs[i][j+1][k][n];
					}
				}
			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void y_solve_cell(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
c-------------------------------------------------------------------*/

	int i, j, k, jsize;

	jsize = grid_points[1]-1;

#pragma omp for  
	for (i = 1; i < grid_points[0]-1; i++) {
		for (k = 1; k < grid_points[2]-1; k++) {

			/*--------------------------------------------------------------------
c     multiply c(i,0,k) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c-------------------------------------------------------------------*/
			binvcrhs( lhs[i][0][k][BB],
					lhs[i][0][k][CC],
					rhs[i][0][k] );
		}
	}

	/*--------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c-------------------------------------------------------------------*/
	for (j = 1; j < jsize; j++) {
#pragma omp for
		for (i = 1; i < grid_points[0]-1; i++) {
			for (k = 1; k < grid_points[2]-1; k++) {

				/*--------------------------------------------------------------------
c     subtract A*lhs_vector(j-1) from lhs_vector(j)
c     
c     rhs(j) = rhs(j) - A*rhs(j-1)
c-------------------------------------------------------------------*/
				matvec_sub(lhs[i][j][k][AA],
						rhs[i][j-1][k], rhs[i][j][k]);

				/*--------------------------------------------------------------------
c     B(j) = B(j) - C(j-1)*A(j)
c-------------------------------------------------------------------*/
				matmul_sub(lhs[i][j][k][AA],
						lhs[i][j-1][k][CC],
						lhs[i][j][k][BB]);

				/*--------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
c-------------------------------------------------------------------*/
				binvcrhs( lhs[i][j][k][BB],
						lhs[i][j][k][CC],
						rhs[i][j][k] );

			}
		}
	}

#pragma omp for  
	for (i = 1; i < grid_points[0]-1; i++) {
		for (k = 1; k < grid_points[2]-1; k++) {

			/*--------------------------------------------------------------------
c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
c-------------------------------------------------------------------*/
			matvec_sub(lhs[i][jsize][k][AA],
					rhs[i][jsize-1][k], rhs[i][jsize][k]);

			/*--------------------------------------------------------------------
c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
c     call matmul_sub(aa,i,jsize,k,c,
c     $              cc,i,jsize-1,k,c,BB,i,jsize,k)
c-------------------------------------------------------------------*/
			matmul_sub(lhs[i][jsize][k][AA],
					lhs[i][jsize-1][k][CC],
					lhs[i][jsize][k][BB]);

			/*--------------------------------------------------------------------
c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
c-------------------------------------------------------------------*/
			binvrhs( lhs[i][jsize][k][BB],
					rhs[i][jsize][k] );

		}
	}
}


/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void z_solve(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     Performs line solves in Z direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c-------------------------------------------------------------------*/

	lhsz();
	z_solve_cell();
	z_backsubstitute();
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void z_backsubstitute(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
c     else assume U(ksize) is loaded in un pack backsub_info
c     so just use it
c     after call u(kstart) will be sent to next cell
c-------------------------------------------------------------------*/

	int i, j, k, m, n;

#pragma omp for  
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {
			for (k = grid_points[2]-2; k >= 0; k--) {
				for (m = 0; m < BLOCK_SIZE; m++) {
					for (n = 0; n < BLOCK_SIZE; n++) {
						rhs[i][j][k][m] = rhs[i][j][k][m]
													   - lhs[i][j][k][CC][m][n]*rhs[i][j][k+1][n];
					}
				}
			}
		}
	}
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

static void z_solve_cell(void) {

	/*--------------------------------------------------------------------
--------------------------------------------------------------------*/

	/*--------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
c-------------------------------------------------------------------*/

	int i,j,k,ksize;

	ksize = grid_points[2]-1;

	/*--------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c-------------------------------------------------------------------*/
#pragma omp for  
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {

			/*--------------------------------------------------------------------
c     multiply c(i,j,0) by b_inverse and copy back to c
c     multiply rhs(0) by b_inverse(0) and copy to rhs
c-------------------------------------------------------------------*/
			binvcrhs( lhs[i][j][0][BB],
					lhs[i][j][0][CC],
					rhs[i][j][0] );

		}
	}

	/*--------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c-------------------------------------------------------------------*/
	for (k = 1; k < ksize; k++) {
#pragma omp for
		for (i = 1; i < grid_points[0]-1; i++) {
			for (j = 1; j < grid_points[1]-1; j++) {

				/*--------------------------------------------------------------------
c     subtract A*lhs_vector(k-1) from lhs_vector(k)
c     
c     rhs(k) = rhs(k) - A*rhs(k-1)
c-------------------------------------------------------------------*/
				matvec_sub(lhs[i][j][k][AA],
						rhs[i][j][k-1], rhs[i][j][k]);

				/*--------------------------------------------------------------------
c     B(k) = B(k) - C(k-1)*A(k)
c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,BB,i,j,k)
c-------------------------------------------------------------------*/
				matmul_sub(lhs[i][j][k][AA],
						lhs[i][j][k-1][CC],
						lhs[i][j][k][BB]);

				/*--------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
c-------------------------------------------------------------------*/
				binvcrhs( lhs[i][j][k][BB],
						lhs[i][j][k][CC],
						rhs[i][j][k] );

			}
		}
	}

	/*--------------------------------------------------------------------
c     Now finish up special cases for last cell
c-------------------------------------------------------------------*/
#pragma omp for  
	for (i = 1; i < grid_points[0]-1; i++) {
		for (j = 1; j < grid_points[1]-1; j++) {

			/*--------------------------------------------------------------------
c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
c-------------------------------------------------------------------*/
			matvec_sub(lhs[i][j][ksize][AA],
					rhs[i][j][ksize-1], rhs[i][j][ksize]);

			/*--------------------------------------------------------------------
c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
c     call matmul_sub(aa,i,j,ksize,c,
c     $              cc,i,j,ksize-1,c,BB,i,j,ksize)
c-------------------------------------------------------------------*/
			matmul_sub(lhs[i][j][ksize][AA],
					lhs[i][j][ksize-1][CC],
					lhs[i][j][ksize][BB]);

			/*--------------------------------------------------------------------
c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
c-------------------------------------------------------------------*/
			binvrhs( lhs[i][j][ksize][BB],
					rhs[i][j][ksize] );

		}
	}
}
