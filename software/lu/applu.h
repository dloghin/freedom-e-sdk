/*--------------------------------------------------------------------
c---------------------------------------------------------------------
c---  applu.h
c---------------------------------------------------------------------
c-------------------------------------------------------------------*/

/*--------------------------------------------------------------------
c   npbparams.h defines parameters that depend on the class and 
c   number of nodes
c-------------------------------------------------------------------*/

#include "npbparams.h"

typedef float element_t;
typedef int boolean;

#define TRUE	1
#define FALSE	0

/*--------------------------------------------------------------------
c   grid
c-------------------------------------------------------------------*/

/* common /cgcon/ */
static int nx, ny, nz;
static int nx0, ny0, nz0;
static int ist, iend;
static int jst, jend;
static int ii1, ii2;
static int ji1, ji2;
static int ki1, ki2;
static element_t dxi, deta, dzeta;
static element_t tx1, tx2, tx3;
static element_t ty1, ty2, ty3;
static element_t tz1, tz2, tz3;

/*--------------------------------------------------------------------
c   dissipation
c-------------------------------------------------------------------*/

/* common /disp/ */
static element_t dx1, dx2, dx3, dx4, dx5;
static element_t dy1, dy2, dy3, dy4, dy5;
static element_t dz1, dz2, dz3, dz4, dz5;
static element_t dssp;

/*--------------------------------------------------------------------
c   field variables and residuals
c   to improve cache performance, second two dimensions padded by 1 
c   for even number sizes only.
c   Note: corresponding array (called "v") in routines blts, buts, 
c   and l2norm are similarly padded
c-------------------------------------------------------------------*/

/* common /cvar/ */
static element_t u[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5];
static element_t rsd[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5];
static element_t frct[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5];
static element_t flux[ISIZ1][ISIZ2/2*2+1][ISIZ3/2*2+1][5];

/*--------------------------------------------------------------------
c   output control parameters
c-------------------------------------------------------------------*/

/* common /cprcon/ */
static int ipr, inorm;

/*--------------------------------------------------------------------
c   newton-raphson iteration control parameters
c-------------------------------------------------------------------*/

/* common /ctscon/ */
static int itmax, invert;
static element_t dt, omega, tolrsd[5], rsdnm[5], errnm[5], frc, ttotal;
  
/* common /cjac/ */
static element_t a[ISIZ1][ISIZ2][5][5];
static element_t b[ISIZ1][ISIZ2][5][5];
static element_t c[ISIZ1][ISIZ2][5][5];
static element_t d[ISIZ1][ISIZ2][5][5];

/*--------------------------------------------------------------------
c   coefficients of the exact solution
c-------------------------------------------------------------------*/

/* common /cexact/ */
static element_t ce[5][13];

/*--------------------------------------------------------------------
c   multi-processor common blocks
c-------------------------------------------------------------------*/

/* common /timer/ */
static element_t maxtime;

/*--------------------------------------------------------------------
c   end of include file
c-------------------------------------------------------------------*/



