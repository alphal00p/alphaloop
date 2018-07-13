/*
  demo-c.c
  test program for the Cuba library
  last modified 13 Mar 15 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif


extern int Integrand(const int *ndim, const cubareal xx[],
		     const int *ncomp, cubareal ff[], void *userdata) ;


/*********************************************************************/
#define CUBACORES 1
#define NDIM 4
#define NCOMP 2
#define USERDATA NULL
#define NVEC 1
#define EPSREL 10*1e-3
#define EPSABS 1e-12
#define VERBOSE 2 //2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 5000000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

#ifndef CUBA_WHAT
#define CUBA_WHAT 1
#endif

#if CUBA_WHAT==1
  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	EPSREL, EPSABS, VERBOSE, SEED,
	MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
	GRIDNO, STATEFILE, SPIN,
	&neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
	 neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
	   (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if CUBA_WHAT==2
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	EPSREL, EPSABS, VERBOSE | LAST, SEED,
	MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
	STATEFILE, SPIN,
	&nregions, &neval, &fail, integral, error, prob);

  printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
	 nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
	   (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if CUBA_WHAT==3
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	  EPSREL, EPSABS, VERBOSE, SEED,
	  MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
	  BORDER, MAXCHISQ, MINDEVIATION,
	  NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
	  STATEFILE, SPIN,
	  &nregions, &neval, &fail, integral, error, prob);

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
	 nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
	   (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if CUBA_WHAT==4
  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
	EPSREL, EPSABS, VERBOSE | LAST,
	MINEVAL, MAXEVAL, KEY,
	STATEFILE, SPIN,
	&nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
	 nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
	   (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

  return 0;
}

