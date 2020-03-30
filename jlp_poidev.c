/**************************************************************
* Poisson random generator
* Returns as a floating -point number an integer value that is a random
* deviate drawn from a Poisson distribution of mean xm, using ran1_poidev(idum)
* as a source of uniform random deviates.
*
* From Numerical recipees Chap. 7
*
* JLP Version 05-07-2008
* (revised for ANSI C syntax in 2008)
**************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <jlp_poidev.h>

#ifndef PI
#define PI 3.141592654
#endif

/*
int JLP_POISSON(float *in_val, float *out_val, long *iseed);
float POIDEV(float xm, long *idum);
float ran1_for_poidev(long *idum);
float gammln(float xx);
*/

/********************************************************************
* Fortran interface with poidev
*
********************************************************************/
int JLP_POISSON(float *in_val, float *out_val, long *iseed)
{
*out_val = POIDEV(*in_val, iseed);
return(0);
}

float POIDEV(float xm, long *idum)
{
static float sq,alxm,g,oldm=(-1.0);
float em,t,y;

if (xm < 12.0) {
	if (xm != oldm) {
		oldm=xm;
		g=exp(-xm);
	}
	em = -1;
	t=1.0;
	do {
		em += 1.0;
		t *= ran1_for_poidev(idum);
	} while (t > g);
} else {
	if (xm != oldm) {
		oldm=xm;
		sq=sqrt(2.0*xm);
		alxm=log(xm);
		g=xm*alxm-gammln(xm+1.0);
	}
	do {
		do {
			y=tan(PI*ran1_for_poidev(idum));
			em=sq*y+xm;
		} while (em < 0.0);
		em=floor(em);
		t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
	} while (ran1_for_poidev(idum) > t);
}
return em;
}

/********************************************************************* 
* Natural log of the gamma function (paragraph 6.2)
*********************************************************************/
float gammln(float xx)
{
double x,tmp,ser;
static double cof[6]={76.18009173,-86.50532033,24.01409822,
	-1.231739516,0.120858003e-2,-0.536382e-5};
int j;

x=xx-1.0;
tmp=x+5.5;
tmp -= (x+0.5)*log(tmp);
ser=1.0;
for (j=0;j<=5;j++) {
	x += 1.0;
	ser += cof[j]/x;
}
return -tmp+log(2.50662827465*ser);
}


/* ran1_for_poidev = ran1.c Random generator between 0 and 1 used
* by poidev (Numerical recipees in C p280)
*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1_for_poidev(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
