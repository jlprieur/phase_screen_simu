/***************************************************************************
* poidev.h
* Include file for poidev.c
*
* JLP
* Version 05/08/2008
***************************************************************************/
#ifndef _poisson_include
#define _poisson_include

#include <jlp_ftoc.h>

#define JLP_POISSON RENAME_(jlp_poisson) 
#define POIDEV      RENAME(poidev) 

int JLP_POISSON(float *in_val, float *out_val, long *iseed);
float POIDEV(float xm, long *idum);
float ran1_for_poidev(long *idum);
float gammln(float xx);
#endif
