/******************************************************************************
* jlp_complex.c
* Set of routines working with complex numbers
*
* JLP
* Version 21/07/2008
*******************************************************************************/
#ifndef _JLP_COMPLEX /* sentry */
#define JLP_COMPLEX 
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>

typedef struct{
float re, im;
} complex;

complex float_to_cplx(float re0, float im0);
complex cplx_product(complex c1, complex c2);
float cplx_sqmodulus(complex c1);

#ifndef PI
#define PI 3.141592653
#endif
/************************************************************************
* Loads two float values to a complex variable
************************************************************************/
complex float_to_cplx(float re0, float im0)
{
complex u;
u.re = re0;
u.im = im0;
return(u);
}
/************************************************************************
* Computes a complex product 
************************************************************************/
complex cplx_product(complex c1, complex c2)
{
complex u;
u.re = c1.re * c2.re - c1.im * c2.im;
u.im = c1.re * c2.im + c1.im * c2.re;
return(u);
}
/************************************************************************
* Computes the square modulus of a complex variable 
************************************************************************/
float cplx_sqmodulus(complex c1)
{
return(c1.re * c1.re + c1.im * c1.im);
}
#endif  /* EOF sentry */
