/*fourn.c 
  This subroutine calculates an n-dimensional FFT. Its taken almost directly
  from Press et al.s Numerical Recipes (pgs 468-470). The only
  difference is that this routine divides by the normalization factor of
  1/(nn[1]*nn[2]*...) after doing an inverse transform.
 
  Note the format of data[] is as follows: a real array of length twice
  the product of the lengths of each dimension in which the data are
  stored with real and imaginary parts of each element in consecutive
  locations. The rightmost index of the array increases most rapidly.
  For a 2-D array this is equivalent to storing the array by rows.

  Pleas note that indices of arrays vary from 1 to n
  instead of from 0 to n-1 in "conventionnal" C programs...

JLP Version 06/07/98
*/
#include <stdio.h>
#include <math.h>
 
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define TWOPI 6.28318530717959

void fourn_(float *data, int *nn, int *ndim, int *isign);
void fourn(float *data, int *nn, int ndim, int isign);
void fourn_for_C_arrays(float *data, int *nn, int ndim, int isign);
void nrerror(char error_text[]);
 
/************************************************************
* Same as fourn, but interface for Fortran calls (by reference)
*
************************************************************/
void fourn_(float *data, int *nn, int *ndim, int *isign)
{
fourn(data, nn, *ndim, *isign);
return;
}
/************************************************************
* Same as fourn, but interface for C arrays 
*
************************************************************/
void fourn_for_C_arrays(float *data, int *nn, int ndim, int isign)
{
fourn(&data[-1], &nn[-1], ndim, isign);
}

/***********************************************************************
* WARNING! arrays start at [1], not at [0]!
*
* data[]:	complex data in row order
*               (re[1], im[1], re[2], im[2], etc...)
* nn[]:	        array containing lengths of each dimension (these must all
*		be powers of 2)
* ndim:	        total number of dimensions
* isign:	direction of transform, 1=forward, -1=reverse
************************************************************************/
void fourn(float *data, int *nn, int ndim, int isign)
{
int i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
int ibit,idim,k1,k2,n,nprev,nrem,ntot;
float tempi,tempr;
double theta,wi,wpi,wpr,wr,wtemp;	/*double prec for trig*/
double norm;				/*recurrences*/
	
ntot = 1;
/* Compute total number of complex values */
for(idim=1; idim<=ndim; idim++)	ntot *= nn[idim];	
nprev = 1;

/*main loop over dimensions*/
for(idim=ndim;idim>=1;idim--) {
	n = nn[idim];
	nrem = ntot/(n * nprev);
	ip1 = nprev << 1;
	ip2 = ip1 * n;
	ip3 = ip2 * nrem;
	i2rev = 1;
/*the bit reversal section*/
	for(i2=1;i2<=ip2;i2+=ip1) {
	   if(i2 < i2rev) {
		for(i1=i2;i1<=i2+ip1-2;i1+=2) {
			for(i3=i1;i3<=ip3;i3+=ip2) {
				i3rev = i2rev + i3 - i2;
				SWAP(data[i3],data[i3rev]);
				SWAP(data[i3+1],data[i3rev+1]);
			}
		}
	    }
	   ibit = ip2 >> 1;
	   while(ibit >= ip1 && i2rev > ibit) {
		i2rev -= ibit;
		ibit >>= 1;
	   }
	   i2rev += ibit;
	}
 
/*Here begins the Danielson-Lanczos section of the routine*/
	ifp1 = ip1;	 
	while(ifp1 < ip2) {
		ifp2 = ifp1 << 1;
		theta = isign * TWOPI / (ifp2 / ip1);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for(i3=1;i3<=ifp1;i3+=ip1) {
	           for(i1=i3;i1<=i3+ip1-2;i1+=2) {
			for(i2=i1;i2<=ip3;i2+=ifp2) {
/* Danielson-Lanczos formula */
			   k1 = i2;	
			   k2 = k1 + ifp1;
			   tempr = wr * data[k2] - wi * data[k2+1];
			   tempi = wr * data[k2+1] + wi * data[k2];
			   data[k2] = data[k1] - tempr;
			   data[k2+1] = data[k1+1] - tempi;
			   data[k1] += tempr;
			   data[k1+1] += tempi;
			}
		  }
/* Trigonometric recurrence: */
		   wr = (wtemp=wr) * wpr - wi * wpi + wr;
		   wi = wi * wpr + wtemp * wpi + wi;
		}
		ifp1 = ifp2;
	}
	nprev *= n;
}

/* Check if inverse FFT: */
if(isign == -1) {
	norm = ntot;
	for(idim=1;idim<=2*ntot;idim++) data[idim] /= norm;
}

return;
}
/*******************************************************************/
/* BOF TOTO, just for info: */
#ifdef TOTO
typedef struct FCOMPLEX {float r,i;} fcomplex;
typedef struct IMMENSE {unsigned long l,r;} immense;
typedef struct GREAT {unsigned short l,c,r;} great;
 
#ifdef ANSI
	void adi(double **a, double **b, double **c, double **d, double **e,
		double **f, double **g, double **u, int jmax, int k,
		double alpha, double beta, double eps);
	void amoeba(float **p, float *y, int ndim, float ftol,
		float (*funk)(float *), int *iter);
	float bessj0(float x);
	void four1(float *data, int nn, int isign);
	void fourn(float *data, int *nn, int ndim, int isign);
#else
	void adi();
	void amoeba();
	float bessj0();
	void four1();
	void fourn();
#endif
free_ivector(v,nl,nh)
int *v;
int nl,nh;
{
	free((char *) (v+nl));
}
/*Frees a matrix allocated with matrix()*/
 
void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl;
int nrh;
int ncl;
int nch;
{
	int i;
 
	for(i=nrh;i>=nrl;i--)
		free((char *) (m[i] + ncl));
	free((char *) (m + nrl));
}
free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char *) (v+nl));
}
#include <malloc.h>
 
int *ivector(nl,nh)
int nl,nh;
{
	int *v;
 
	v = (int *) malloc((unsigned) (nh - nl +1) * sizeof(int));
	if(!v)
		nrerror("allocation failure in ivector()");
	return(v - nl);
}
/*allocate a float matrix with range [nrl..nrh][ncl..nch]*/
/* #include <malloc.h> */
 
float **matrix(nrl,nrh,ncl,nch)
int nrl;	/*low range of rows*/
int nrh;	/*high range of rows*/
int ncl;	/*low range of cols*/
int nch;	/*high range of cols*/
{
	int i;
	float **m;
 
	/*allocate pointers to rows*/
	m = (float **) malloc((unsigned) (nrh - nrl + 1) * sizeof(float *));
	if(!m)
		nrerror("allocation failure 1 in matrix()");
	m -= nrl;
 
	/*allocate rows and set pointers to them*/
	for(i=nrl;i<=nrh;i++)
	{
		m[i] = (float *) malloc((unsigned) (nch - ncl + 1) *
			sizeof(float));
		if(!m[i])
			nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	/*return pointer to array of pointers to rows*/
	return m;
}
 
void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error:\n");
	fprintf(stderr,"%s\n",error_text);
	exit(1);
}
/* #include <malloc.h> */
 
float *vector(nl,nh)
int nl,nh;
{
	float *v;
 
	v= (float *) malloc((unsigned) (nh -nl + 1) * sizeof(float));
	if(!v)
		nrerror("allocation failure in vector()");
	return(v - nl);
}
#ifdef ANSI
	void nrerror(char *error_text);
	float *vector(int nl, int nh);
	int *ivector(int nl, int nh);
	double *dvector(int nl, int nh);
	float **matrix(int nrl, int nrh, int ncl, int nch);
	int **imatrix(int nrl, int nrh, int ncl, int nch);
	double **dmatrix(int nrl, int nrh, int ncl, int nch);
	float **submatrix(float **a, int oldrl, int oldrh, int oldcl,
		int oldch, int newrl, int newcl);
	float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
	void free_vector(float *v, int nl, int nh);
	void free_ivector(int *v, int nl, int nh);
	void free_dvector(double *v, int nl, int nh);
	void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
	void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
	void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
	void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch);
	void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch);
#else
/*	void nrerror(); */
	float *vector();
	int *ivector();
	double *dvector();
	float **matrix();
	int **imatrix();
	double **dmatrix();
	float **submatrix();
	float **convert_matrix();
/*	void free_vector();
	void free_ivector(); */
	void free_matrix();
	void free_imatrix();
	void free_dmatrix();
	void free_submatrix();
	void free_convert_matrix();
#endif
/* EOF TOTO */
#endif
