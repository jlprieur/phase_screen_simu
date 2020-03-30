/*************************************************************
* zernike 
* Set of routines with Zernike polynomials
*
* JLP
* Version 07-08-2008
************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>
#include <jlp_fftw.h>    /* For FFTs */
#ifdef FINAL_FORM
#include <zernike.h>     /* Prototypes of routines defined in zernike.c */
#endif

#ifndef PI
#define PI 3.14159
#endif

/* Defined in bessel_j.c: */
double bessj(int n, double x);

/* Defined in fourn.c: */
void fourn_for_C_arrays(float *data, int *nn, int ndim, int isign);

/* Contained here: */
int factorial(int j);
double zernike_Fourier_R(int nn_rad, int mm_azim, double rr);
double zernike_R(int nn_rad, int mm_azim, double rr);
double zernike_Z(int mm_azim, double theta);
static int zernike_pupil0(double *array0, int nx0, int ny0, 
                          int patch_radius, int nn_rad, int mm_azim);
static int zernike_fourier0(double *array0, int nx0, int ny0, 
                            int patch_radius, int nn_rad, int mm_azim);
static int zernike_demo(char *out_prefix);
static int copy_patch_to_image(double *array0, int nx0, int ny0, 
                               double *zernike1, int Nx, int Ny, int ixcent, 
                               int iycent, int patch_radius, double backgd);
static int psf_from_pupil(double *array0, int nx0, int ny0);

#define MAIN_PROGRAM
#ifdef MAIN_PROGRAM
int main(int argc, char *argv[])
{
float *phase_screen1; 
int istat;
INT4 nx1, ny1;
INT_PNTR pntr_ima;
char psc_infile[60], zernike_outfile[60], comments[80];
char out_prefix[20];

if(argc == 7) {
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
}
if(argc != 3) {
  fprintf(stderr, "argc=%d\n", argc);
  fprintf(stderr, "zernike/Fatal error\n");
  fprintf(stderr, " Syntax: zernike in_phase_screen out_prefix \n");
  fprintf(stderr, " Example: zernike my_screen tt\n");
  return(-1);
  }

/* Read input parameters: */
strcpy(psc_infile, argv[1]);
strcpy(out_prefix, argv[2]);
sprintf(zernike_outfile, "%s_kernike.dat", out_prefix);

/* Output parameters to check it is OK: */
printf("OK: infile=%s out_prefix=>%s<\n", psc_infile, out_prefix);
printf("Will output profiles to: >%s<\n", zernike_outfile);

JLP_INQUIFMT();

/* Read input phase array: */
   istat=JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, psc_infile, comments);
   phase_screen1 = (float *)pntr_ima;
    if(istat != 0) {
      fprintf(stderr, " Fatal error reading %s istat = %d \n", 
              psc_infile, istat);
      return(-1);
     }

zernike_demo(out_prefix);

return(0);
}
#endif

/***************************************************************************
* Compute an image with patches corresponding to the Zernike polynomials
*
****************************************************************************/
static int zernike_demo(char *out_prefix)
{
int Nx, Ny, ixcent, iycent, nn_max, nn_rad, mm_azim, k, nx0, ny0;
double *zernike1, *zernike2, patch_radius, *array0, backgd;
char out_file[80], out_comments[80];
register int i;

  Nx = 350;
  Ny = 280;
  nn_max = 7;
  nx0 = 32;
  ny0 = nx0;
  patch_radius = nx0/2;
  iycent = 0.94 * Ny - patch_radius;

  array0 = (double *)malloc(nx0 * ny0 * sizeof(double));
  zernike1 = (double *)malloc(Nx * Ny * sizeof(double));
  zernike2 = (double *)malloc(Nx * Ny * sizeof(double));

  backgd = -3.;
  for(i = 0; i < Nx * Ny; i++) zernike1[i] = backgd; 
  for(i = 0; i < Nx * Ny; i++) zernike2[i] = 0.; 

/* Compute Zernike polynomial: */
/* patch_radius = 16 */
  for(nn_rad = 0; nn_rad < nn_max; nn_rad++){
    for(mm_azim = - nn_rad; mm_azim <= nn_rad; mm_azim+=2){
    printf(" index=%d nn_rad=%d mm_azim=%d\n", k, nn_rad, mm_azim);
    ixcent = Nx/2 + patch_radius * mm_azim * 1.1;
    zernike_pupil0(array0, nx0, ny0, patch_radius, nn_rad, mm_azim);
    copy_patch_to_image(array0, nx0, ny0, zernike1, Nx, Ny, ixcent, iycent,
                        patch_radius, backgd);
    psf_from_pupil(array0, nx0, ny0);
    copy_patch_to_image(array0, nx0, ny0, zernike2, Nx, Ny, ixcent, iycent,
                        1000, 0.);
    k++;
    }
    iycent -= patch_radius * 2.;
  }

/* Output the pupil patches to a FITS file:
*/
  sprintf(out_file, "%s_z", out_prefix);
  strcpy(out_comments, "Pupil patches corresponding to Zernike polynomial");
  JLP_D_WRITEIMAG(zernike1, &Nx, &Ny, &Nx, out_file, out_comments);

/* Output the PSF patches to a FITS file:
*/
  sprintf(out_file, "%s_zpsf", out_prefix);
  strcpy(out_comments, "PSF patches corresponding to Zernike polynomial");
  JLP_D_WRITEIMAG(zernike2, &Nx, &Ny, &Nx, out_file, out_comments);

  free(zernike1);
  free(zernike2);
  free(array0);
return(0);
}
/************************************************************************
* Load a circular patch to an image corresponding to zernike polynomial
* of parameters nn_rad and mm_azim
*
************************************************************************/
static int zernike_pupil0(double *array0, int nx0, int ny0, 
                          int patch_radius, int nn_rad, int mm_azim)
{
double rad, angl, dx, dy;
int ixc, iyc;
register int ix, iy;

ixc = nx0/2;
iyc = nx0/2;

for(iy = 0; iy < ny0; iy++) {
  for(ix = 0; ix < nx0; ix++) {
   dx = ix - ixc;
   dy = iy - iyc; 
   rad = sqrt(SQUARE(dx) + SQUARE(dy));
   if(rad < patch_radius) { 
/* atan2(y,x) returns arctang of Y/X: */
     if(dx != 0)
        angl = atan2(dy/rad, dx/rad);
     else if (dy > 0)
        angl = PI/2.;
     else
        angl = -PI/2.;

     rad /= patch_radius; 
     array0[ix + iy * nx0] = 
             zernike_R(nn_rad, mm_azim, rad) * zernike_Z(mm_azim, angl);
   } else {
     array0[ix + iy * nx0] = 0.; 
   }
  }
}

return(0);
}
/************************************************************************
* Add an image corresponding to the Fourier transform of a zernike polynomial
* of parameters nn_rad and mm_azim
*
************************************************************************/
static int zernike_fourier0(double *array0, int nx0, int ny0, 
                            int patch_radius, int nn_rad, int mm_azim)
{
double rad, angl, dx, dy, ww;
int ixc, iyc;
register int ix, iy;


ixc = nx0/2;
iyc = nx0/2;

for(iy = 0; iy < ny0; iy++) {
  for(ix = 0; ix < nx0; ix++) {
   dx = ix - ixc;
   dy = iy - iyc; 
   rad = sqrt(SQUARE(dx) + SQUARE(dy));
   if(rad < patch_radius) { 
/* atan2(y,x) returns arctang of Y/X: */
     if(dx != 0)
        angl = atan2(dy/rad, dx/rad);
     else if (dy > 0)
        angl = PI/2.;
     else
        angl = -PI/2.;

     rad /= patch_radius; 
/* JLP2008: I multiply with 5 to be able to see something ... */
     rad *= 5.;
     ww = zernike_Fourier_R(nn_rad, mm_azim, rad) * zernike_Z(mm_azim, angl);
     array0[ix + iy * nx0] = ww * ww; 
   } else {
     array0[ix + iy * nx0] = 0.; 
   }
  }
}

return(0);
}
/************************************************************************
* Compute the radial component of the Fourier transform 
* of the Zernike polynomial
*
* INPUT:
* nn_rad, mm_azim: indices of FT of the R_nn^mm polynomial
* kk_rad: radial fraquency where the FT of the R_nn^mm polynomial has 
*         to be computed
*
* OUTPUT:
* return the value the FT of R_nn^mm(kk_rad)
************************************************************************/
double zernike_Fourier_R(int nn_rad, int mm_azim, double kk_rad)
{
double ww, sqrt_np1, m1;

ww = 0.;
m1 = pow(-1.,(double)(nn_rad/2 + ABS(mm_azim)));
sqrt_np1 = sqrt((double)(nn_rad + 1));

if(kk_rad == 0.) {
/* ? JLP: to be investigated... ? */
  ww = 0.;
  } else {
/* Bessel_J: the nn th-order Bessel function of the first kind
*/
  ww = m1 * sqrt_np1 * bessj(nn_rad + 1, 2. * PI * kk_rad) / kk_rad;
  }
return(ww);
}
/************************************************************************
* Compute the nn th-order Bessel function of the first kind.
*
* INPUT:
*   nn: order of the Bessel function
*   xx: value for which J_nn(xx) has to be computed
*
* OUTPUT:
*   J_nn(xx)
*
************************************************************************/
/************************************************************************
* Compute the Zernike radial polynomial
*
* INPUT:
* nn_rad, mm_azim: indices of R_nn^mm polynomial
* rr: radius value where the R_nn^mm polynomial has to be evaluated
*
* OUTPUT:
* return the value R_nn^mm(rr)
************************************************************************/
double zernike_R(int nn_rad, int mm_azim, double rr)
{
double ww, sqrt_np1, m1;
int ss_maxi;
register int ss;

if(rr < 0 || rr > 1.){
  fprintf(stderr, "zernike_R/Fatal error wrong value for radius: rr=%f\n", rr);
  exit(-1);
  }

ww = 0.;
m1 = -1.;
sqrt_np1 = sqrt((double)(nn_rad + 1));

ss_maxi = (nn_rad - ABS(mm_azim))/2;
for(ss = 0; ss <= ss_maxi; ss++) {
  m1 *= (-1); 
  ww += (m1 * sqrt_np1 * factorial(nn_rad - ss) 
            * pow(rr, (double)(nn_rad - 2*ss))) / 
        (factorial(ss) * factorial((nn_rad + mm_azim)/2 - ss) 
           * factorial((nn_rad - mm_azim)/2 - ss)); 
  }
return(ww);
}
/************************************************************************
* Compute the Zernike azimuthal polynomial
*
* INPUT:
* mm_azim: index of Z^mm polynomial
* theta: angle value for which the Z^mm polynomial has to be evaluated
*        (in radians)
*
* OUTPUT:
* return the value Z_nn^mm(rr)
************************************************************************/
double zernike_Z(int mm_azim, double theta)
{
double ww;
if(mm_azim == 0) {
  ww = 1.;
} else if(mm_azim < 0) {
  ww = - sqrt(2) * sin((double)(mm_azim * theta);
} else {
  ww = sqrt(2) * cos((double)(mm_azim * theta);
}
return(ww);
}
/**************************************************************************
* Compute the factorial
* j * (j-1) * ... * 2 * 1
**************************************************************************/
int factorial(int j)
{
int k;
register int i;

if(j < 1) return(1);

k= 1;
for(i = j; i > 1; i--) k *= i;

return(k);
}
/****************************************************************************
* Copy subarray array0 to big array zernicke1. 
* The subarray is centered at (ixcent, iycent)
*
* INPUT:
*  patch_radius: radius of the selected area of array0 to be copied to zernicke1
*  backgd: value of the background (outside the circle)
****************************************************************************/
static int copy_patch_to_image(double *array0, int nx0, int ny0, 
                               double *zernike1, int Nx, int Ny, int ixcent, 
                               int iycent, int patch_radius, double backgd)
{
int ix_start, iy_start, patch_radius2, rad2;
register int ix, iy;

ix_start = ixcent - nx0/2;
iy_start = iycent - ny0/2;
if((ix_start < 0) || (iy_start < 0) || (ix_start + nx0 > Nx) 
   || (ix_start + nx0 > Nx)) {
   fprintf(stderr, "copy_patch_to_image/Fatal error: ixcent=%d, iycent=%d\n",
           ixcent, iycent);
   exit(-1);
   } 

patch_radius2 = SQUARE(patch_radius);
for(iy = 0; iy < ny0; iy++) {
  for(ix = 0; ix < nx0; ix++) {
   rad2 = SQUARE(ix - nx0/2) + SQUARE(iy - ny0/2);
   if(rad2 <= patch_radius2) 
     zernike1[ix + ix_start + (iy + iy_start)*Nx] = array0[ix + iy * nx0];
   else
     zernike1[ix + ix_start + (iy + iy_start)*Nx] = backgd;
  }
}
return(0);
}
/****************************************************************************
* Compute psf_from_pupil
*
*  INPUT
*   re1, im1: used as work arrays
*
*  OUTPUT
*   array0: psf is loaded to array0 
*
****************************************************************************/
static int psf_from_pupil(double *array0, int nx0, int ny0)
{
double *phase, ww;
float *data;
int nn[2], nx1, ny1;
register int i, j, k;

nx1 = 2 * nx0;
ny1 = 2 * ny0;

  phase = (double *)malloc(nx1 * ny1 * sizeof(double));
  data = (float *)malloc(2 * nx1 * ny1 * sizeof(float));
  for(i = 0; i < nx1 * ny1; i++) phase[i] = 0.; 

/* Copy to a larger array (to avoid problems at the edges): */
for(j = 0; j < ny0; j++) 
  for(i = 0; i < nx0; i++) 
/* JLPPP ?
    phase[(i + nx0/2) + (j + ny0/2) * nx1] = 2. * PI * array0[i + j * nx0]; 
*/
    phase[(i + nx0/2) + (j + ny0/2) * nx1] = array0[i + j * nx0]; 

RECENT_FFT_DOUBLE(phase, phase, &nx1, &ny1, &nx1);

/*
void fourn_for_C_arrays(float *data, int *nn, int ndim, int isign)
*/
k = 0;
for(i = 0; i < nx1 * ny1; i++) {
   data[k++] = cos(phase[i]);
   data[k++] = sin(phase[i]);
  }

/* Compute the complex amplitude in image plane from that of the pupil plane: */
nn[0] = nx1;
nn[1] = ny1;
fourn_for_C_arrays(data, nn, 2, 1);

/* Compute intensity of image plane from the complex amplitude of the field: */
k = 0;
ww = SQUARE(nx1 * ny1);
for(i = 0; i < nx1 * ny1; i++) {
   phase[i] = SQUARE(data[k]) + SQUARE(data[k+1]) / ww;
   k += 2;
  }

RECENT_FFT_DOUBLE(phase, phase, &nx1, &ny1, &nx1);

/* Copy to array0 before exit: */
for(j = 0; j < ny0; j++) 
  for(i = 0; i < nx0; i++) 
    array0[i + j * nx0] = phase[(i + nx0/2) + (j + ny0/2) * nx1]; 

free(data);
free(phase);
return(0);
}
