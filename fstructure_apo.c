/*************************************************************
*  fstructure.c 
*  Set of routines to compute the structure function of a phase array 
*  called by simu3.c and by s_pscorg.for
*
*  JLP
*  Version 11-07-2008
************************************************************/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>         /* For file handling */
#include <jlp_ftoc.h>
#include <jlp_fftw.h>       /* For FFTs */
#include <fstructure.h>     /* Prototypes of routines defined in fstructure.c */

static float blackman(float argx);
static float hamming(float argx);
int APODIZE_FLOAT_ARRAY(float *in_aa, float *out_aa, int *nx, int *ny, 
                        int *idim, float *width1, float *filter_effic);
int extract_subarray(float *in_aa, int *nx, int *ny, int *idim);

#ifdef MAIN_PROGRAM
int main(int argc, char *argv[])
{
float *phase_array, *fstruct_array, r0_pix, width1, filter_effic;
int istat, fourn_format;
INT4 nx, ny, idim;
INT_PNTR pntr_ima;
char phase_infile[60], fstruct_outfile[60], comments[120];
register int i, j;

if(argc == 7 && *argv[4]) argc = 5;
if(argc == 7 && *argv[3]) argc = 4;
if(argc == 7 && *argv[2]) argc = 3;
if(argc == 7 && *argv[1]) argc = 2;
if(argc != 3) {
  fprintf(stderr, "argc=%d\n", argc);
  fprintf(stderr, "fstructure/Fatal error\n");
  fprintf(stderr, " Syntax: fstructure in_phase_screen out_structure_file\n");
  return(-1);
  }

/* Read input parameters: */
strcpy(phase_infile, argv[1]);
strcpy(fstruct_outfile, argv[2]);
printf("OK: infile=%s outfile=%s\n", phase_infile, fstruct_outfile);

JLP_INQUIFMT();

/* Read input phase array: */
   istat=JLP_VM_READIMAG1(&pntr_ima, &nx, &ny, phase_infile, comments);
   phase_array = (float *)pntr_ima;
    if(istat != 0) {
      fprintf(stderr, " Fatal error reading %s istat = %d \n", 
              phase_infile, istat);
      return(-1);
     }
  idim = nx;

/* Case of simulations with arroyo: 276x281....
 and the input image is not a power of two: 
*/
  if(nx == 276 && ny == 281) {
   idim = nx;
   nx = 256;
   ny = 256; 
   } 

/* Allocation of memory for structure function: */
  fstruct_array = (float *)malloc(nx * ny * sizeof(float));
 
/* First extract the central part of the array: */
 extract_subarray(phase_array, &nx, &ny, &idim);
 nx /= 2;
 ny /= 2;

/* Then apodize the phase array: */
  width1 = 1.1 * (float)nx;
  APODIZE_FLOAT_ARRAY(phase_array, phase_array, &nx, &ny, &idim, 
                      &width1, &filter_effic);
  printf("With apodization\n");

/* Compute the structure function of the real phase array: 
*/
   fourn_format = 0;
/* Two (identical) ways of computing the structure function:
* either directly (but very slow):
   STRUCTURE_REAL(phase_array, &nx, &ny, &nx, fstruct_array, &fourn_format);

or with Fourier Transform (fast):
(I checked with a mas 0,1 that the average value as about 0.5 which is OK)
*/
   STRUCTURE_REAL_WITH_FT(phase_array, &nx, &ny, &idim, fstruct_array, 
                          &fourn_format);

/* Correct the structure function taking into account the filter factor: */
  for(j = 0; j < ny; j++)
     for(i = 0; i < nx; i++)
            fstruct_array[i + j * idim] /= filter_effic;

/* Compute r_0 from structure function: 
  R0_FROM_STRUCTURE(fstruct_array, &nx, &ny, &nx, &r0_x, &r0_y);
*/
/*
* JLP2008: I use r0=10cm, fov=9" et 256x256 pixels in simu_close_binary.input...
*/
  r0_pix = 6.7128; 
  CHECK_STRUCTURE_FUNCT(fstruct_array, &nx, &ny, &idim, &r0_pix);

/* Save output to file: */
   printf(" Output of structure function of the phase to %s \n", 
          fstruct_outfile);
   sprintf(comments,"Structure function of %s", phase_infile);
   JLP_WRITEIMAG(fstruct_array, &nx, &ny, &idim, fstruct_outfile, comments);

   strcpy(fstruct_outfile,"test");
   printf(" Output of apodized screen to %s \n", 
          fstruct_outfile);
   sprintf(comments,"Apodized version of %s", phase_infile);
   JLP_WRITEIMAG(phase_array, &nx, &ny, &idim, fstruct_outfile, comments);

free(fstruct_array);
return(0);
}
#endif

/**********************************************************************
* Compute r0 along X and Y axes
*
***********************************************************************/
int  R0_FROM_STRUCTURE(float *fstruct, int *nx, int *ny, int *idim, 
                       float *r0_x, float *r0_y)
{
double sum, D_x, D_y;
float *profile_x, *profile_y;
int ixc, iyc, npts;
FILE *fp;
register int i;

profile_x = (float *)malloc(*nx * sizeof(float));
profile_y = (float *)malloc(*ny * sizeof(float));

/* Compute the cuts along the 0x and Oy axes: */
ixc = (*nx)/2;
iyc = (*ny)/2;
for(i = 0; i < ixc; i++) profile_x[i] = fstruct[i + ixc + (*idim)*iyc];
for(i = 0; i < iyc; i++) profile_y[i] = fstruct[ixc + (*idim)*(i + iyc)];

/* Save this cut to a profile (for a diagnostic) */
if((fp = fopen("fstruct_profile.dat","w")) == NULL){
  fprintf(stderr,"RO_FROM_STRUCTURE/Error opening output profile file\n");
  free(profile_x);
  free(profile_y);
  return(-1);
  }

/* Computes r0_x: */
sum = 0.;
npts = 0;
for(i = 1; i < ixc/5; i++){ 
   sum += log10(profile_x[i]) - (5./3.) * log10((double)i) - log10(6.88);
   npts++;
   }
sum /= (float)npts;
*r0_x = pow(10.,sum / (-5./3.));
printf("R0_FROM_STRUCTURE:  r0_x=%f ", *r0_x);  

/* Computes r0_y: */
sum = 0.;
npts = 0;
for(i = 1; i < iyc/5; i++){ 
   sum += log10(profile_y[i]) - (5./3.) * log10((double)i) - log10(6.88);
   npts++;
   }
sum /= (float)npts;
*r0_y = pow(10.,sum / (-5./3.));
printf(" r0_y=%f\n", *r0_y);  

fprintf(fp,"%% X, cut_x, cut_y, Kolmo_x kolmo_y\n");

/* Don't save the first point (0,0) since pb with logarithm... */
for(i = 1; i < ixc; i++) { 
  D_x = log10(6.88) + (5./3.) * log10((double)i/(*r0_x));
  D_x = pow(10.,D_x);
  D_y = log10(6.88) + (5./3.) * log10((double)i/(*r0_y));
  D_y = pow(10.,D_y);
  fprintf(fp,"%d %f %f %f %f\n", i, profile_x[i], profile_y[i], D_x, D_y);
}

fclose(fp);
free(profile_x);
free(profile_y);
return(0);
}
/**********************************************************************
* Check if r0 is compatible with cuts of structure function along X and Y axes
* 
* INPUT:
*  r0: theoretical value of r0 (value in pixels)
***********************************************************************/
int  CHECK_STRUCTURE_FUNCT(float *fstruct, int *nx, int *ny, int *idim, 
                           float *r0)
{
double sum, D_x;
float *profile_x, *profile_y, r0_x, r0_y;
int ixc, iyc, npts;
FILE *fp;
register int i;

profile_x = (float *)malloc(*nx * sizeof(float));
profile_y = (float *)malloc(*ny * sizeof(float));

/* Compute the cuts along the 0x and Oy axes: */
ixc = (*nx)/2;
iyc = (*ny)/2;
for(i = 0; i < ixc; i++) profile_x[i] = fstruct[i + ixc + (*idim)*iyc];
for(i = 0; i < iyc; i++) profile_y[i] = fstruct[ixc + (*idim)*(i + iyc)];

/* Save this cut to a profile (for a diagnostic) */
if((fp = fopen("fstruct_profile.dat","w")) == NULL){
  fprintf(stderr,"RO_FROM_STRUCTURE/Error opening output profile file\n");
  free(profile_x);
  free(profile_y);
  return(-1);
  }

/* Computes r0_x: */
sum = 0.;
npts = 0;
for(i = 1; i < ixc/5; i++){ 
   sum += log10(profile_x[i]) - (5./3.) * log10((double)i) - log10(6.88);
   npts++;
   }
sum /= (float)npts;
r0_x = pow(10.,sum / (-5./3.));
printf("R0_FROM_STRUCTURE:  r0_x=%f ", r0_x);  

/* Computes r0_y: */
sum = 0.;
npts = 0;
for(i = 1; i < iyc/5; i++){ 
   sum += log10(profile_y[i]) - (5./3.) * log10((double)i) - log10(6.88);
   npts++;
   }
sum /= (float)npts;
r0_y = pow(10.,sum / (-5./3.));
printf(" r0_y=%f (theoretical value=%f)\n", r0_y, *r0);  

fprintf(fp,"%% X, cut_x, cut_y, Kolmo_profile (with r0=%f)\n", *r0);

/* Don't save the first point (0,0) since pb with logarithm... */
for(i = 1; i < ixc; i++) { 
  D_x = log10(6.88) + (5./3.) * log10((double)i/(*r0));
  D_x = pow(10.,D_x);
  fprintf(fp,"%d %f %f %f\n", i, profile_x[i], profile_y[i], D_x);
}

fclose(fp);
free(profile_x);
free(profile_y);
return(0);
}
/****************************************************************************
* Compute structure function with Fourier Transform (to go faster)
*
****************************************************************************/
int  STRUCTURE_REAL_WITH_FT(float *phase_array, int *nx, int *ny, 
                            int *idim, float *fstruct, int *fourn_format)
{
float *ree, *imm, scale;
int size, ifourn, kod;
register int i, j, k;

ifourn = (*fourn_format == 1) ? 2 : 1;
size = (*nx) * (*ny);

ree = (float *)malloc(size * sizeof (float));
imm = (float *)malloc(size * sizeof (float));

for(k = 0; k < size; k++) imm[k] = 0.;

k = 0;
for(j = 0; j < *ny; j++)
 for(i = 0; i < *nx; i++)
  ree[k++] = phase_array[(i + j * (*idim)) * ifourn];

/* Fourier Transform of phase screen: */
/* Warning: FFT_2D calls FOURN2 that calls FOURN1... */
   kod = 1;
   FFT_2D(ree, imm, nx, ny, nx, &kod);

/* Compute square modulus: */
  for(k = 0; k < size; k++) ree[k] = SQUARE(ree[k]) + SQUARE(imm[k]);

/* Set imaginary part to zero before inverse FT: */
  for(k = 0; k < size; k++) imm[k] = 0.;

/* Autocorrelation is inverse Fourier Transform of square modulus:
*/
   kod = -1;
   FFT_2D(ree, imm, nx, ny, nx, &kod);

/* Normalization by nx * ny: */
  scale = (*nx)*(*ny);
  for(k = 0; k < size; k++) ree[k] /= scale;

/* Fstruct = 2 (autoc(0) - autoc(x)) */
k = 0;
for(j = 0; j < *ny; j++)
 for(i = 0; i < *nx; i++)
  fstruct[i + j * (*idim)] = 2. * (ree[0] - ree[k++]);

/* Recenter image: */
 RECENT_FFT(fstruct, fstruct, nx, ny, idim);

free(ree);
free(imm);
return(0);
}

/*********************************************************
* Compute structure function of the phase of the real part of a phase_array
* (stored in the format compatible with "fourn1")
*
* INPUT:
* fourn_format: flag set to one if fourn format is used
*               (i.e. real part is 2 * index, 
*                     imag part is 2 * index + 1 )
*
* OUTPUT:
* fstruct: struct function
**********************************************************/
int  STRUCTURE_REAL(float *phase_array, int *nx, int *ny, int *idim, 
                    float *fstruct, int *fourn_format)
{
int i, j, k, kj, nx1, ny1;
int i1, j1, k1, kj1, ii1, jj1;
int k2, kj2, ifourn;
float npoints;
double sum;

ifourn = (*fourn_format == 1) ? 2 : 1;

nx1 = *nx;
ny1 = *ny;
/***** main loop: computing the upper quadrant ****/
   for(j=0; j<ny1/2; j++)   
   {
    kj = (j + ny1/2) * (*idim);
    for(i=-nx1/2; i<nx1/2; i++)
     {
       k = kj + i + nx1/2;
       sum = 0.;
/**** Internal loop: compute the (i,j) term */
         for(j1=0; j1<ny1; j1++)   
          {
           kj1 = j1 * (*idim);
           jj1 = j1 + j;
/* Makes the input real phase_array periodic with period (nx,ny): */
           jj1 = (jj1 >= ny1) ? jj1-ny1 : jj1;
           jj1 = (jj1 < 0) ? jj1+ny1 : jj1;
           kj2 = jj1 * (*idim);

           for(i1=0; i1<nx1; i1++)
             {
              k1 = kj1 + i1;
              ii1 = i1 + i;
/* Makes the input real phase_array periodic with period (nx,ny): */
              ii1 = (ii1 >= nx1) ? ii1-nx1 : ii1;
              ii1 = (ii1 < 0) ? ii1+nx1 : ii1;
              k2 = kj2 + ii1;
/* Structure term for (i1,j1) and (i1+i,j1+j) : */
/* Real part is found in 2*k1, 2*k2 if fourn format: */
/* Real part is found in k1, k2 if plain real array: */
              sum = sum + SQUARE(phase_array[ifourn*k1] - phase_array[ifourn*k2]); 
             }
           }
       if(i == 0) printf("fstruct[%d,%d]=%f\n", i, j, sum);
       fstruct[k] = (float)sum; 
     }
   }

/* Symmetry relative to the center, so fills the lower quadrant */
   for(j=0; j<ny1/2; j++)   
   {
    for(i=-nx1/2; i<nx1/2; i++)
     {
       k = (i + nx1/2) + (j + ny1/2) * (*idim);
       k1 = (nx1/2 - i) + (ny1/2 - j) * (*idim);
       fstruct[k1] = fstruct[k];
     }
   }

/* Division by npoints = nx*ny : */
 npoints = (float)(nx1 * ny1);
 for(j=0; j<ny1; j++)
   for(i=0; i<nx1; i++)
     fstruct[ i + j * (*idim)] /= npoints;

return(0);
}
/***************************r******************************
* Compute structure function of the phase stored in a  complex array
* (stored in the format compatible with "fourn1")
*
**********************************************************/
int  STRUCTURE_CPLX(float *phase_array, int *nx, int *ny, int *idim, float *fstruct)
{
int i, j, k, kj;
int i1, j1, k1, kj1, ii1, jj1;
int k2, kj2, nx1, ny1;
float npoints;
double sum;

nx1 = *nx;
ny1 = *ny;
/***** main loop: ****/
   for(j=-ny1/2; j<ny1/2; j++)   
   {
    kj = (j + ny1/2) * (*idim);
    for(i=-nx1/2; i<nx1/2; i++)
     {
       k = kj + i + nx1/2;
       sum = 0.;
/**** Internal loop: */
         for(j1=0; j1<ny1; j1++)   
          {
           kj1 = j1 * (*idim);
           jj1 = j1 + j;
/* Makes the complex phase_array periodic with period (nx,ny): */
           jj1 = (jj1 >= ny1) ? jj1-ny1 : jj1;
           jj1 = (jj1 < 0) ? jj1+ny1 : jj1;
           kj2 = jj1 * (*idim);

           for(i1=0; i1<nx1; i1++)
             {
              k1 = kj1 + i1;
              ii1 = i1 + i;
/* Makes the complex phase_array periodic with period (nx,ny): */
              ii1 = (ii1 >= nx1) ? ii1-nx1 : ii1;
              ii1 = (ii1 < 0) ? ii1+nx1 : ii1;
              k2 = kj2 + ii1;
/* Structure term for (i1,j1) and (i1+i,j1+j) : */
              sum = sum + SQUARE(phase_array[2*k1] - phase_array[2*k2]) 
                        + SQUARE(phase_array[2*k1+1] - phase_array[2*k2+1]); 
             }
           }
       fstruct[k] = (float)sum; 
     }
   }

/* Division by npoints = nx*ny : */
 npoints = (float)(nx1 * ny1);
 for(j=0; j<ny1; j++)
   for(i=0; i<nx1; i++)
    fstruct[i + j * (*idim)] /= npoints;

return(0);
}
int extract_subarray(float *in_aa, int *nx, int *ny, int *idim)
{
float *out;
int ixc, iyc;
register int i, j, k;
out = (float*)malloc((*nx) * (*ny) * sizeof(float));

ixc = (*nx)/2;
iyc = (*ny)/2;
k = 0;
for(j = 0; j < iyc; j++)
  for(i = 0; i < ixc; i++){
    out[k++] = in_aa[i + ixc/2 + (j + iyc/2) * (*idim)];
    }

k = 0;
for(j = 0; j < iyc; j++)
  for(i = 0; i < ixc; i++){
    in_aa[i + j * (*idim)] = out[k++];
    }

free(out);
return(0);
}
/*
#define APODI hamming 
#define APODI_NAME "Hamming" 
*/
#define APODI blackman 
#define APODI_NAME "Blackman" 
#ifndef PI
#define PI 3.14159
#endif
/**********************************************************************
* Check if r0 is compatible with cuts of structure function along X and Y axes
* 
* INPUT:
*  in_aa: input float array
*  r0: theoretical value of r0 (value in pixels)
*  width1: width of filter
*
* OUTPUT:
*  out_aa: output apodized array
*  filter_effic: ratio of output intensity/input intensity assuming input=1
***********************************************************************/
int APODIZE_FLOAT_ARRAY(float *in_aa, float *out_aa, int *nx, int *ny, 
                        int *idim, float *width1, float *filter_effic)
{
float width_x, width_y, work;
float argx, argy, flux;
int icent, jcent;
register int i, j;

width_x = MAXI(1., (*width1)/2.);
width_y = MAXI(1., (*width1)/2.);

icent = (*nx)/2;
jcent = (*ny)/2;

flux = 0.;
  for(j = 0; j < (*ny); j++) {
     argy = (float)(j - jcent) / width_y;
      if(argy < -1. || argy > 1.)
        work = 0.;
      else
        work = APODI(PI * argy);
     for(i = 0; i < (*nx); i++)
        {
         argx = (float)(i - icent) / width_x;
          if(argx < -1. || argx > 1.)
            out_aa[i + j * (*idim)] = 0.;
          else
            {
            out_aa[i + j * (*idim)] = in_aa[i + j * (*idim)] 
                                      * work * APODI(PI * argx);
            flux += work * APODI(PI * argx);
            }
        }
     }

*filter_effic = flux / (float)((*nx)*(*ny));
printf("Efficieny factor=%f\n", *filter_effic);

return(0);
}
/******************************************
* Hamming Filter
******************************************/
static float hamming(float argx)
{
return(0.54 + 0.46 * cos((double)argx));
}
/******************************************
* Blackman Filter
******************************************/
static float blackman(float argx)
{
return(0.42 + 0.5 * cos((double)argx) + 0.08 * cos( 2 * (double)argx));
}
