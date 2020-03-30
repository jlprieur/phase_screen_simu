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

#ifndef PI
#define PI 3.14159
#endif

static int compute_fstruct(float *phase_screen1, float *phase_screen2, 
                           float *fstruct_array, float *prof_x, float *prof_y,
                           int nx1, int ny1, int idim, int *nprof, 
                           int i_apodization, int i_subarray, 
                           int i_output_images, char *psc_infile, 
                           char *psc_outfile, char *fstruct_outfile, 
                           char *profile_outfile);
static int profiles_fstruct(float *fstruct, int nx, int ny, int idim, 
                            float *prof_x, float *prof_y, int nprof);
static int extract_subarray(float *in_aa, int in_nx, int in_ny, int in_idim,
                            float *out_aa, int *out_nx, int *out_ny, 
                            int out_idim);
static int copy_array(float *in_aa, int in_nx, int in_ny, int in_idim,
                      float *out_aa, int *out_nx, int *out_ny, 
                      int out_idim);
static int r0_from_fstructure(float *prof_x, float *prof_y, int nprof, 
                              float *r0_x, float *r0_y, float r0_theor_pix,
                              char *profile_outfile);
static float blackman(float argx);
static float hamming(float argx);
int APODIZE_FLOAT_ARRAY(float *in_aa, float *out_aa, int *nx, int *ny, 
                        int *in_idim,int *out_idim, float *width1, 
                        float *filter_effic);

#ifdef MAIN_PROGRAM
int main(int argc, char *argv[])
{
float *phase_screen1, *phase_screen2, *fstruct_array; 
float *prof_x, *prof_y, ww;
int istat, i_subarray, i_apodization, i_output_images, n_pow2, nprof;
INT4 nx1, ny1, idim;
INT_PNTR pntr_ima;
char psc_infile[60], fstruct_outfile[60], comments[80];
char psc_outfile[60], profile_outfile[60], out_prefix[60];
register int i, j;

if(argc == 7) {
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
}
if(argc != 4) {
  fprintf(stderr, "argc=%d\n", argc);
  fprintf(stderr, "fstructure/Fatal error\n");
  fprintf(stderr, " Syntax: fstructure in_phase_screen out_prefix [Options]\n");
  fprintf(stderr, "Options: i_subarray,i_apodization,i_output_images\n");
  fprintf(stderr, " Example: fstructure my_screen tt 1,0,1 \n");
  return(-1);
  }

/* Read input parameters: */
strcpy(psc_infile, argv[1]);
strcpy(out_prefix, argv[2]);
sscanf(argv[3], "%d,%d,%d", &i_subarray, &i_apodization, &i_output_images);
sprintf(fstruct_outfile, "%s_struct.fits", out_prefix);
sprintf(psc_outfile, "%s_psc.fits", out_prefix);
sprintf(profile_outfile, "%s_struct_prof.dat", out_prefix);

/* Output parameters to check it is OK: */
printf("OK: infile=%s out_prefix=>%s<\n", psc_infile, out_prefix);
printf(" i_subarray=%d i_apodization=%d i_output_images=%d\n",
       i_subarray,i_apodization,i_output_images);
printf("Will output profiles to: >%s<\n", profile_outfile);
if(i_output_images){
printf("Will output phase screen and structure function to:\n >%s< and >%s<\n", 
        psc_outfile, fstruct_outfile);
}

JLP_INQUIFMT();

/* Read input phase array: */
   istat=JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, psc_infile, comments);
   phase_screen1 = (float *)pntr_ima;
    if(istat != 0) {
      fprintf(stderr, " Fatal error reading %s istat = %d \n", 
              psc_infile, istat);
      return(-1);
     }
  idim = nx1;

/* Case of simulations with arroyo: 276x281....
 and the input image is not a power of two: 
*/
  ny1 = MINI(nx1, ny1);
  nx1 = ny1;
  ww = log((double)nx1) / log(2.);
  printf(" log2(nx1) = %f\n", ww);
  n_pow2 = (int)ww;
  n_pow2 = pow(2.,(double)n_pow2);
  printf(" log2(nx1) = %f n_pow2=%d\n", ww, n_pow2);
  if(nx1 != n_pow2) {
   nx1 = n_pow2;
   ny1 = n_pow2; 
/* Conversion from index fluctuations to phase fluctuations: 
* phi = k n = 2 PI n / lambda
*/
    for(j = 0; j < ny1; j++)
      for(i = 0; i < nx1; i++)
        phase_screen1[i + j * idim] *= (2. * PI / 0.5e-6); 
   }

/* Allocation of memory for structure function and profiles: */
  fstruct_array = (float *)malloc(idim * idim * sizeof(float));
  phase_screen2 = (float *)malloc(idim * idim * sizeof(float));
  prof_x = (float *)malloc(idim * sizeof(float));
  prof_y = (float *)malloc(idim * sizeof(float));

  compute_fstruct(phase_screen1, phase_screen2, fstruct_array, prof_x, prof_y,
                  nx1, ny1, idim, &nprof, i_apodization, i_subarray, 
                  i_output_images, psc_infile, psc_outfile, fstruct_outfile,
                  profile_outfile);

free(fstruct_array);
free(phase_screen2);
free(prof_x);
free(prof_y);
return(0);
}
#endif

/**********************************************************************
* Compute r0 along X and Y axes
*
* INPUT:
*  r0_theor_pix: theoretical value of r0 (value in pixels)
*  profile_outfile: name of output profile
***********************************************************************/
static int r0_from_fstructure(float *prof_x, float *prof_y, int nprof, 
                              float *r0_x, float *r0_y, float r0_theor_pix,
                              char *profile_outfile)
{
double sum, D_x, D_y, D_z;
int npts;
FILE *fp;
register int i;

/* Save this cut to a profile (for a diagnostic) */
if((fp = fopen(profile_outfile,"w")) == NULL){
  fprintf(stderr,"RO_FROM_FSTRUCTURE/Fatal error opening output profile file\n");
  exit(-1);
  }

/* Computes r0_x: */
sum = 0.;
npts = 0;
for(i = 1; i < MAXI(10,nprof/20); i++){ 
   sum += log10(prof_x[i]) - (5./3.) * log10((double)i) - log10(6.88);
   npts++;
   }
sum /= (float)npts;
*r0_x = pow(10.,sum / (-5./3.));
printf("r0_from_fstructure:  r0_x=%f ", *r0_x);  

/* Computes r0_y: */
sum = 0.;
npts = 0;
for(i = 1; i < MAXI(10,nprof/20); i++){ 
   sum += log10(prof_y[i]) - (5./3.) * log10((double)i) - log10(6.88);
   npts++;
   }
sum /= (float)npts;
*r0_y = pow(10.,sum / (-5./3.));
printf(" r0_y=%f\n", *r0_y);  

fprintf(fp,"%% X, cut_x, cut_y, Kolmo_x kolmo_y Kolmo_theo\n");

/* Don't save the first point (0,0) since pb with logarithm... */
for(i = 1; i < nprof; i++) { 
  D_x = log10(6.88) + (5./3.) * log10((double)i/(*r0_x));
  D_x = pow(10.,D_x);
  D_y = log10(6.88) + (5./3.) * log10((double)i/(*r0_y));
  D_y = pow(10.,D_y);
  D_z = log10(6.88) + (5./3.) * log10((double)i/r0_theor_pix);
  D_z = pow(10.,D_z);
  fprintf(fp,"%d %f %f %f %f %f\n", i, prof_x[i], prof_y[i], D_x, D_y, D_z);
}

fclose(fp);
return(0);
}
/**********************************************************************
* Extract profiles (here cuts) of the structure function along X and Y axes
* 
* INPUT:
*  fstruct: structure function
*  nx, ny: size of structure function
*  idim: first dimension of structure array, i.e. fstruct(idim,*)
*  nprof: number of points required for the profiles 
*
* OUTPUT:
*  prof_x: X profile corresponding to a cut in the middle line
*  prof_x: Y profile corresponding to a cut in the middle column 
*
***********************************************************************/
static int profiles_fstruct(float *fstruct, int nx, int ny, int idim, 
                            float *prof_x, float *prof_y, int nprof)
{
int ixc, iyc;
register int i;

/* Compute the cuts along the 0x and Oy axes: */
ixc = nx/2;
iyc = ny/2;
if(nprof > MAXI(ixc, iyc)){
  fprintf(stderr, "profiles_fstruct/Fatal error: nprof=%d whereas ixc=%d iyc=%d\n", 
          nprof, ixc, iyc);
  exit(-1);
  }
for(i = 0; i < nprof; i++) prof_x[i] = fstruct[i + ixc + idim * iyc];
for(i = 0; i < nprof; i++) prof_y[i] = fstruct[ixc + idim * (i + iyc)];

return(0);
}
/****************************************************************************
* Compute structure function with Fourier Transform (to go faster)
*
****************************************************************************/
int  FSTRUCTURE_REAL_WITH_FT(float *phase_screen, int *nx, int *ny, 
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
  ree[k++] = phase_screen[(i + j * (*idim)) * ifourn];

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
* Compute structure function of the phase of the real part of a phase_screen
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
int  FSTRUCTURE_REAL(float *phase_screen, int *nx, int *ny, int *idim, 
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
/* Makes the input real phase_screen periodic with period (nx,ny): */
           jj1 = (jj1 >= ny1) ? jj1-ny1 : jj1;
           jj1 = (jj1 < 0) ? jj1+ny1 : jj1;
           kj2 = jj1 * (*idim);

           for(i1=0; i1<nx1; i1++)
             {
              k1 = kj1 + i1;
              ii1 = i1 + i;
/* Makes the input real phase_screen periodic with period (nx,ny): */
              ii1 = (ii1 >= nx1) ? ii1-nx1 : ii1;
              ii1 = (ii1 < 0) ? ii1+nx1 : ii1;
              k2 = kj2 + ii1;
/* Structure term for (i1,j1) and (i1+i,j1+j) : */
/* Real part is found in 2*k1, 2*k2 if fourn format: */
/* Real part is found in k1, k2 if plain real array: */
              sum = sum + SQUARE(phase_screen[ifourn*k1] - phase_screen[ifourn*k2]); 
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
int  FSTRUCTURE_CPLX(float *phase_screen, int *nx, int *ny, int *idim, float *fstruct)
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
/* Makes the complex phase_screen periodic with period (nx,ny): */
           jj1 = (jj1 >= ny1) ? jj1-ny1 : jj1;
           jj1 = (jj1 < 0) ? jj1+ny1 : jj1;
           kj2 = jj1 * (*idim);

           for(i1=0; i1<nx1; i1++)
             {
              k1 = kj1 + i1;
              ii1 = i1 + i;
/* Makes the complex phase_screen periodic with period (nx,ny): */
              ii1 = (ii1 >= nx1) ? ii1-nx1 : ii1;
              ii1 = (ii1 < 0) ? ii1+nx1 : ii1;
              k2 = kj2 + ii1;
/* Structure term for (i1,j1) and (i1+i,j1+j) : */
              sum = sum + SQUARE(phase_screen[2*k1] - phase_screen[2*k2]) 
                        + SQUARE(phase_screen[2*k1+1] - phase_screen[2*k2+1]); 
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
/***********************************************************************
* Extract the central part of the input array
*
* INPUT:
*  in_aa: input array
*  in_nx, in_ny: size of input array
*  in_idim, out_idim: first dimensions of input and output arrays
*
* OUTPUT:
*  out_aa: central part of the input array
*  out_nx, out_ny: size of ouput array (i.e in_nx/2, in_ny/2)
*
***********************************************************************/
static int extract_subarray(float *in_aa, int in_nx, int in_ny, int in_idim,
                            float *out_aa, int *out_nx, int *out_ny, 
                            int out_idim)
{
int ixc, iyc;
register int i, j;

ixc = in_nx/2;
iyc = in_ny/2;

for(j = 0; j < iyc; j++)
  for(i = 0; i < ixc; i++){
    out_aa[i + j * out_idim] = in_aa[i + ixc/2 + (j + iyc/2) * in_idim];
    }

*out_nx = in_nx / 2;
*out_ny = in_ny / 2;
return(0);
}
/***********************************************************************
* Copy input array to output array
*
* INPUT:
*  in_aa: input array
*  in_nx, in_ny: size of input array
*  in_idim, out_idim: first dimensions of input and output arrays
*
* OUTPUT:
*  out_aa: copy of the input array
*  out_nx, out_ny: size of ouput array (i.e in_nx, in_ny)
*
***********************************************************************/
static int copy_array(float *in_aa, int in_nx, int in_ny, int in_idim,
                      float *out_aa, int *out_nx, int *out_ny, int out_idim)
{
register int i, j;

if(in_nx > out_idim){
  fprintf(stderr,"copy_array/Fatal error: in_nx=%d out_idim=%d!\n", 
          in_nx, out_idim);
  exit(-1);
  }
for(j = 0; j < in_nx; j++)
  for(i = 0; i < in_ny; i++){
    out_aa[i + j * out_idim] = in_aa[i + (j) * in_idim];
    }

*out_nx = in_nx;
*out_ny = in_ny;
return(0);
}
/*
#define APODI hamming 
#define APODI_NAME "Hamming" 
*/
#define APODI blackman 
#define APODI_NAME "Blackman" 
/**********************************************************************
* Check if r0 is compatible with cuts of structure function along X and Y axes
* 
* INPUT:
*  in_aa: input float array
*  width1: width of filter
*  in_idim: first dimension of input array
*  out_idim: first dimension of output array
*
* OUTPUT:
*  out_aa: output apodized array
*  filter_effic: ratio of output intensity/input intensity assuming input=1
***********************************************************************/
int APODIZE_FLOAT_ARRAY(float *in_aa, float *out_aa, int *nx, int *ny, 
                        int *in_idim,int *out_idim, float *width1, 
                        float *filter_effic)
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
            out_aa[i + j * (*out_idim)] = 0.;
          else
            {
            out_aa[i + j * (*out_idim)] = in_aa[i + j * (*in_idim)] 
                                      * work * APODI(PI * argx);
            flux += work * APODI(PI * argx);
            }
        }
     }

*filter_effic = flux / (float)((*nx)*(*ny));
printf("APODIZE_FLOAT_ARRAY/Efficiency factor=%f\n", *filter_effic);

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
/******************************************
* Compute structure function with various options 
* 
* INPUT:
*  phase_screen1: original phase screen (dimension idim) 
*  idim: first dimension of arrays phase_screen1 and phase_screen2
*  NB: idim is also the first dimension of fstruct array, i.e. fstruct(idim,*)
*  nx1, ny1: size of the phase_screen1
*  i_apodization: flag set to one if apodization is required
*  i_subarray: flag set to one if extraction of a sub-array is required
*  i_output_images: flag set to one if image files have to be output 
*  profile_outfile: name of output profile
*
* OUTPUT:
*  prof_x, prof_y: X and Y cuts of the structure function
*  nprof: number of points of prof_x, prof_y
*
* WORK arrays:
*  fstruct: used to store the structure function
*  phase_screen2: used to store the transformed version of phase_screen1,
*                 according to the options i_apodization and i_subarray 
******************************************/
static int compute_fstruct(float *phase_screen1, float *phase_screen2, 
                           float *fstruct_array, float *prof_x, float *prof_y,
                           int nx1, int ny1, int idim, int *nprof, 
                           int i_apodization, int i_subarray, 
                           int i_output_images, char *psc_infile, 
                           char *psc_outfile, char *fstruct_outfile,
                           char *profile_outfile)
{
float width1, filter_effic, r0_theor_pix, r0_x, r0_y;
INT4 nx2, ny2;
int fourn_format;
char comments[120];
register int i, j;

/* Extract the central part of the array if needed: */
if(i_subarray){
    extract_subarray(phase_screen1, nx1, ny1, idim, 
                     phase_screen2, &nx2, &ny2, idim);
    printf("Sub-frame extracted from input phase screen \n");
  } else {
    copy_array(phase_screen1, nx1, ny1, idim,  
               phase_screen2, &nx2, &ny2, idim);
  }

/* Set width value for apodization: */
    width1 = 1.1 * (float)nx2;

/* Then apodize the phase array, if needed: */
  if(i_apodization) {
    APODIZE_FLOAT_ARRAY(phase_screen2, phase_screen2, &nx2, &ny2, &idim, 
                        &idim, &width1, &filter_effic);
    printf("Apodization performed on phase screen\n");
  }

/* Compute the structure function of the real phase array: 
*/
   fourn_format = 0;
/* Two (identical) ways of computing the structure function:
* either directly (but very slow):
   FSTRUCTURE_REAL(phase_screen, &nx, &ny, &nx, fstruct_array, &fourn_format);

or with Fourier Transform (fast):
(I checked with a mas 0,1 that the average value as about 0.5 which is OK)
*/
   FSTRUCTURE_REAL_WITH_FT(phase_screen2, &nx2, &ny2, &idim, fstruct_array, 
                          &fourn_format);

/* Correct the structure function taking into account the filter factor: */
  if(i_apodization){
    for(j = 0; j < ny2; j++)
      for(i = 0; i < nx2; i++)
        fstruct_array[i + j * idim] /= filter_effic;
  }

  *nprof = MINI(nx2/2, ny2/2);
  profiles_fstruct(fstruct_array, nx2, ny2, idim, prof_x, prof_y, *nprof);

/* Compute r_0 from structure function: 
*/
/*
* JLP2008: I use r0=10cm, fov=9" et 256x256 pixels in simu_close_binary.input...
*/
  r0_theor_pix = 6.7128; 
  r0_from_fstructure(prof_x, prof_y, *nprof,  &r0_x, &r0_y, r0_theor_pix,
                     profile_outfile);

if(i_output_images){
  if (i_apodization && !i_subarray)
   sprintf(comments,"apodized frame from %s", psc_infile);
  else if (!i_apodization && i_subarray)
   sprintf(comments,"subarray from %s", psc_infile);
  else
   sprintf(comments,"Phase screen copied from %s", psc_infile);

/* Save apodized/extracted frame: */
   printf(" Output transformed phase screen to %s \n", psc_outfile);
   printf("Comments: %s\n", comments);
   JLP_WRITEIMAG(phase_screen2, &nx2, &ny2, &idim, psc_outfile, comments);

/* Save structure function to file: */
   printf(" Output of structure function of the phase to %s \n", 
          fstruct_outfile);
   sprintf(comments,"Structure function of %s", psc_infile);
   JLP_WRITEIMAG(fstruct_array, &nx2, &ny2, &idim, fstruct_outfile, comments);
}

return(0);
}
