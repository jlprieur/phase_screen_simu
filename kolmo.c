/*************************************************************
*  kolmo.c 
*  Set of routines to simulate a Kolmogorov turbulence screens
*  called by simu3.c
*
*  Used as a main program, it generates simulated speckle data 
*  (IDIMxIDIM array).
*
*  JLP
*  Version 15-05-2008
************************************************************/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>         /* For file handling */
#include <jlp_ftoc.h>
#include <kolmo.h>         /* Prototypes of routines defined in kolmo.c */
#include <fstructure.h>    /* Prototypes of routines defined in fstructure.c */

/* In "poidev.c": */
#define POIDEV poidev_
float POIDEV(float xm, int *idum);
 
/*
#define DEBUG
*/
/* If Poisson random generator:
#define POISSON
*/

/* Use fourn1 Fortran routine contained in fft/jlp_fft.for
* or  fourn C routine in fft/fourn.c */
/*
#ifdef ibm
#define FOURN1 fourn1
#else
#define FOURN1 fourn1_
#endif
*/

#define IDIM 128 
#define TWOPI 6.2831853 

/* L0 is the outer scale in meters, L0 = N * Delta_L
   Coordinates in pupil plane are proportional to L0/N.
     Their range is  1 * (L0/N)   to     N * (L0/N) = L0   
     They are measured in meters 
   Coordinates in plane conjugate to pupil are proportional to 1/L0   
     Their range is  -(N/2) * 1/L0   to    (N/2 -1) * (1/L0)
     They are measured in meters-1 
   Coordinates in image plane are proportional to lambda/L0  in radians
    or (lambda/L0) * 360. * 3600. / TWOPI  in arcseconds
     Their range is  1 * lambda/L0 * constant  to  N * lambda/L0 * constant 
     They are measured in arcseconds */ 
/************************************************************/

#ifdef MAIN_PROGRAM 
int main(int argc, char *argv[])  
{
/* Data storage arrays */
float ree[2*IDIM*IDIM];
/* For autocorrelation: 
float ima[2*IDIM*IDIM];
*/

/* Packed arrays for FFT routine */
float array[2*IDIM*IDIM];

/* input variables */
float r0, L0, scale, lambda, diam, cent_obscu, norm_fourier;
int iseed, nphot, ndim, nn[2], nx, ny, idim;
char fname[61], comments[81], logfile[81];
char re_name[61], im_name[61], re_comments[81], im_comments[81];
register int i;
#ifdef FOURN1 
int ioption;
#endif
#ifdef POISSON
int nphr;
#endif
FILE *fp_log;
 
/* r0 is the Fried parameter in meters */
r0=0.1;
if(argc > 1) sscanf(argv[1],"%f",&r0);

/* Seed is needed for random generator: */
iseed = 100;
if(argc > 2) sscanf(argv[2],"%d",&iseed);

/* Mean number of photons per frame */ 
nphot = 1000;
if(argc > 3) sscanf(argv[3],"%d",&nphot);

if( argc < 3) {
  printf(" Syntax is :    kolmo [r_0 iseed nphot] \n"); 
  printf(" r_0: in meters\n iseed: integer used for random generator \n");
  printf(" nphot: number of photons/frame\n");
  printf(" Example:    kolmo 0.1 100 nphot] \n"); 
  }

/* *********************************************************** */ 
/* JLP format: */
JLP_INQUIFMT();
 
nx = IDIM; ny = IDIM;
idim = IDIM;

strcpy(logfile, "kolmo_simu.dat");
if((fp_log = fopen(logfile, "w")) == NULL) {
  fprintf(stderr, "Fatal error opening logfile: %s\n", logfile);
  return(-1);
  }
/* Lambda is in nm */
lambda = 550.;

/* diam is the telescope diameter in meters */
diam = 2.0; 
L0=2*diam;

/* cent_obscu is the central obscuration, i.e. ratio of secondary
mirror diameter with main mirror diameter */
cent_obscu = 0.2;

 printf(" Parameters of the simulation: \n");
 fprintf(fp_log, " Parameters of the simulation: \n");
 printf(" Diameter = %.3f m  r0 = %.3f m lambda = %.2f nm\n",
        diam, r0, lambda);
 fprintf(fp_log, " Diameter = %.3f m r0 = %.3f m  lambda = %.2f nm\n",
         diam, r0, lambda);
 printf(" Central obscuration = %f \n",cent_obscu); 
 fprintf(fp_log, " Central obscuration = %f \n",cent_obscu); 
 printf(" Random generator seed = %d \n",iseed); 
 fprintf(fp_log, " Random generator seed = %d \n",iseed); 
 printf(" Mean number of photons per frame = %d \n",nphot); 
 fprintf(fp_log, " Mean number of photons per frame = %d \n",nphot); 

/* scale is in arcseconds per pixel */
 scale = (lambda * 1.E-9 / L0) * 360. * 3600. / TWOPI;
 printf(" Scale in image plane is %f arcsecond/pixel \n",scale);
 fprintf(fp_log, " Scale in image plane is %f arcsecond/pixel \n",scale);

/************************************************************ */ 
JLP_RANDOM_INIT(&iseed);
nn[0] = nx; nn[1] = ny;
ndim = 2;

/************************************************************ */ 

/* First generate complex Gaussian random numbers (sigma=1 mean=0)
* in Fourier domain: */
 get_gaussft(array,nx,ny,idim);

/* Multiply with square root of Kolmogorov spectrum */ 
/* Zero frequency is in the center */
 get_phift(array,nx,ny,idim,L0,r0);

/*************************************************************/

/* Shifting center of array (to have zero frequency at (0,0) */
recent_complex(array,nx,ny,idim);

/* Inverse FFT to compute the phase term: */
#ifdef FOURN1
   ioption=-1; FOURN1(&array[0],&nn[0],&ndim,&ioption);
#else
   fourn(array-1,nn-1,2,-1); 
#endif

/* WARNING!!!!! 
 FOURN1 (and fourn) divide by nx*ny when ioption=-1 (inverse FFT).
*/
norm_fourier = sqrt((double)(nx * ny));
for(i=0; i<nx*ny*2; i++) array[i] *= norm_fourier;

/*************************************************************/
/* Autocorrelation of argument Re(phi) of pupil exp(i Re(phi)) 
*/
   autocor_real(array,nx,ny,idim,ree);

/* Output the autocorrelation of the phase of the pupil: */
   strcpy(fname,"autoc.fits");
   printf(" Output of autocorrelation (of the phase) in %s \n",fname);
   fprintf(fp_log, " Output of autocorrelation (of the phase) in %s \n",fname);
   for(i = 0; i < 80; i++) comments[i]=' ';
   comments[80]='\0';
   strcpy(comments," Autocorrelation of the phase of the pupil //");
   JLP_WRITEIMAG(ree,&nx,&ny,&idim,fname,comments);

/*************************************************************/
/* Structure function of argument Re(phi) of pupil exp(i Re(phi)) 
*/
   structure_real(array,&nx,&ny,&idim,ree);

   strcpy(fname,"struct_phase.fits");
   printf(" Output of structure function of the phase in %s \n",fname);
   fprintf(fp_log, " Output of structure function of the phase in %s \n",fname);
   for(i = 0; i < 80; i++) comments[i]=' ';
   comments[80]='\0';
   strcpy(comments," Structure function of the phase //");
   JLP_WRITEIMAG(ree,&nx,&ny,&idim,fname,comments);

/*************************************************************/
/* We only keep the real part, and compute cos(phi) and sin(phi): */
   get_phasor(array,nx,ny,idim);

/* Multiply with pupil: */
   diam = (float)(nx/2);
   mask_pupil(array,nx,ny,idim,diam,cent_obscu);

/*************************************************************/
/* Autocorrelation of pupil exp(i phi)
   autocor_cplx(array,nx,ny,idim,ree,ima);
*/

/* Structure function of pupil exp(i Re(phi)) 
*/
   structure_cplx(array,&nx,&ny,&idim,ree);

   strcpy(fname,"struct_phasor.fits");
   printf(" Output of structure function of the pupil in %s \n",fname);
   fprintf(fp_log, " Output of structure function of the pupil in %s \n",fname);
   for(i = 0; i < 80; i++) comments[i]=' ';
   comments[80]='\0';
   strcpy(comments," Structure function of the pupil //");
   JLP_WRITEIMAG(ree,&nx,&ny,&idim,fname,comments);

/*************************************************************/
/* Output pupil: */

  for(i = 0; i < 80; i++) re_comments[i]=' ';
  re_comments[80]='\0';
  strcpy(re_name,"phare.fits");
  strcpy(re_comments," Phase -real- //");

  for(i = 0; i < 80; i++) im_comments[i]=' ';
  im_comments[80]='\0';
  strcpy(im_name,"phaim.fits");
  strcpy(im_comments," Phase -imaginary- //");

  output_pupil(array,nx,ny,idim,re_name,im_name,re_comments,im_comments);
/*************************************************************/

/* FFT to go from pupil to image plane: */
#ifdef FOURN1
   ioption=1; FOURN1(&array[0],&nn[0],&ndim,&ioption);
#else
   fourn(array-1,nn-1,2,1); 
#endif

/* From field to intensity, and normalisation to nphot */
   get_intensity(array,ree,nx,ny,idim,nphot); 

/* Poisson noise to simulate photo-counting detector
*/
#ifdef POISSON
   kolmo_clipping(ree,nx,ny,idim,&nphr,&iseed);
   printf(" Final number of photons %d \n",nphr);
   fprintf(fp_log, " Final number of photons %d \n",nphr);
#else
   printf("No Poisson noise generator.\n");
   fprintf(fp_log, "No Poisson noise generator.\n");
#endif

/* Shifting center of array (To have the patch in the middle)*/
   recent_real(ree,nx,ny,idim);

/*************************************************************/
/* Output of intensity: */
  strcpy(fname,"intens1.fits");
  printf(" Output of intensity in %s \n",fname);
  fprintf(fp_log, " Output of intensity in %s \n",fname);
/* add a '\0' to prevent any accident... */
   for(i = 0; i < 80; i++) comments[i]=' ';
   comments[80]='\0';
  strcpy(comments," Intensity in image plane//");
  JLP_WRITEIMAG(ree,&nx,&ny,&idim,fname,comments);
/*************************************************************/

  printf("See details in logfile: %s\n", logfile);
  fclose(fp_log);

return(0);
}
/* End of set of instructions when MAIN_PROGRAM is defined*/
#endif
/*********************************************************
* Generate random numbers following
* a Gaussian law with sigma=1. and mean=0. 
*
*
* OUTPUT:
*  array[]: Gaussian random numbers stored as real, imag, real, imag...
*           ("fourn" format but starting at index=0)
**********************************************************/
int get_gaussft(float *array, int nx, int ny, int idim)
{
register int i, j, kj, k;

 for(j = 0; j < ny; j++)
  {
   kj = j * idim;
   for(i = 0; i < nx; i++)
    {
      k = 2 * (kj + i);
/*
* Complex Gaussian random numbers are generated by the inverse
* transform method.
float gmod, arg, rad;
      gmod = 0.; 
      while(gmod <= 0.) JLP_RANDOM(&gmod);
      rad = sqrt(-2.0 * log((double)gmod)); 
      JLP_RANDOM(&arg);
      arg = TWOPI * arg;
      array[k] = rad * cos((double)arg); 
      array[k+1] = rad * sin((double)arg); 
*/
      JLP_RANDOM_GAUSS(&array[k]);
      JLP_RANDOM_GAUSS(&array[k+1]);
    }
  }

return(0);
}
/*********************************************************
* Multiply with square root of Kolmogorov spectrum 
*    Phi(r) =   0.023 * nx * ny * (L0/r0)^5/3 * r^{11/3}
*
* INPUT:
*   array: Fourier data (real, imag, real, imag, ...)
*          (Gaussian random numbers)
*   L0: external scale in meters
*   r0: Fried radius in meters
*
* OUTPUT:
*   array: Fourier data multiplied by Kolmo spectrum 
**********************************************************/
int get_phift(float *array, int nx, int ny, int idim, float L0, float r0)
{
double rad2;
float wkol, wkol0;
register int i, j, kj, k, icent, jcent, ii, jj;

/* Coordinates of center: */
icent = nx/2; 
jcent = ny/2; 

/* 1.666 is 5/3 */
wkol0 = 0.023 * nx * ny * pow((double)(L0/r0),1.6666);
wkol0 = sqrt(wkol0);

 for(j=0; j<ny; j++)
  {
   kj=j*idim;
   jj = j-jcent;
   for(i=0; i<nx; i++)
    {
      ii = i-icent;
      rad2 = (double)(SQUARE(ii) + SQUARE(jj));
/* To avoid division by zero and too large value at 0,0: */
      if(rad2 == 0.) rad2 = 0.5;
/* wkol: sqrt of Kolmogorov spectrum at frequency (ii,jj) */
/* 1.83333 is 11/6 */
      wkol = sqrt( wkol0 /pow(rad2,1.83333));
      k = 2*(kj + i);
      array[k] = array[k] * wkol; 
      array[k+1] = array[k+1] * wkol; 
    }
  }

return(0);
}
/*********************************************************
* Compute phasor term by keeping only real part of phase. 
*
**********************************************************/
int get_phasor(float *array, int nx, int ny, int idim)
{
double arg;
register int i, j, kj, k;

 for(j=0; j<ny; j++)
  {
   kj=j*idim;
   for(i=0; i<nx; i++)
    {
      k = 2*(kj + i);
      arg = (double)array[k];
      array[k] = cos(arg);
      array[k+1] = sin(arg);
    }
  }

return(0);
}
/***************************r******************************
* Shifting center of real array (for FFT) 
*
**********************************************************/
int recent_real(float *array, int nx, int ny, int idim)
{
int i, j, icent, jcent, kj, kjcent;
float w1, w2;

icent = nx/2;
jcent = ny/2;
kjcent = jcent * idim;

   for(j=0; j<jcent; j++)   
   {
    kj = j*idim;
    for(i=0; i<icent; i++)   
      {
       w1=array[i+kj];
       w2=array[i+icent+kj];
       array[i+kj]=array[i+icent+kj+kjcent];
       array[i+icent+kj]=array[i+kj+kjcent];
       array[i+icent+kj+kjcent]=w1;
       array[i+kj+kjcent]=w2;
      }
   }    

return(0);
}
/***************************r******************************
* Shifting center of complex array (for FFT) 
*
**********************************************************/
int recent_complex(float *array, int nx, int ny, int idim)
{
int i, j, icent, jcent, kj, kjcent;
float w1, w2;
float z1, z2;

icent = nx/2;
jcent = ny/2;
kjcent = jcent * idim;

   for(j=0; j<jcent; j++)   
   {
    kj = j*idim;
    for(i=0; i<icent; i++)   
      {
/* Real part: */
       w1=array[2*(i+kj)];
       w2=array[2*(i+icent+kj)];
       array[2*(i+kj)]=array[2*(i+icent+kj+kjcent)];
       array[2*(i+icent+kj)]=array[2*(i+kj+kjcent)];
       array[2*(i+icent+kj+kjcent)]=w1;
       array[2*(i+kj+kjcent)]=w2;
/* Imaginary part: */
       z1=array[2*(i+kj)+1];
       z2=array[2*(i+icent+kj)+1];
       array[2*(i+kj)+1]=array[2*(i+icent+kj+kjcent)+1];
       array[2*(i+icent+kj)+1]=array[2*(i+kj+kjcent)+1];
       array[2*(i+icent+kj+kjcent)+1]=z1;
       array[2*(i+kj+kjcent)+1]=z2;
      }
   }    

return(0);
}
/***************************r******************************
* Pupil mask 
* Assume that the telescope diameter is half of the pupil...
*
**********************************************************/
int mask_pupil(float *array, int nx, int ny, int idim, float diam, 
               float cent_obscu)
{
int icent, jcent, ii, jj, i, j, k, kj;
float rad2, radmin2, radmax2;

icent = nx/2;
jcent = ny/2;
radmax2 = (float)nx/4;
radmin2 = radmax2 * cent_obscu;
radmax2 = radmax2 * radmax2;
radmin2 = radmin2 * radmin2;

   for(j=0; j<ny; j++)   
   {
    kj = j*idim;
    jj = j-jcent;
    for(i=0; i<nx; i++)
     {
      ii = i-icent;
      rad2 = (float)(ii*ii + jj*jj);
      if(rad2 > radmax2 || rad2 < radmin2) 
         {
          k = 2*(kj + i);
          array[k] = 0.; 
          array[k+1] = 0.; 
         } 
     }
   }
return(0);
}
/***************************r******************************
* Output pupil 
*
**********************************************************/
int output_pupil(float *array, int nx, int ny, int idim, char *re_name, 
                 char *im_name, char *re_comments, char *im_comments)
{
float *ree, *ima;
int isize;
register int i, j, k, kj;
 
isize = nx * ny * sizeof(float);

JLP_GVM(&ree, &isize);
JLP_GVM(&ima, &isize);

/* Unpacking : */ 
     for(j=0; j<ny; j++)  
     {
     kj = j*idim;
     for(i=0; i<nx; i++)  
        {
        k = kj + i;
        ree[k]=array[2*k];
        ima[k]=array[2*k+1];
        }
     }
 

/* Create output files */
 printf(" Output of %s \n",re_name);
 JLP_WRITEIMAG(ree,&nx,&ny,&idim,re_name,re_comments);

 printf(" Output of %s \n",im_name);
 JLP_WRITEIMAG(ima,&nx,&ny,&idim,im_name,im_comments);

 JLP_FVM(&ree);
 JLP_FVM(&ima);

return(0);
}
/***************************r******************************
* Compute intensity from array to intens 
*
**********************************************************/
int get_intensity(float *array, float *intens, int nx, int ny, int idim, 
                  int nphot) 
{
double sum1;
int status, i, j, k, kj;

status = 0;

   for(j=0; j<ny; j++)   
   {
    kj = j*idim;
    for(i=0; i<nx; i++)
     {
       k = kj + i;
       intens[k] = array[2*k] * array[2*k] + array[2*k+1] * array[2*k+1];
       sum1 = sum1 + intens[k];
     }
   }

/* Normalisation to nphot: */
if(sum1 != 0.)
   {
   sum1 = (float)nphot / sum1;
   for(j=0; j<ny; j++)   
   {
    kj = j*idim;
    for(i=0; i<nx; i++)
     {
       k = kj + i;
       intens[k] = intens[k] * sum1; 
     }
   }
   }
else
   {
   printf(" get_intensity/ error: number of photons is zero !\n");
   status = -1;
   }

return(status);
}
/***************************r******************************
*
* Compute autocorrelation of a complex array
*
**********************************************************/
int  autocor_cplx(float *array, int nx, int ny, int idim, float *ree, 
                  float *ima)
{
int i, j, k, kj;
int i1, j1, k1, kj1, ii1, jj1;
int k2, kj2;
double sumr, sumi;

/***** main loop: ****/
   for(j=-ny/2; j<ny/2; j++)   
   {
    kj = (j + ny/2) * idim;
    for(i=-nx/2; i<nx/2; i++)
     {
       k = kj + i + nx/2;
       sumr = 0.;
       sumi = 0.;
/**** Internal loop: */
         for(j1=0; j1<ny; j1++)   
          {
           kj1 = j1 * idim;
           jj1 = j1 + j;
/* Makes the complex array periodic with period (nx,ny): */
           jj1 = (jj1 >= ny) ? jj1-ny : jj1;
           jj1 = (jj1 < 0) ? jj1+ny : jj1;
           kj2 = jj1 * idim;

           for(i1=0; i1<nx; i1++)
             {
              k1 = kj1 + i1;
              ii1 = i1 + i;
/* Makes the complex array periodic with period (nx,ny): */
              ii1 = (ii1 >= nx) ? ii1-nx : ii1;
              ii1 = (ii1 < 0) ? ii1+nx : ii1;
              k2 = kj2 + ii1;
/* Real part of f*(i1,j1) * f(i1+i,j1+j) : */
              sumr = sumr + array[2*k1]*array[2*k2] 
                        + array[2*k1+1]*array[2*k2+1];
/* Imaginary part of f*(i1,j1) * f(i1+i,j1+j) : */
              sumi = sumi - array[2*k1+1]*array[2*k2] 
                        + array[2*k1]*array[2*k2+1];
             }
           }
       ree[k] = (float)sumr; 
       ima[k] = (float)sumi; 
     }
   }
return(0);
}
/***************************r******************************
*
* Compute autocorrelation of the real part of a complex array
*
**********************************************************/
int  autocor_real(float *array, int nx, int ny, int idim, float *ree)
{
int i, j, k, kj;
int i1, j1, k1, kj1, ii1, jj1;
int k2, kj2;
double sumr;

/***** main loop: ****/
   for(j = -ny/2; j < ny/2; j++)   
   {
    kj = (j + ny/2) * idim;
    for(i = -nx/2; i < nx/2; i++)
     {
       k = kj + i + nx/2;
       sumr = 0.;
/**** Internal loop: */
         for(j1 = 0; j1 < ny; j1++)   
          {
           kj1 = j1 * idim;
           jj1 = j1 + j;
/* Makes the complex array periodic with period (nx,ny): */
           jj1 = (jj1 >= ny) ? jj1-ny : jj1;
           jj1 = (jj1 < 0) ? jj1+ny : jj1;
           kj2 = jj1 * idim;

           for(i1 = 0; i1 < nx; i1++)
             {
              k1 = kj1 + i1;
              ii1 = i1 + i;
/* Makes the complex array periodic with period (nx,ny): */
              ii1 = (ii1 >= nx) ? ii1-nx : ii1;
              ii1 = (ii1 < 0) ? ii1+nx : ii1;
              k2 = kj2 + ii1;
/* Real part of f*(i1,j1) * f(i1+i,j1+j) : */
              sumr = sumr + array[2*k1]*array[2*k2]; 
             }
           }
       ree[k] = (float)sumr; 
     }
   }

return(0);
}
/*************************************************************
* Clipping routine to generate simulated photo-counting "frame"
* with Poisson noise.
* Calls poidev() to generate Poisson noise,
* from Numerical recipees Chap 7, in "poidev.c".
*
* Before calling this routine, 
* the total flux should be normalised to nphot!
*
*************************************************************/
int kolmo_clipping(float *intens, int nx, int ny, int idim, int *nphr,
                   int *iseed)
{
int i, j, k, kj;

   *nphr=0;
   for(j=0; j<ny; j++)   
   {
    kj = j*idim;
    for(i=0; i<nx; i++)
     {
       k = kj + i;
/* Poisson random generator */
       if(intens[k] > 0.) 
                intens[k] = (int)POIDEV(intens[k],iseed);
/* Computing total number of photons: */
       *nphr += intens[k];
     }
   }

return(0);
}
