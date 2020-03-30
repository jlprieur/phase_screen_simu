/******************************************************************************
* s_simu1.c
*
* from s_simu.for version 26/07/2008) 
* (from s_pscorg.for version 2004) is a phase screen simulator which produces 
* artificial fringe
* or speckle patterns and obtain integrated power spectrum or(/and)
* optionally specified bispectral components.
*
* The atmospheric disturbance is given by the Kolmogorov spectrum.
* Seeing condition, aperture configuration and the object ( at this
* point, only multiple point sources ) and the usage of the phase screen
* must be specified by the user in an ascii input file.
*
* Example:
*   runs s_simu simu_triple.input triple phot_gv7_128.fits
*
* Following is a sample input file. (free format, ! indicates a comment line.)
*
* ! Size of images
* Size
* 128
* ! mean number of photons per frame
* Photons
* 200
* ! Field of view: side of image in arcseconds
* Field of view
* 12
* ! seed for random number generator (large odd integer)
* Random
* 7597551
* ! outer scale, inner scale of the turbulence(m) - recommended values
* Turbulence scale
* 10.0  0.01
* ! Fried radius (r_zero) in m 
* Fried 
* 0.1
* ! central wavelength(nm), bandwidth (nm), sampling points
* Wavelength
* 650 70 3
* ! source: # of point sources
* !         intensity(arbitrary unit), position offset(x,y) in arcsecond
* Source
* 2
* 1.0  -0.7 -0.6
* 0.5   0.8 0.4
* ! # of circular apertures (including voids)
* ! configuraton: center of aperture(x,y), diameter(m) and type
* ! 1 is an aperture,0 is a void. Voids should follow ordinary apertures.
* ! Example: annular mask
* Aperture
* 2
* 0.0 0.0 2.0  1
* 0.0 0.0 1.9  0
* ! total number of frames to be processed
* Frames
* 100
* ------ as an option `no atmosphere' can be specified by ------
* ! obtaining beam pattern
* No_disturbance
* --------------------------------------------------------------
*
* JLP
* Version 21/07/2008
*******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <fcntl.h>        /* For file handling */
#include <jlp_ftoc.h>
#include <jlp_fftw.h>     /* For FFTs */
#include <../hrsa/jlp_cover_mask.h> /* For COVERA and BISP1 */
#include <fstructure.h>   /* Prototypes of routines defined in fstructure.c */
#include <jlp_complex.h>  /* Complex routines */
#include <jlp_poidev.h>   /* Poisson routines */

#ifndef PI
#define PI 3.141592653
#endif
/* rapsec: radians per arcsecond */
#define RAD_PER_ARCSECOND (PI/(180. * 3600.0))

/* MAXCOMP: maximum number of components in source 
* MAXAP: maximum number of apertures
*/
#define MAXCOMP 50
#define MAXAP 15

#define DEBUG

void fourn_for_C_arrays(float *data, int *nn, int ndim, int isign);

/**************************************************************************
* Contained here:
**************************************************************************/
int JLP_RECENT_CMPLX(complex *psc, int Ns, int idim);
int JLP_RECENT_FLOAT(float *psc, int Ns, int idim);
int JLP_RECENT_DOUBLE(double *psc, int Ns, int idim);
static int read_param_file(char *parameter_fname, int *Np, long *seed, int *nph,
                           float *L0, float *lo, float *lambda_cent, 
                           float *bandwidth, int *nsamp, int *ncomp, 
                           float *s_int, float s_pos[][2], int *ncircle, 
                           float c_pos[][2], float *c_rad, int *c_typ, 
                           float *fov, float *r_zero, int *nframes, 
                           int *no_dist, char *out_prefix);
static int read_photon_response(char *photon_modsq_fname, int Np, 
                                complex **DetectorTransf);
static int compute_scale_and_resolution(float lambda_cent, float r_zero,
                                        int Np, int Ns, float fov, float *LCn2, 
                                        float *Ls, float *Scale_s, 
                                        float *image_scale);
static int create_object(float **object, complex **FT_object, int Np, 
                         float fov, int ncomp, float *s_int, float s_pos[][2], 
                         char *out_prefix);
static int compute_sampling(float Ls, int Ns, int Np, 
                            int *frm_per_screen_width, 
                            float *x_min, float *x_max, float *xstep, 
                            float *y_min, float *y_max, float *ystep, 
                            int ncircle, float c_pos[][2], float *c_rad);
static int init_arrays(complex **psc, double **pwr, double **modsq, 
                       double **snrm, double **long_int, double **bisp1, 
                       int Ns, int Np, int ngamma);
static int free_arrays(complex *psc, double *pwr, double *modsq, 
                       double *snrm, double *long_int, double *bisp1);
static int jlp_phase_screen(complex *psc, int Ns, float Ls, 
                            float r_zero, int no_dist, char *out_prefix, 
                            int *output_demo);
static int process_frames(complex *psc, int Ns, complex *FT_object, 
                          complex *DetectorTransf, double *pwr, double *modsq, 
                          double *long_int, double *snrm, double *bisp1, int Np,
                          int *iframe, int nframes, float Ls, long *seed, 
                          int nph, int ir, int nbeta, int ngamma, float x_min, 
                          float x_max, float xstep, float y_min, float y_max, 
                          float ystep, float lambda_cent, float bandwidth, 
                          int nsamp, int ncomp, float *s_int, float s_pos[][2],
                          int ncircle, float c_pos[][2], float *c_rad, 
                          int *c_typ, char *out_prefix, int *output_demo);
static int compute_elementary_PSFs(complex *psc, float xx, float yy, float Ls,
                                   float lambda_cent, float bandwidth, 
                                   int nsamp, int ncircle, float c_pos[][2], 
                                   float *c_rad, int *c_typ, 
                                   complex *PupilFunc1, complex *PupilFunc2,
                                   float *PSF1, float *PSF2, int Ns, int Np, 
                                   int *output_demo, char *out_prefix);
static int process_elementary_frame(double *pwr, double *modsq, 
                               double *long_int, double *snrm, double *bisp1, 
                               float *PSF1, complex *FT_object, 
                               complex *DetectorTransf, int Np, long *seed, 
                               int nph, int ir, int nbeta, int ngamma, 
                               int *output_demo, char *out_prefix);
static int photon_clipping(float *in_array, float *out_array, int Np, int nph, 
                           int *nphr, long *seed);
static int output_results(double *pwr, double *modsq, double *long_int, 
                          double *snrm, double *bisp1, int Np, int ngamma, 
                          int nframes, char *parameter_fname, char *out_prefix);

/**************************************************************************/
int main(int argc, char *argv[])
{
/* Ns = size of phase screen 
*  Np = size of images */
int Ns, Np;
/* Maximum size for Np and Ns: */
int idim, idim2;

/*
* Integrated structure constant of the phase fluctuation
* LCn2 = \int Cn2(h) dh
* And also: 
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
*/
float LCn2, r_zero;

/* Ls: width of phase screen, Scale_s: resolution of phase screen (m/pix)
*/
float Ls, Scale_s;
 
/* Wave number corresponding to lambda_cent, scale length of turbulence
* fov: field of view (for the simulated images)
*/
float L0, lambda_cent, lo, fov;
 
/* psc: complex phase screen, real and imag are independent.
* FT_object(Np,Np): Fourier Transform of the object
*/
complex *psc, *FT_object;

/* Detector transfer function: */
complex *DetectorTransf;
 
/* pwr: mean power spectrum of analog elementary frames 
* modsq: mean power spectrum with photon clipping and detector photon response
* long_int: long integration with photon clipping and detector photon response 
*/
double *pwr, *modsq, *long_int, *snrm, *bisp1;
 
/* nframes: number of frames
* frm_per_screen_width: number of frames per phase screen in X or in Y
* nph: mean number of photons per frame
*/
int  nframes, frm_per_screen_width, nph;
 
/* Source information (intensity and position in arcsec) */
float s_int[MAXCOMP], s_pos[MAXCOMP][2], *object;
 
/* Aperture information (coordinate, radius, type)
*/
float c_pos[MAXAP][2], c_rad[MAXAP];
int c_typ[MAXAP];

/* File names: */
char parameter_fname[60], photon_modsq_fname[60], psc_fname[60], out_prefix[20];

/* Miscellaneous: */
float x_min, x_max, xstep, y_min, y_max, ystep;
float bandwidth, image_scale;
int ncircle, nsamp, iframe;
int ncomp, no_dist, output_demo = 0;
long seed;
int ir, max_nclosure, nbeta, ngamma;
char *pc;
register int k;

if(argc == 7) {
 if(*argv[4]) argc = 5;
 else if(*argv[3]) argc = 4;
 else if(*argv[2]) argc = 3;
 else if(*argv[1]) argc = 2;
}
/* Example:
*   s_simu1 simu_triple.input triple phot_gv7_128.fits
*/
if(argc != 4) {
  fprintf(stderr,"argc=%d\n", argc);
  fprintf(stderr,"fstructure/Fatal error\n");
  fprintf(stderr," Syntax: s_simu1 parameter_file out_prefix photon_modsq\n");
  fprintf(stderr," (photon_modsq = Power spectrum of detector photon response)\n");
  fprintf(stderr," Example: s_simu1 simu_triple.input triple phot_gv7_128.fits\n");
  return(-1);
  }

/* Read input parameters: */
strcpy(parameter_fname, argv[1]);
strcpy(out_prefix, argv[2]);
/* Removes blanks at the end if necessary: */
pc = out_prefix;
while(*pc && *pc != ' ') pc++;
*pc = '\0';
strcpy(photon_modsq_fname, argv[3]);
sprintf(psc_fname, "%s_psc.fits", out_prefix);

/* Output parameters to check it is OK: */
printf("OK: photon_modsq=%s psc_outfile=%s out_prefix=>%s<\n", 
        photon_modsq_fname, psc_fname, out_prefix);
printf("Will output phase screen to: >%s<\n", psc_fname);

JLP_INQUIFMT();

/* First read the parameter file (in order to obtain the value of Np):
*/
 read_param_file(parameter_fname, &Np, &seed, &nph, &L0, &lo, 
                 &lambda_cent, &bandwidth, &nsamp, &ncomp, 
                 s_int, s_pos, &ncircle, 
                 c_pos, c_rad, c_typ, &fov, &r_zero, 
                 &nframes, &no_dist, out_prefix);
  Ns = 2 * Np;
  idim2 = Ns;
  idim = Np;

/* Then read the photon response of the detector and check that nx=Np:
*/
 read_photon_response(photon_modsq_fname, Np, &DetectorTransf);

/* Compute the uv coverage with ir=30 and nclosure=200: */
  ir = 30;
  max_nclosure = 200;
  printf(" Radius of uv-coverage (IR) in pixels: %d\n", ir);
  printf(" Maximum number of closure relations: max_nclosure=%d\n", 
           max_nclosure);
  COVERA(&ir, &max_nclosure, &nbeta, &ngamma);

/* Initialization of random generator: */
  JLP_RANDOM_INIT(&seed);

/* Compute resolution of image and phase screen: */
 compute_scale_and_resolution(lambda_cent, r_zero, Np, Ns, fov, 
                              &LCn2, &Ls, &Scale_s, &image_scale);

/* Create an image of the source from ncomp, s_int, s_pos
*/
  create_object(&object, &FT_object, Np, fov, ncomp, s_int, s_pos, out_prefix);

/* frm_per_screen_width: number of frames per phase screen in X or in Y
*/

/* Compute sampling parameters and frm_per_screen_width: 
*/
 compute_sampling(Ls, Ns, Np, &frm_per_screen_width, &x_min, &x_max, &xstep, &y_min, 
                  &y_max, &ystep, ncircle, c_pos, c_rad);

/*
* Allocate memory and initialize arrays:
*/
  init_arrays(&psc, &pwr, &modsq, &snrm, &long_int, &bisp1, Ns, Np, ngamma);

/* Two real phase screens per complex random number generation in Fourier
* domain.
*/
  printf(" ******** Will now simulate %d elementary frames ************\n", 
         nframes);
  printf(" Number of frames per phase screen: %d\n", SQUARE(frm_per_screen_width)); 
  printf(" Number of calls of phase_screen routine needed: %d\n", 
         MAXI(1,(nframes/(2 * SQUARE(frm_per_screen_width)))));
 
  iframe = 0;
/* Main loop (maximum number of iterations will be nframes)
*/
  for(k = 0; k < nframes; k++) {
/*------------------------------------
* Create a complex phase screen
*/
      printf("+ generating %dth complex phase screen (iframe=%d)\n",
              k, iframe+1);
/* Generation of a phase screen: 
*/
      jlp_phase_screen(psc, Ns, Ls, r_zero, no_dist, out_prefix, &output_demo);
 
/* Sample phase screen, form images, calculate power spectrum
*/
      process_frames(psc, Ns, FT_object, DetectorTransf, 
                    pwr, modsq, long_int, snrm, bisp1, Np, &iframe, 
                    nframes, Ls, &seed, nph, ir, nbeta, ngamma,
                    x_min, x_max, xstep, y_min, y_max, ystep,
                    lambda_cent, bandwidth, nsamp, ncomp,
                    s_int, s_pos, ncircle, c_pos, c_rad, c_typ,
                    out_prefix, &output_demo);
 
/* Exit from loop when nframes have been generated */
        if(iframe >= nframes - 1) break; 
   } 


/* Output FITS files:
*/
  output_results(pwr, modsq, long_int, snrm, bisp1, Np, ngamma, nframes, 
                 parameter_fname, out_prefix);

/* Free memory: */
free_arrays(psc, pwr, modsq, snrm, long_int, bisp1);
free(object);
free(FT_object);
free(DetectorTransf);

return(0);
}
/***************************************************************************
* Creates a complex phase screen (JLP's version)
*
* INPUT:
* Ns: size of phase screen 
* output_demo = 0 changed to 1 in exit
* r_zero: Fried's parameter (in meters) for lambda_cent
* Scale_s: Screen scale (in meter/pixel)
* Ls: size of screen in meters (corresponding to Ns pixels)
* no_dist: flag set to one if no atmospheric disturbance
* Np: size of image
* fov: field of view
*
* OUPUT:
*
***************************************************************************/
static int jlp_phase_screen(complex *psc, int Ns, float Ls, 
                            float r_zero, int no_dist, char *out_prefix, 
                            int *output_demo)
{
/* nn[]: dimension array for "fourn" FFT routine
*/
int fourn_format, ix_cent, iy_cent, nn[2];
complex iiu, u_zero;
float urad, utheta, radi, Fs, Fs0, *work, *structure_funct;
double rad2, sum, sumsq, mean, sigma;
register int i, ix, iy;
char out_file[60], out_comments[80];

/* Constants:
* iiu = sqrt(-1) or "i"
*/
  iiu = float_to_cplx(0., 1.);

/* Case of no atmospheric disturbance
*/
  u_zero = float_to_cplx(0., 0.);
  if (no_dist) { 
    for(i = 0; i < Ns * Ns; i++) psc[i] = u_zero;
    return(0);
   }

/* r_0 = (0.423 kk^2 LCn2)^{-3/5}
* Hence Fs0 is proportional to kk^2
*/
/* Formula from Porro (2000, Applied Optics 39, 10, p 1644)
*   Fs0 = sqrt(0.023) * (Ls / r_zero)**(5./6.)
*/
/* Formula derived by JLP:
*/
   Fs0 = sqrt(0.023 * Ns * Ns * pow((Ls/ r_zero),(5./3.)));
/* Will need a JLP's ad'hoc correction (see below):
   Fs0 = Fs0 * Ns
*/

/*------ generation of random phase in spatial domain -------------
* Complex Gaussian random numbers are generated by
* the inverse transform method.
*/
 
ix_cent = Ns / 2;
iy_cent = Ns / 2;
 
for(iy = 0; iy < Ns; iy++){
  for(ix = 0; ix < Ns; ix++) {
    rad2 = SQUARE(ix - ix_cent) + SQUARE(iy - iy_cent);
 
/* rad**(-11/3) is undefined at rad2=0.
*/
    if(rad2 == 0.0) rad2 = 0.5; 

/* 11/6=1.833333
*/
    Fs = pow(rad2, -11.0/6.0);
 
    while(1) {
      JLP_RANDOM(&urad);
      if(urad > 0) {
       radi = sqrt(-2. * log(urad));
       break;
      }
    }
    JLP_RANDOM(&utheta);
/* iiu = sqrt(-1) or "i"
* Warning: psc is then set with a complex number, since cexp is
* a complex exponential!
    psc[ix + iy * Ns] = Fs0 * sqrt(Fs) * radi * cexp(iiu*2.0*pi*utheta);
*/
    psc[ix + iy * Ns].re = Fs0 * sqrt(Fs) * radi * cos((double)2.0*PI*utheta);
    psc[ix + iy * Ns].im = Fs0 * sqrt(Fs) * radi * sin((double)2.0*PI*utheta);
/*
* For DEBUG (calibration)
*   psc[ix + iy * Ns] = radi*cexp(iiu*2.0*pi*utheta)
* or:
*   psc[ix + iy * Ns].re = radi * cos(2.0*PI*utheta);
*   psc[ix + iy * Ns].im = radi * sin(2.0*PI*utheta);
*/
  } /* EOF loop on ix */
} /* EOF loop on iy */
 
/*---------- Fourier transform ------------------------------
* Shift origin of Fourier space to (0,0) before inverse Fourier transform
*/
  JLP_RECENT_CMPLX(psc, Ns, Ns);

/* SHOULD USE Numerical recipee's fourn routine (in this directory) 
*/
  nn[0] = Ns;
  nn[1] = Ns;
  fourn_for_C_arrays((float*)psc, nn, 2, -1);
 
/* Correction since "fourn" has divided by Ns when computing the inverse FFT:
* (Calibrated in July 2008 with random Gauss, and obtained sigma=1 in Fourier
* domain from sigma=1 in direct domain)
*/
  for(i = 0; i < Ns * Ns; i++) {
    psc[i].re *= Ns;
    psc[i].im *= Ns;
    }

/* Now psc is a complex random phase screen
* Center array to avoid a discontinuity of phase screen at (Ns/2,Ns/2)
*/
  JLP_RECENT_CMPLX(psc, Ns, Ns);

  printf("OK6: psc 10, 20: %f %f\n", psc[9 + 19 * Ns].re, psc[9 + 19 * Ns].im);

/* Output phase screen:
*/
  if(*output_demo == 0){
     work = (float *)malloc(Ns * Ns * sizeof(float));
     structure_funct = (float *)malloc(Ns * Ns * sizeof(float));
     sum = 0.;
     sumsq = 0.;
     for(i = 0; i < Ns * Ns; i++) {
        work[i] = psc[i].re;
        sum += work[i];
        sumsq += SQUARE(work[i]);
        }
     mean = sum / (float)(Ns * Ns);
     sigma = sqrt(sumsq / (float)(Ns * Ns) - SQUARE(mean));
     printf("Phase screen: mean=%e sigma=%e \n", mean, sigma);
     sprintf(out_file, "%s_phase_screen", out_prefix);
     strcpy(out_comments, "Real part of phase screen");
     JLP_WRITEIMAG(work, &Ns, &Ns, &Ns, out_file, out_comments);

/* Structure function of the real part of the complex phase array:
*/
     fourn_format = 1;
     FSTRUCTURE_REAL_WITH_FT((float *)psc, &Ns, &Ns, &Ns, structure_funct, 
                             &fourn_format);
     sprintf(out_file, "%s_phase_struct", out_prefix);
     strcpy(out_comments, "Struct. func. of real part of cplx phase screen");
     JLP_WRITEIMAG(structure_funct, &Ns, &Ns, &Ns, out_file, out_comments);

/* Free memory: */
     free(work);
     free(structure_funct);
     *output_demo = 1;
  } /* EOF if output_demo == 0 */

return(0);
}
/*******************************************************************************
* process_frames
* 
* Produces fringe/speckle patterns of multiple point sources.
* Finite bandwidth effect is taken into account.
*
* INPUT:
* psc(idim2,idim2): phase screen
* Ns: size of phase screen
* FT_object: Fourier Transform of the source
* DetectorTransf: detector transfer function (i.e., modulus of Fourier 
*                 Tranform of PSF)
* Np: size of elementary frames (=Ns/2)
* iframe: index of current frame to be processed
* nframes: number of frames to be simulated
* ncomp: number of components for the source
* nph: mean number of photons per frame
* s_int, s_pos: source information (intensity and position in arcsec)
* float s_int[MAXCOMP], s_pos[MAXCOMP][2]
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* float c_pos[MAXAP][2], c_rad[MAXAP];
* Ls: size of phase screen (in m)
* lambda_cent, bandwidth, nsamp: wavelength parameters: central wavelength (m),
*                          bandwidth (m) and sampling points
* x_min, x_max, xstep : sampling parameters to extract a phase screen 
*                       from a larger one 
* y_min, y_max, ystep : sampling parameters to extract a phase screen 
*                       from a larger one 
* output_demo: 1 at first call, and it is incremented in exit
*
* OUTPUT:
* pwr: mean power spectrum of analog elementary frames 
* modsq: mean power spectrum with photon clipping and detector photon response
* long_int: long integration with photon clipping and detector photon response 
*******************************************************************************/
static int process_frames(complex *psc, int Ns, complex *FT_object, 
                          complex *DetectorTransf, double *pwr, double *modsq, 
                          double *long_int, double *snrm, double *bisp1, int Np,
                          int *iframe, int nframes, float Ls, long *seed, 
                          int nph, int ir, int nbeta, int ngamma, float x_min, 
                          float x_max, float xstep, float y_min, float y_max, 
                          float ystep, float lambda_cent, float bandwidth, 
                          int nsamp, int ncomp, float *s_int, float s_pos[][2],
                          int ncircle, float c_pos[][2], float *c_rad, 
                          int *c_typ, char *out_prefix, int *output_demo)
{
complex iiu;
complex *PupilFunc1, *PupilFunc2;
float *PSF1, *PSF2, xx, yy, Scale_s;

PSF1 = (float *)malloc(Np * Np * sizeof(float));
PSF2 = (float *)malloc(Np * Np * sizeof(float));
PupilFunc1 = (complex *)malloc(Np * Np * sizeof(complex));
PupilFunc2 = (complex *)malloc(Np * Np * sizeof(complex));
Scale_s = Ls / Ns;

/* iiu = sqrt(-1) or "i"
*/
   iiu = float_to_cplx(0., 1.);

/* Setting apertures on Np sized window in Ns sized complex phase screen
*/
 
for(xx = x_min; xx <= x_max; xx += xstep){
  for(yy = y_min; yy <= y_max; yy += ystep){
/**************************** major loop ***************************************
*/
#ifdef DEBUG
    printf("+ creating a pair of speckle/Centre of subwindow at x=%.3f y=%.3f (m) (or %d %d pixels)\n",
            xx, yy, (int)(xx/Scale_s), (int)(yy/Scale_s)); 
#endif
 
 
/* Compute two elementary PSFs (PSF1 and PSF2)
*/
    compute_elementary_PSFs(psc, xx, yy, Ls, lambda_cent, bandwidth, nsamp, 
                            ncircle, c_pos, c_rad, c_typ, PupilFunc1, 
                            PupilFunc2, PSF1, PSF2, Ns, Np, output_demo, 
                            out_prefix);

/* Convolution of the object with the elementary PSFs 
*/
        
/* Processing a first frame with real part of phase screen:
*/
    if(*iframe < 10 || (*iframe % 20) == 0) 
            printf("+ processing %dth frame\n", *iframe+1);

    process_elementary_frame(pwr, modsq, long_int, snrm,
                             bisp1, PSF1, FT_object, DetectorTransf,
                             Np, seed, nph, ir, nbeta, ngamma,
                             output_demo, out_prefix);
    (*iframe)++;
    if(*iframe >= nframes){
       free(PupilFunc1); free(PupilFunc2);
       free(PSF1); free(PSF2);
       return(0);
       }

/* Processing a second frame with imag part of phase screen:
*/
    process_elementary_frame(pwr, modsq, long_int, snrm,
                             bisp1, PSF2, FT_object, DetectorTransf,
                             Np, seed, nph, ir, nbeta, ngamma,
                             output_demo, out_prefix);
    (*iframe)++;
    if(*iframe >= nframes){
       free(PupilFunc1); free(PupilFunc2);
       free(PSF1); free(PSF2);
       return(0);
       }
 
/************************* end of major loop **********************************/
  } /* EOF loop on yy */
} /* EOF loop on xx */

free(PupilFunc1); 
free(PupilFunc2);
free(PSF1);
free(PSF2);
return(0);
}
/***************************************************************************
* read_param_file: read input parameter file
*
* INPUT:
* parameter_fname: name of input parameter file 
*
* OUTPUT:
* Np: size of elementary frames (=Ns/2)
* seed: seed for random number generator (large odd integer)
* L0, lo: external and internal scale lengths of turbulence
* lambda_cent, bandwidth, nsamp: wavelength parameters: central wavelength (m),
*                          bandwidth (m) and sampling points
* r_zero: Fried's parameter (m) for the central wavelength
* ncomp: number of components for the source
* s_int, s_pos: source information (intensity and position in arcsec)
* float s_int[MAXCOMP], s_pos[MAXCOMP][2]
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* float c_pos[MAXAP][2], c_rad[MAXAP];
* fov: field of view (in arcseconds)
* nframes: number of frames to be simulated
* no_dist: flag set to one if no atmospheric disturbance
* MAXCOMP: maximum number of components in source 
* MAXAP: maximum number of apertures
*
***************************************************************************/
static int read_param_file(char *parameter_fname, int *Np, long *seed, int *nph,
                           float *L0, float *lo, float *lambda_cent, 
                           float *bandwidth, int *nsamp, int *ncomp, 
                           float *s_int, float s_pos[][2], int *ncircle, 
                           float c_pos[MAXAP][2], float *c_rad, int *c_typ, 
                           float *fov, float *r_zero, int *nframes, 
                           int *no_dist, char *out_prefix) 
{
char buffer[80];
int italk = 1;
register int i;
FILE *fp_in;
 
/* Initialize no_dist to default value:
*/
*no_dist = 0;
*nph = 200;
*Np = 128;

/* Open the input file:
*/
if((fp_in = fopen(parameter_fname,"r")) == NULL){
  fprintf(stderr,"read_param_file/Fatal error opening input file: %s\n", 
          parameter_fname);
  exit(-1);
  }
 
while(!feof(fp_in)) {
  fgets(buffer,80,fp_in);
  if(buffer[0] != '!'){
/* Keyword "Size of images", parameter: Np 
*/
    if(!strncmp(buffer, "Size", 4)){
      if(fscanf(fp_in, "%d", Np) != 1) {
        fprintf(stderr, "read_param/Fatal error reading Size!\n");
        exit(-1);
        }
      if(italk) printf("Size (of images):  Np = %d pixels\n", *Np);
    }
/* Keyword "Photons per frame", parameter: nph
*/
    if(!strncmp(buffer, "Photons", 7)){
      if(fscanf(fp_in, "%d", nph) != 1) {
        fprintf(stderr, "read_param/Fatal error reading nph!\n");
        exit(-1);
        }
      if(italk) printf("Photons per frame :  nph = %d photons/frame\n", *nph);
    }
/* Keyword "Field of view", parameter: fov 
*/
    if(!strncmp(buffer, "Field", 5)){
      if(fscanf(fp_in, "%f", fov) != 1) {
        fprintf(stderr, "read_param/Fatal error reading fov!\n");
        exit(-1);
        }
      if(italk) printf("Field of view:  fov = %.4f arcseconds\n", *fov);
    }
/* Keyword "Random" (random generator), parameter: seed
*/
    if(!strncmp(buffer, "Random", 6)){
      if(fscanf(fp_in, "%ld", seed) != 1) {
        fprintf(stderr, "read_param/Fatal error reading seed!\n");
        exit(-1);
        }
      if(italk) printf("Random generator:  seed = %ld\n", *seed);
    }
/* Keyword "Turbulence scale", parameters: outer and inner scales
*/
    if(!strncmp(buffer, "Turbulence", 10)){
      if(fscanf(fp_in, "%f %f", L0, lo) != 2) {
        fprintf(stderr, "read_param/Fatal error reading L0 and lo!\n");
        exit(-1);
        }
      if(italk) printf("Turbulence scales:  L0 = %.4f  lo = %.4f\n", *L0, *lo);
    }
/* Keyword "Fried radius", parameter: r0 
*/
    if(!strncmp(buffer, "Fried", 5)){
      if(fscanf(fp_in, "%f", r_zero) != 1) {
        fprintf(stderr, "read_param/Fatal error reading Fried (i.e. r0)!\n");
        exit(-1);
        }
      if(italk) printf("Fried radius:  r0 = %.4f\n", *r_zero);
    }
/* Keyword "Wavelength", parameters: 
* lambda_cent, bandwidth, nsamp: central wavelength (nm), bandwidth (nm) 
* and sampling points
*/
    if(!strncmp(buffer, "Wavelength", 10)){
      if(fscanf(fp_in, "%f %f %d", lambda_cent, bandwidth, nsamp) != 3) {
        fprintf(stderr, "read_param/Fatal error wavelength parameters!\n");
        exit(-1);
        }
/* Conversion from nm to m: */
      *lambda_cent *= 1.e-9;
      *bandwidth *= 1.e-9;
      if(italk) printf("lambda_cent=%e (m) bandwidth=%e (m) n_samples=%d\n", 
             *lambda_cent, *bandwidth, *nsamp);
    }
/* Keyword "Source", parameters: number of points, followed by
*   Intensity1, X1, Y1 
*   Intensity2, X2, Y2... 
*/
    if(!strncmp(buffer, "Source", 6)){
      if(fscanf(fp_in, "%d", ncomp) != 1) {
        fprintf(stderr, "read_param/Fatal error reading ncomp!\n");
        exit(-1);
        }
      if(italk) printf("Source with %d components\n", *ncomp);
       for(i = 0; i < *ncomp; i++) {
          if(fscanf(fp_in, "%f %f %f", 
                    &(s_int[i]), &(s_pos[i][0]), &(s_pos[i][1])) != 3) {
             fprintf(stderr, "read_param/Fatal error reading ncomp!\n");
             exit(-1);
             }
          if(italk) printf("Intensity, x, y: %.4f, %4f %.4f\n", s_int[i],
                           s_pos[i][0], s_pos[i][1]);
          }
    }
/* Keyword "Aperture", parameters: number of apertures, followed by
* Xcenter, Ycenter, diameter, type (0 or 1)
*/
    if(!strncmp(buffer, "Aperture", 8)){
      if(fscanf(fp_in, "%d", ncircle) != 1) {
        fprintf(stderr, "read_param/Fatal error reading ncircle!\n");
        exit(-1);
        }
      if(italk) printf("Telescope aperture with %d circles\n", *ncircle);
       for(i = 0; i < *ncircle; i++) {
          if(fscanf(fp_in, "%f %f %f %d", &(c_pos[i][0]), &(c_pos[i][1]), 
                    &(c_rad[i]), &(c_typ[i])) != 4) {
             fprintf(stderr, "read_param/Fatal error reading aperture parameters!\n");
             exit(-1);
             }
/* Conversion from diameter to radius: */
          c_rad[i] /= 2.;
          if(italk) printf("x, y, rad, type: %.4f, %4f %.4f %d\n", c_pos[i][0],
                           c_pos[i][1], c_rad[i], c_typ[i]);
          }
    }
/* Keyword "Frames", parameters: nber of frames 
*/
    if(!strncmp(buffer, "Frames", 6)){
      if(fscanf(fp_in, "%d", nframes) != 1) {
        fprintf(stderr, "read_param/Fatal error reading nframes!\n");
        exit(-1);
        }
      if(italk) printf("Nber of frames  nframes = %d\n", *nframes);
    }
/* Keyword "No disturbance" (no atmosphere)
*/
    if(!strncmp(buffer, "No disturb", 10)){
      *no_dist = 1;
    }
  }  /* EOF buffer[0] != '!' */
}  /* EOF while */

 if(italk && *no_dist) printf("OK: no atmospheric disturbance\n");

fclose(fp_in);
return(0);
}
/*****************************************************************
* Recenter Fourier transform
* and shift origin from (0,0) to (Ns/2+1, Ns/2+1) or inversely
*
* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*****************************************************************/
int JLP_RECENT_CMPLX(complex *psc, int Ns, int idim)
{
complex tmp;
register int ix, iy;
/* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*/
for(ix = 0; ix < Ns/2; ix++) {
  for(iy = 0; iy < Ns/2; iy++) {
/* Flip quadrants 2 and 3:
*/
     tmp = psc[ix + Ns/2  + iy * idim];
     psc[ix + Ns/2 + iy * idim] = psc[ix + (iy + Ns/2) * idim];
     psc[ix + (iy + Ns/2) * idim] = tmp;
/* Flip quadrants 1 and 4:
*/
     tmp = psc[ix + Ns/2  + (iy + Ns/2) * idim];
     psc[ix + Ns/2 + (iy + Ns/2) * idim] = psc[ix + iy * idim];
     psc[ix + iy * idim] = tmp;
    }
  }

 return(0);
}
/*****************************************************************
* Recenter Fourier transform
* and shift origin from (0,0) to (Ns/2+1, Ns/2+1) or inversely
*
* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*****************************************************************/
int JLP_RECENT_FLOAT(float *psc, int Ns, int idim)
{
float tmp;
register int ix, iy;
/* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*/
for(ix = 0; ix < Ns/2; ix++) {
  for(iy = 0; iy < Ns/2; iy++) {
/* Flip quadrants 2 and 3:
*/
     tmp = psc[ix + Ns/2  + iy * idim];
     psc[ix + Ns/2 + iy * idim] = psc[ix + (iy + Ns/2) * idim];
     psc[ix + (iy + Ns/2) * idim] = tmp;
/* Flip quadrants 1 and 4:
*/
     tmp = psc[ix + Ns/2  + (iy + Ns/2) * idim];
     psc[ix + Ns/2 + (iy + Ns/2) * idim] = psc[ix + iy * idim];
     psc[ix + iy * idim] = tmp;
    }
  }

 return(0);
}
/*****************************************************************
* Recenter Fourier transform
* and shift origin from (0,0) to (Ns/2+1, Ns/2+1) or inversely
*
* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*****************************************************************/
int JLP_RECENT_DOUBLE(double *psc, int Ns, int idim)
{
double tmp;
register int ix, iy;
/* Quadrants 
*     3 4        2 1 
*     1 2    to  4 3 
*/
for(ix = 0; ix < Ns/2; ix++) {
  for(iy = 0; iy < Ns/2; iy++) {
/* Flip quadrants 2 and 3:
*/
     tmp = psc[ix + Ns/2  + iy * idim];
     psc[ix + Ns/2 + iy * idim] = psc[ix + (iy + Ns/2) * idim];
     psc[ix + (iy + Ns/2) * idim] = tmp;
/* Flip quadrants 1 and 4:
*/
     tmp = psc[ix + Ns/2  + (iy + Ns/2) * idim];
     psc[ix + Ns/2 + (iy + Ns/2) * idim] = psc[ix + iy * idim];
     psc[ix + iy * idim] = tmp;
    }
  }

 return(0);
}
/***********************************************************************
* To create an image corresponding to the source
*
* INPUT:
* Np: size of elementary frames (=Ns/2)
* fov: Field of view
* ncomp: number of components for the source
* s_int, s_pos: source information (intensity and position in arcsec)
* out_prefix: prefix used to generate names of output files
*
* OUTPUT:
* object(Np,Np): array with ones at the location of the components
* FT_object(Np,Np): Fourier Transform of the object
***********************************************************************/
static int create_object(float **object, complex **FT_object, int Np, 
                         float fov, int ncomp, float *s_int, float s_pos[][2],
                         char *out_prefix)
{
float Scale_i;
int icent, ix, iy, nn[2];
char out_file[60], out_comments[80];
register int i, k;

*object = (float *)malloc(Np * Np * sizeof(float));
*FT_object = (complex *)malloc(Np * Np * sizeof(complex));

/* Scale_i: scale of image in pixels/arcsecond
*/
  Scale_i = (float)Np / fov; 
  icent = Np / 2;
  printf(" create_object/Scale of image: %.3f pixels/arcsec\n", 
          Scale_i); 
 
/* Initializing the arrays:
*/
   for(i = 0; i < Np * Np; i++) (*object)[i] = 0.;
   for(i = 0; i < Np * Np; i++) (*FT_object)[i].im = 0.;

/* Loop on all components:
* (there are ncomp components in the description of the source)
*
*/
for(k = 0; k < ncomp; k++) {
          ix = icent + (int)(Scale_i * s_pos[k][0] + 0.5);
          iy = icent + (int)(Scale_i * s_pos[k][1] + 0.5);
          ix = MINI(Np, MAXI(1, ix));
          iy = MINI(Np, MAXI(1, iy));
          printf("create_object/Component #%d ix=%d iy=%d s_int=%f\n", k, ix, iy, s_int[k]);
          (*object)[ix + iy * Np] += s_int[k]; 
  }

/*
 Transfer object to FT_object for Fourier transform
*/
   for(i = 0; i < Np * Np; i++) (*FT_object)[i].re = (*object)[i];

    nn[0] = Np;
    nn[1] = Np; 
    fourn_for_C_arrays((float*)(*FT_object), nn, 2, 1);
 
/* Output object for a diagnostic:
*/
    sprintf(out_file, "%s_object", out_prefix);
    strcpy(out_comments, "Source");
    JLP_WRITEIMAG(*object, &Np, &Np, &Np, out_file, out_comments);

return(0);
}
/**************************************************************************
* Compute two elementary PSFs PSF1 and PSF2
* from a subwindow of the phase screen starting at x,y
* Take into account wavelength smearing due to finite bandpass
*
* INPUT:
* psc: phase screen to be used for computing the elementary PSFs 
* lambda_cent: central wavelength (in meters)
* bandwidth: half width of filter used (in meters) 
* nsamp: number of sampling points to take the finite bandwidth into account
* xx, yy : position of the center of the subwindow used to extract 
*          a phase screen from psc 
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* float c_pos[MAXAP][2], c_rad[MAXAP];
*  
* OUTPUT:
* PSF1, PSF2: elementary PSFs corresponding to real and imaginary
*                  parts of the phase screen
**************************************************************************/
static int compute_elementary_PSFs(complex *psc, float xx, float yy, float Ls,
                                   float lambda_cent, float bandwidth, 
                                   int nsamp, int ncircle, float c_pos[][2], 
                                   float *c_rad, int *c_typ, 
                                   complex *PupilFunc1, complex *PupilFunc2,
                                   float *PSF1, float *PSF2, int Ns, int Np, 
                                   int *output_demo, char *out_prefix)
{
float *work;
float factor, xa, ya, Scale_s, ww;
float lambda_min, l_step, xcent, ycent;
int ix, iy, lx, ly, nn[2], ix_start, iy_start;
char out_file[60], out_comments[80];
complex  u_zero;
register int i, isamp, icirc;

Scale_s = Ls/(float)Ns;

/* Reset image arrays
*/
for(i = 0; i < Np*Np; i++) {
    PSF1[i] = 0.;
    PSF2[i] = 0.;
  }
 
/* Maximum variation of lambda is bandwidth 
*/
lambda_min = lambda_cent - bandwidth / 2.;
l_step = bandwidth / (float)nsamp;

ix_start = (int)(xx / Scale_s);
iy_start = (int)(yy / Scale_s);

#ifdef DEBUG
printf("compute_elementary_PSF/ix_start=%d iy_start=%d (xx=%f yy=%f)\n", 
       ix_start, iy_start, xx, yy);
#endif

/******************* BEGIN OF LOOP1 ON ISAMP *************************
******* on all the (nsamp) wavelengths belonging to the bandwidth ****
* nsamp data points within the bandpass
*/
for(isamp = 0; isamp < nsamp; isamp++) {
/*
* JLP's formula since psc is proportional to kk_cent = 2 pi / lambda_cent:
*/
    factor =  lambda_cent / (lambda_min + isamp * l_step);
#ifdef DEBUG_
printf("isamp=%d factor=%f lstep=%f\n", isamp, factor, l_step);
#endif

/* Reset temporary arrays
*/
  u_zero = float_to_cplx(0., 0.);
  for(i = 0; i < Np*Np; i++){
    PupilFunc1[i] = u_zero; 
    PupilFunc2[i] = u_zero; 
  }
 
/* Loop on all the circular apertures
*/
  for(icirc = 0; icirc < ncircle; icirc++){
      xcent = xx + c_pos[icirc][0];
      ycent = yy + c_pos[icirc][1];
#ifdef DEBUG_
      printf("Circle #%d centered at xx+cposx=%.3f, yy+cposy=%.3f\n", 
             icirc, xcent, ycent);
#endif
 
/*
*  Filling Pupil data with exp (i phase) with aperture phases
*/
     for(iy = 0; iy < Np; iy++) {
        for(ix = 0; ix < Np; ix++) {
          xa = Scale_s * (ix - Np/2) - c_pos[icirc][0]; 
          ya = Scale_s * (iy - Np/2) - c_pos[icirc][1]; 
/* Check if current point lies inside the circle: */
          if( SQUARE(xa)+SQUARE(ya) <= SQUARE(c_rad[icirc]) ) { 
/* NB: psc is proportional to kk = 2 pi / lambda
* correction by factor = lambda_cent / (lambda_cent + dlambda):
*/
            lx = (ix - Np/2) + ix_start;
            ly = (iy - Np/2) + iy_start;
            PupilFunc1[ix + iy * Np].re = c_typ[icirc] 
                                  * cos((double)factor * psc[lx + ly * Ns].re);
            PupilFunc1[ix + iy * Np].im = c_typ[icirc] 
                                  * sin((double)factor * psc[lx + ly * Ns].re);
            PupilFunc2[ix + iy * Np].re = c_typ[icirc] 
                                  * cos((double)factor * psc[lx + ly * Ns].im);
            PupilFunc2[ix + iy * Np].im = c_typ[icirc] 
                                  * sin((double)factor * psc[lx + ly * Ns].im);
/* Possibility of outputing the mask used to extract the phase screen: 
*/
#ifdef OUTPUT_PSC_MASK
          psc[lx + ly * Ns].re =0.;
          psc[lx + ly * Ns].im =0.;
#endif
          } /*EOF if inside the circle */
       }  /* EOF loop on ix */
     }  /* EOF loop on iy */
   } /* EOF loop on icirc (to synthesize the aperture ... ) */

/* Output real part of complex pupil function: 
*/
     if(*output_demo == 1){
       work = (float *)malloc(Np * Np * sizeof(float));
       for(i = 0; i < Np * Np; i++) work[i] = PupilFunc1[i].re;
       sprintf(out_file, "%s_PupilFunc", out_prefix);
       strcpy(out_comments, "Real part of exp (i phase screen)");
       JLP_WRITEIMAG(work, &Np, &Np, &Np, out_file, out_comments);
       free(work);
       *output_demo = 2;
       }

#ifdef OUTPUT_PSC_MASK
       work = (float *)malloc(Ns* Ns * sizeof(float));
       for(i = 0; i < Ns * Ns; i++) work[i] = psc[i].re;
       sprintf(out_file, "%s_psc_test", out_prefix);
       strcpy(out_comments, "Real part of phase screen");
       JLP_WRITEIMAG(work, &Ns, &Ns, &Ns, out_file, out_comments);
       free(work);
#endif

/* Then simulate propagation from aperture plane to detector plane
*/
 
/* SHOULD USE Numerical recipee's fourn routine (in this directory) 
*/
/* fourn(data,nn,ndim,isign)
*/
     nn[0] = Np;
     nn[1] = Np;
     fourn_for_C_arrays((float*)PupilFunc1, nn, 2, 1);
     fourn_for_C_arrays((float*)PupilFunc2, nn, 2, 1);
 
/* The intensity of the Fourier Transform of the Pupil function
* corresponds to an elementary image
*
* Recenter images (necessary for the subsequent convolution ...) 
*/
     JLP_RECENT_CMPLX(PupilFunc1, Np, Np);
     JLP_RECENT_CMPLX(PupilFunc2, Np, Np);
 
/* Form speckle/fringe patterns
* i.e., intensity of the complex field (PupilFunc1 or PupilFunc2)
*/
     ww = 1. / (float)(Np * Np);
     for(i = 0; i < Np * Np; i++) {
        PSF1[i] += ww * cplx_sqmodulus(PupilFunc1[i]);
        PSF2[i] += ww * cplx_sqmodulus(PupilFunc2[i]);
     }
 
/******************* END OF LOOP1 ON ISAMP ********************/
} /* EOF loop on isamp (to emulate a finite bandpass with many wavelengths) */

/* Output PSF1:
*/
  if(*output_demo == 2){
     sprintf(out_file, "%s_psf", out_prefix);
     strcpy(out_comments, "Intensity of complex field in image plane");
     JLP_WRITEIMAG(PSF1, &Np, &Np, &Np, out_file, out_comments);
  }

return(0);
}
/********************************************************************
* Compute an elementary frame by convolving PSF1 with source
* and pile up the corresponding power spectrum and long integration 
*
* INPUT:
* PSF1: elementary PSF 
* FT_object: Fourier Transform of the source
* DetectorTransf: detector transfer function (i.e., modulus of Fourier 
*                 Tranform of PSF)
* nph: number of photons per frame
* Np: size of elementary frames (=Ns/2)
*
* INPUT/OUTPUT:
* pwr: mean power spectrum of analog elementary frames 
* modsq: mean power spectrum with photon clipping and detector photon response
* long_int: long integration with photon clipping and detector photon response 
* output_demo: value of a dummy variable used to know whether
*                an image has to be output
********************************************************************/
static int process_elementary_frame(double *pwr, double *modsq, 
                               double *long_int, double *snrm, double *bisp1, 
                               float *PSF1, complex *FT_object, 
                               complex *DetectorTransf, int Np, long *seed, 
                               int nph, int ir, int nbeta, int ngamma, 
                               int *output_demo, char *out_prefix)
{
/* Temporary arrays
*/
float *ary1, *ary1_ph;
double *re1, *im1;
complex *cary1, *cary1_ph, *TransferFunc;
/* Miscellaneous: */
int nn[2], i, nphr;
char out_file[60], out_comments[80];

ary1 = (float *)malloc(Np * Np * sizeof(float));
ary1_ph = (float *)malloc(Np * Np * sizeof(float));
re1 = (double *)malloc(Np * Np * sizeof(double));
im1 = (double *)malloc(Np * Np * sizeof(double));
cary1 = (complex *)malloc(Np * Np * sizeof(complex));
cary1_ph = (complex *)malloc(Np * Np * sizeof(complex));
TransferFunc = (complex *)malloc(Np * Np * sizeof(complex));

/* Compute the transfer function corresponding to the input PSF: */
     for(i = 0; i < Np * Np; i++){ 
          TransferFunc[i].re = PSF1[i]; 
          TransferFunc[i].im = 0.; 
          }
     nn[0] = Np;
     nn[1] = Np;
     fourn_for_C_arrays((float *)TransferFunc, nn, 2, 1);

/* Convolution is done in Fourier domain:
* Product is a complex product here!
*/
     for(i = 0; i < Np * Np; i++) 
             cary1[i] = cplx_product(TransferFunc[i], FT_object[i]); 

/* Pile up the analog power spectrum:
* (NB: pwr is set to zero in init_arrays...)
*/
     for(i = 0; i < Np * Np; i++) 
            pwr[i] += SQUARE(cary1[i].re) + SQUARE(cary1[i].im); 
 
/* Inverse Fourier Transform to create an analog image:
*/
     nn[0] = Np;
     nn[1] = Np;
     fourn_for_C_arrays((float *)cary1, nn, 2, -1);

     for(i = 0; i < Np * Np; i++) ary1[i] = cary1[i].re; 
     JLP_RECENT_FLOAT(ary1, Np, Np);

/* Output one of those elementary frames:
*/
     if(*output_demo == 2) {
       sprintf(out_file, "%s_elem_frame_an", out_prefix);
       strcpy(out_comments, "PSF1 convolved by source");
       JLP_WRITEIMAG(ary1, &Np, &Np, &Np, out_file, out_comments);
       }

/* Photon clipping:
*/
     photon_clipping(ary1, ary1_ph, Np, nph, &nphr, seed);
 
/* Output one of those photon clipped images:
*/
     if(*output_demo == 2) {
       sprintf(out_file, "%s_elem_frame_ph", out_prefix);
       strcpy(out_comments, "PSF1 convolved by source and photon clipped");
       JLP_WRITEIMAG(ary1_ph, &Np, &Np, &Np, out_file, out_comments);
       }

/* Fourier transform to compute mosq:
*/
     for(i = 0; i < Np * Np; i++){ 
       cary1_ph[i].re = ary1_ph[i];
       cary1_ph[i].im = 0.;
       } 

     nn[0] = Np;
     nn[1] = Np;
     fourn_for_C_arrays((float *)cary1_ph, nn, 2, 1);

/* Convolution with detector response (i.e product in Fourier domain):
* Warning: product is a complex product here!
*/
     for(i = 0; i < Np * Np; i++) 
             cary1[i] = cplx_product(cary1_ph[i], DetectorTransf[i]); 

/* Loading Fourier transform to re,im arrays:
*/
     for(i = 0; i < Np * Np; i++) {
       re1[i] = cary1[i].re;
       im1[i] = cary1[i].im;
       }

/*
* Processing the spectrum of this image now:
* int bispec3(double *re, double *im, double *modsq, double *snrm,
*            INT4 *nx, INT4 *ny, double *bispp, INT4 *ir,
*            INT4 *nbeta, INT4 *ngamma);
*/
     BISPEC3(re1, im1, modsq, snrm, &Np, &Np, bisp1, &ir, &nbeta, &ngamma);
 
/* Back to direct domain:
*/
     nn[0] = Np;
     nn[1] = Np;
     fourn_for_C_arrays((float *)cary1, nn, 2, -1);

/* Pile up the long integration 
* (NB: long_int is set to zero in init_arrays...)
*/
     for(i = 0; i < Np * Np; i++) {
          ary1[i] = cary1[i].re; 
          long_int[i] += ary1[i];
        }

/* Output one of those elementary frames (photon clipped 
* and convolved with detector response) 
*/
     if(*output_demo == 2) {
       sprintf(out_file, "%s_elem_frame", out_prefix);
       strcpy(out_comments, "Simulation of output from detector");
       JLP_WRITEIMAG(ary1, &Np, &Np, &Np, out_file, out_comments);
       *output_demo = 3;
       }

/* Free memory */
free(ary1);
free(ary1_ph);
free(re1);
free(im1);
free(cary1);
free(cary1_ph);
free(TransferFunc);

return(0);
}
/**************************************************************************
* Before calling poidev, the sum of in_array has to be normalised to nph
*
* INPUT:
* in_array: analog elementary frame according to Poisson statistics
* Np: size of elementary frames (=Ns/2)
* nph: number of photons/frame wanted
* seed: seed used for Poisson random generator
*
* OUTPUT:
* out_array: elementary frame that was digitized according to Poisson statistics
* nphr: number of photons of the output frame
**************************************************************************/
static int photon_clipping(float *in_array, float *out_array, int Np, int nph, 
                           int *nphr, long *seed)
{
double sum;
float w_scale;
register int i;

/* WARNING: before calling poidev, the sum of prf has to be normalised to nph
*/
/* First normalise the sum of image (ary1) to nph photons (but still analogic) 
*/
   sum = 0.;
   for(i = 0; i < Np * Np; i++) sum += in_array[i];
    
   w_scale = (float)nph / sum;
   for(i = 0; i < Np * Np; i++) in_array[i] *= w_scale;


/* Now using Poisson generator
*/
  *nphr = 0;
   for(i = 0; i < Np * Np; i++) {
     JLP_POISSON(&in_array[i], &out_array[i], seed);
/* Computing total number of photons: 
*/
     *nphr += out_array[i];
     }

#ifdef DEBUG_
   printf(" # of photons after clipping: %d\n", *nphr); 
#endif

return(0);
}
/**************************************************************************
* INPUT:
* Ls: width of phase screen (in m)
* Ns: width of phase screen (in pixels)
* Np: width of image (in pixels)
* ncirc: number of circles for the aperture
* c_pos, c_rad, c_typ: aperture information (coordinate, radius, type)
* float c_pos[MAXAP][2], c_rad[MAXAP];
* MAXAP: maximum number of apertures
*
* OUTPUT:
* frm_per_screen_width: number of frames per phase screen in X or in Y
* x_min, x_max, xstep : sampling parameters to extract a phase screen 
*                       from a larger one 
* y_min, y_max, ystep : sampling parameters to extract a phase screen 
*                       from a larger one 
**************************************************************************/
static int compute_sampling(float Ls, int Ns, int Np, 
                            int *frm_per_screen_width, 
                            float *x_min, float *x_max, float *xstep, 
                            float *y_min, float *y_max, float *ystep, 
                            int ncircle, float c_pos[][2], float *c_rad)
{
float aperture_diameter;
register int i;

aperture_diameter = 0.;
 for(i = 0; i < ncircle; i++) {
    aperture_diameter = MAXI(aperture_diameter, 2. * c_rad[i]);
  }

/* Np is Ns/2 hence screen size is Ls/2; 
*/
 if(aperture_diameter > Ls/2.) {
  fprintf(stderr, "Fatal error, maximum aperture diameter= %.3f m\n", 
           aperture_diameter);
  fprintf(stderr, " which is larger than the half width of the phase screens: %.3f m ! \n", Ls/2.);
  exit(-1);
  }

/* Example of an annular aperture (in the parameter file):
Aperture
2
0.0 0.0 2.0  1
0.0 0.0 0.4  0
*/

/* x_min,ymin: lower left coordinates of the subwindow to be used
* to extract a small phase screen from the large phase screen
* First compute x_min, and y_min, by looking at the minimum value found
* when using all the apertures:
*/
   *x_min = Ls; 
   *y_min = Ls;
   for(i = 0; i < ncircle; i++) {
        *x_min = MINI(*x_min, c_pos[i][1] - c_rad[i]);
        *y_min = MINI(*y_min, c_pos[i][2] - c_rad[i]);
    }

/* 
* If x_min is negative, it means that the center should be taken
* large enough to avoid extracting data point outside of the frame 
* (same for y_imin)
* Tolerance is 0.1 m to avoid zero at the edge...
*/
   if(*x_min < 0.) 
       *x_min = 0.1 - *x_min;
   else
       *x_min = 0.1;

   if(*y_min < 0.) 
       *y_min = 0.1 - *y_min;
   else
       *y_min = 0.1;

/* Then x_max, y_max:
*/
   *x_max = 0.;
   *y_max = 0.;
   for(i = 0; i < ncircle; i++) {
        *x_max = MAXI(*x_max, c_pos[i][1] + c_rad[i]);
        *y_max = MAXI(*y_max, c_pos[i][2] + c_rad[i]);
    }

   if(*x_max > 0.) 
       *x_max = Ls - 0.1 - *x_max;
   else
       *x_max = Ls - 0.1;

   if(*y_max > 0.) 
       *y_max = Ls - 0.1 - *y_max;
   else
       *y_max = Ls - 0.1;

/* Compute the number of frames to be extracted from the phase screen 
* I assume here that the user wants to extract non-overlapping pupils
*/
 *frm_per_screen_width = (int)(Ls/aperture_diameter);

/* Then xstep, ystep:
*/
  if(*frm_per_screen_width > 1){
     *xstep = (*x_max - *x_min) / (*frm_per_screen_width - 1);
     *ystep = (*y_max - *y_min) / (*frm_per_screen_width - 1);
  } else {
     *x_max = *x_min;
     *y_max = *y_min;
     *xstep = 1000.;
     *ystep = 1000.;
  }

/* Error analysis:
*/
 if(*xstep <= 0. || *ystep <= 0. || *x_min >= Ls || *y_min >= Ls
      || *x_max < 0. || *y_max < 0. || *frm_per_screen_width == 0) {
   fprintf(stderr, "Fatal error computing sampling parameter\n");
   fprintf(stderr, "x_min, x_max, y_min, y_max: %.2f %.2f %.2f %.2f (m)\n", 
            *x_min, *x_max, *y_min, *y_max);
   fprintf(stderr, "xstep, ystep: %.2f %.2f (m)\n", *xstep, *ystep);
   fprintf(stderr, "Np=%d Ns=%d frm_per_screen_width=%d (pixels)\n", 
           Np, Ns, *frm_per_screen_width);
   exit(-1);
   }

   printf("x_min, x_max, y_min, y_max: %.2f %.2f %.2f %.2f (m) (Ls=%.2f m)\n", 
            *x_min, *x_max, *y_min, *y_max, Ls);
   printf("xstep, ystep: %.2f %.2f (m)   (aperture diameter=%.2f m)\n", 
           *xstep, *ystep, aperture_diameter);
return(0);
}
/***********************************************************************
* To read the FITS file containing the photon response of the detector 
*
* INPUT:
*  photon_modsq: name of the FITS file containing the detector response
*                (= '0' when no modsq is wanted)
*
* OUTPUT:
*  DetectorTransf: detector transfer function 
***********************************************************************/
static int read_photon_response(char *photon_modsq_fname, int Np, 
                                complex **DetectorTransf)
{
float *detector_modsq, w_scale;
int istat;
INT4 nx1, ny1;
INT_PNTR pntr_ima;
char comments[80];
register int i;

*DetectorTransf = (complex *)malloc(Np * Np * sizeof(complex));
for(i = 0; i < Np*Np; i++)  (*DetectorTransf)[i].im = 0.;

/* No modsq file is available or wanted: 
* I then assume the detector has a uniform response in Fourier domain 
* (i.e. a Delta function in the image plane) */
if((photon_modsq_fname[0] == '0') || !(*photon_modsq_fname)) {

  printf("No modsq file entered, will assume that the dectector is perfect\n");
  for(i = 0; i < Np*Np; i++)  (*DetectorTransf)[i].re = 1.;

} else {

/* Read input file: */
   istat=JLP_VM_READIMAG1(&pntr_ima, &nx1, &ny1, photon_modsq_fname, comments);
   detector_modsq = (float *)pntr_ima;
    if(istat != 0) {
      fprintf(stderr, "read_photon_response/Fatal error reading %s (istat=%d)\n", 
              photon_modsq_fname, istat);
      exit(-1);
     }
    if(nx1 != Np || ny1 != Np) {
      fprintf(stderr, "read_photon_response/Fatal error, inconsistent size\n");
      fprintf(stderr, " nx1=%d ny1=%d whereas Np=%d (from parameter file)\n",
              nx1, ny1, Np);
      exit(-1);
     }

/* Shift zero frequency to [0,0]: */ 
   JLP_RECENT_FLOAT(detector_modsq, Np, Np);

/* Normalization to unity: */
   w_scale = 1./sqrt(detector_modsq[0]);
   if(w_scale == 0.) {
     fprintf(stderr, 
             "read_photon_response/Fatal error modsq[central_freq]=0. !\n"); 
     exit(-1);
     }

/* Transfer function corresponds to sqrt(modsq): 
*/
  for(i = 0; i < Np*Np; i++) 
      (*DetectorTransf)[i].re = w_scale * sqrt(detector_modsq[i]);
  }

return(0);
}
/********************************************************************
********************************************************************/
static int output_results(double *pwr, double *modsq, double *long_int, 
                          double *snrm, double *bisp1, int Np, int ngamma, 
                          int nframes, char *parameter_fname, char *out_prefix)
{
double w_scale, w1, w2, w3, *bisp2;
char out_file[60], out_comments[160];
int i1, i2, i3, i4, ndim;
register int i;

bisp2 = (double *)malloc(ngamma *3 * sizeof(double));

/* Recenter Fourier transforms:
*/
    JLP_RECENT_DOUBLE(pwr, Np, Np);
    JLP_RECENT_DOUBLE(modsq, Np, Np);
    JLP_RECENT_DOUBLE(snrm, Np, Np);

/* Compute the mean of the power spectrum, square modulus, SNR,
* long integration and bispectrum: 
*/
   w_scale = 1./(float)nframes;
   for(i = 0; i < Np * Np; i++) {
      pwr[i] *= w_scale;
      modsq[i] *= w_scale;
      snrm[i] *= w_scale;
      long_int[i] *= w_scale;
     }
   for(i = 0; i < 4 * ngamma; i++) bisp1[i] *= w_scale;

/* Computing the variance of the bispsctrum:
*/
   for(i = 0; i < ngamma; i++) { 
/* sum(re)
*/
            i1 = 4*i ;
/* sum(im)
*/
            i2 = 4*i + 1;
/* sumsq(re)
*/
            i3 = 4*i + 2;
/* sumsq(im)
*/
            i4 = 4*i + 3; 
/* variance of re:
*/
            w1 = bisp1[i3] - SQUARE(bisp1[i1]);
/* variance of im:
*/
            w2 = bisp1[i4] - SQUARE(bisp1[i2]);
/* Full variance (variance of real + variance of imag):
*/
            w1 = w1 + w2;
            w1 = MAXI(w1, 1.e-10);
/* Standard deviation:
*/
            w1 = sqrt(w1);
/* Computing the phase factor on the first two lines:
*/
            w3 = SQUARE(bisp1[i1]) + SQUARE(bisp1[i2]);
            w3 = MAXI(w3, 1.e-10);
            w3 = sqrt(w3);
/*  In old versions, I stored the phase factor, but now I keep the modulus
 * (to be able to correct the photon noise...) 
 *          bisp2(lx,1) = bisp1(i1)/w3
 *          bisp2(lx,2) = bisp1(i2)/w3
*/
/* As in vcrb: I save the bispectrum instead of the phase factor:
* to be able to display the bispectrun with "bisp_to_image.c"
* and to be compatible with other programs:
*/
            bisp2[i] = bisp1[i1];
            bisp2[i + ngamma] = bisp1[i2];
/* SNR of bispectrum on the third line:
*/
            bisp2[i + 2*ngamma] = w3/w1;
        }
 
/* Output the mean long integration to a FITS file:
*/
       sprintf(out_file, "%s_l", out_prefix);
       sprintf(out_comments, "Mean long integ. with %s", parameter_fname);
       JLP_D_WRITEIMAG(long_int, &Np, &Np, &Np, out_file, out_comments);

/* Output the analog power spectrum to a FITS file:
*/
       sprintf(out_file, "%s_power", out_prefix);
       sprintf(out_comments, "Analog power spectrum with %s", parameter_fname);
       JLP_D_WRITEIMAG(pwr, &Np, &Np, &Np, out_file, out_comments);

/* Output the SNR of the power spectrum to a FITS file:
*/
       sprintf(out_file, "%s_snrm", out_prefix);
       sprintf(out_comments, "SNRM of modsq with %s", parameter_fname);
       JLP_D_WRITEIMAG(snrm, &Np, &Np, &Np, out_file, out_comments);

/* Output the fully simulated power spectrum to a FITS file:
*/
       sprintf(out_file, "%s_m", out_prefix);
       sprintf(out_comments, "modsq with %s", parameter_fname);
       JLP_D_WRITEIMAG(modsq, &Np, &Np, &Np, out_file, out_comments);

/* Output the bispectrum to a FITS file:
*/
       sprintf(out_file, "%s_b", out_prefix);
       sprintf(out_comments, "bispectrum with %s", parameter_fname);
       ndim = 3;
       JLP_D_WRITEIMAG(bisp2, &ngamma, &ndim, &ngamma, out_file, out_comments);

free(bisp2);
return(0);
}
/*************************************************************************
* Allocate memory for arrays and initialize them
*
*************************************************************************/
static int init_arrays(complex **psc, double **pwr, double **modsq, 
                       double **snrm, double **long_int, double **bisp1, 
                       int Ns, int Np, int ngamma)
{
complex u_zero;
int isize;
register int i;

*psc = (complex *)malloc(Ns * Ns * sizeof(complex));
u_zero = float_to_cplx(0., 0.); 
for(i = 0; i < Ns * Ns; i++) (*psc)[i] = u_zero;

isize = Np * Np;
*pwr = (double *)malloc(isize * sizeof(double));
*modsq = (double *)malloc(isize * sizeof(double));
*snrm = (double *)malloc(isize * sizeof(double));
*long_int = (double *)malloc(isize * sizeof(double));

for(i = 0; i < isize; i++) {
    (*pwr)[i] = 0.; 
    (*modsq)[i] = 0.; 
    (*snrm)[i] = 0.; 
    (*long_int)[i] = 0.; 
    }

*bisp1 = (double *)malloc(4 * ngamma * sizeof(double));
for(i = 0; i < 4 * ngamma; i++) (*bisp1)[i] = 0.;

return(0);
}
/*************************************************************************
* Free memory that was allocated in init_arrays() 
*
*************************************************************************/
static int free_arrays(complex *psc, double *pwr, double *modsq, 
                       double *snrm, double *long_int, double *bisp1) 
{
free(psc);
free(pwr);
free(modsq);
free(snrm);
free(long_int);
free(bisp1);
return(0);
}
/***********************************************************************
* compute_scale_and_resolution
* Compute and display various parameters for a diagnostic...
*
* INPUT:
* Ns: size of phase screen 
* Np: size of elementary frames (=Ns/2)
* L0, lo: external and internal scale lengths of turbulence
* lambda_cent: central wavelength (m),
* fov: field of view (in arcseconds)
* r_zero: Fried's parameter (m) for the central wavelength
*
* OUTPUT:
* LCn2: integrated CN2 corresponding to r_zero and wavelength
* Ls: width of phase screen (in m) 
* Scale_s: resolution of phase screen (in m/pixel)
* image_scale: resolution of image (in arcsecond/pixel)
***********************************************************************/
static int compute_scale_and_resolution(float lambda_cent, float r_zero,
                                        int Np, int Ns, float fov, float *LCn2, 
                                        float *Ls, float *Scale_s, 
                                        float *image_scale)
{
float kk_cent, r0_pix, fwhm_pix;

/*  kk_cent: wave number corresponding to lambda_cent 
*/
kk_cent = 2.0 * PI / lambda_cent;
 
/* Original version with LCn2 = 5.0e-13
* r_0 = (0.423 kk^2 LCn2)^{-3/5}
*/
*LCn2 = pow((double)r_zero, 3./5.) / (0.423 * kk_cent * kk_cent);
printf("scale_and_resolution routine:\n");
printf(" lambda_cent= %e r_zero=%.3f LCn2=%e\n", lambda_cent, r_zero, *LCn2);

/* Ls: width of phase screen (in m) */
/* rapsec = pi/(180. * 3600.0)
* Scale of image in pixels/arcsecond: Scale = (Ls/lambda_cent)*(Np/Ns)*rapsec
*/
   *image_scale = (float)Np / fov;
/* Field of view = Np/Scale
* Size of phase screen and resolution (m/pix)
*/
  *Ls = Ns / ((fov / lambda_cent) * RAD_PER_ARCSECOND);
  printf("Width of phase screen: Ls = %.3f m\n", *Ls);

  *Scale_s = *Ls / (float)Ns;
  printf("Resolution of phase screen: Scale_s = %.4f m/pixel\n", *Scale_s);

  r0_pix = r_zero / *Scale_s;
  printf("Hence r0 corresponds to %.3f pixels\n", r0_pix);

  printf(" Image scale: %.4f pixel/arcsecond \n", *image_scale);
  printf(" Check:  %.4f (should be the same...)\n", 
          (*Ls/lambda_cent)*((float)Np/(float)Ns)*RAD_PER_ARCSECOND);

/* rapsec = pi/(180. * 3600.0)
* and image_scale in pixel/arcsec
*/
  fwhm_pix = ((lambda_cent / r_zero) / RAD_PER_ARCSECOND) * (*image_scale);
  printf("Mean FWHM of long integration: %.3f pixels\n",fwhm_pix);

return(0);
}
/*
* Ls: size of phase screen (in m)
* lambda_cent: wavelength (in meters)
* MAXCOMP: maximum number of components in source 
* Ns: size of phase screen 
*/

