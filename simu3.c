/*****************************************************************************
* simu3.c
*
* Creates simulated speckle data for a IDIMxIDIM array.
* From Jeff Morgan and Elliott (kiloho::jeff or kiloho::elliott)
* JLP95: add the possibility of no clipping 
*    (i.e. "continous frame", no Poisson clipping to simulate photon detection)
* 
*
* SYNTAX: 
*      simu3 [-n nframes] [-p nph] [-r r_O] [-o ofname] [-O Ofname]
*            [-s seed] [-a] [-q]
*    where the flags have the following functions:
*	    -n  the program makes 'nframes' frames of simulated data
*	    -p  the average number of photons per frame (default is 150, 
*               Enter a negative value for no photo_counting mode) 
*	    -r  Kolmogorov's radius in meters (linked with seeing).
*            -o  the program uses the input image 'ofname' as the object
*		file (default is a point source)
*            -O  the output is written to 'Ofname' (default is
*	        'simu2.data1')
*            -s  'seed' is used as the seed for the random number
*	         generator
*            -a  aperture functions and psf's for frames 1 through 4
*		outputted in images  'aper1.fits' etc. and 
*		'psf1.fits' etc (default is no output)
*            -q  frame information is not printed to the terminal 
*		(default is to print)
*
* It contains calls to an FFT routine called 'fourn',
* which is nearly a literal transcription of the routine of the same
* name found in the book "Numerical Recipes" by Press, et al.
*
* OUTPUT photo-counting file:
* - The first value is a negative integer with the number of the frame 
* - Then the coordinates ix0 iy0 ix1 iy1 ix2 iy2 ... 
*
*  Previous JLP's version: 16-06-95
* (then 2006)
*
*  JLP
*  Version 09/05/2008 
**********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fcntl.h>         /* For file handling */
#include <jlp_ftoc.h>
#include <kolmo.h>         /* Prototypes of routines defined in kolmo.c */
 
/*
#define DEBUG
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

/* In "../fft/fourn.c": 
*/
#ifndef FOURN1
void fourn(float *data, int *nn, int ndim, int isign);
#endif

/* In "poidev.c": */
float poidev(float xm, int *idum);

/* Size of frames:
*/
#define IDIM 128 
#define DIM (IDIM * IDIM) 
/*
*/
#define MAXPHOTONS 24000
#define MAXFRAMES 40000
 
/* Data storage arrays */
static float object[DIM];
static float prf[DIM], psf[DIM];
static float filhatr[DIM], filhati[DIM];
 
/* Packed arrays for FFT routine */
static float psfpc[2*DIM], objpc[2*DIM], prfrm[2*DIM];
static float packhat[2*DIM], packfil[2*DIM];
 
/* Prototypes of routines: */
static int get_filter(float *filter, char *ffname, float fradius, int naper,
                      int idim1);
static int get_object(int objectfile, char *ofname, int naper, int idim1);
static int get_psf(float *psf, float *work_frame, int idim, int naper, 
                   int iframe, float r0);
static int Poisson_clipping(float *in_array, float *out_array, int isize,
                            int *iseed, int nph, int *nphr);
static int compute_and_write_photon_coordinates(float *prf, float *work_frame, 
                     char *pacfrme, int isize, int idim, int *iseed, int nph, 
                     int *nphr, int fd1, int iframe, int quiet);

int main(int argc, char *argv[])  
{
char pacfrme[MAXPHOTONS*2+1];
float work_frame[DIM], nphotons[MAXFRAMES];
int nframes, nph, naper, iseed, quiet, objectfile;
int iframe, fd1, n_written, nphr, photo_counting_mode;
float r0;
char ofname[41],Ofname[41],ffname[41];
char *s, tail[1], comments[81];
int nn[2], nx1, ny1, idim1, isize, istat;
register int i;
/*
short int tail[1];
*/
#ifdef FOURN1
int ndim = 2, ioption; 
#endif
 
 
 /* add a '\0' to prevent any accident... */
 comments[80]='\0'; ofname[41]='\0'; Ofname[41]='\0'; ffname[41]='\0';

    /* Set default command line parameters */
   nframes = 1;                  /* # of frames to be made = 1        */
   nph = 150;                    /* Ave. # of photons per frame = 150 */
   photo_counting_mode = 1;      /* Photon clipping by default */
   r0 = 0.1;                     /* Kolmogorov's radius (10 cm)       */
   strcpy(Ofname,"simu2.data1"); /* output written to simu2.data1     */
   objectfile=0;                 /* No object file,default is pt      */
                                 /* source.                           */
   quiet = 0;                    /* print out frame info              */
   iseed = 102;                  /* seed for rand. # gen set to 102   */
   naper = 0;                    /* no aperture fn or psf's outputted */
/*   naper = 1; */               /* aperture and psf outputted        */   
 
/* Get command line arguments */
     while(argc-- > 1)
        if( **++argv == '-')
           for(s= *argv+1;*s != '\0';s++)
	    switch(*s)
	     {
		 case 'n':
			sscanf(*++argv,"%d",&nframes);
			if(nframes > MAXFRAMES)
			{ nframes = MAXFRAMES;
			printf(" Fatal error: Max number of frames is %d\n",
			nframes); exit(-1);
			}
                        --argc;
	         	break;
                 case 'p':
			sscanf(*++argv,"%d",&nph);
                        if(nph < 0) {
                        photo_counting_mode = 0;
                        nph *= -1;
                        } else {
                        photo_counting_mode = 1;
                        }
			if(nph > MAXPHOTONS*0.9)
			{ nph = MAXPHOTONS*0.9;
			printf(" Fatal error: Max number of photons is %d\n",
			nph); exit(-1);
			}
                        --argc;
			break;
                 case 'r':
			sscanf(*++argv,"%f",&r0);
                        --argc;
			break;
                 case 'o':
			strcpy(ofname,*++argv);
                        --argc;
			objectfile=1;
			break;
                 case 'O':
			strcpy(Ofname,*++argv);
                        --argc;
			break;
                 case 's':
			sscanf(*++argv,"%d",&iseed);
                        --argc;
		        break;
                 case 'q':
			quiet=1;
                        break;
                 case 'a':
			naper=1;
			break;
                 default:
			printf("Illegal option: %s\n",s);
			printf("simu3 [-n nframes]");
			printf(" [-p nph] [-o ofname]");
			printf(" [-r r0] [-O Ofname]");
			printf(" [-s seed] [-a] [-q]\n");
                        exit(-1);
		        break;
                 }

printf(" simu3.c: Version 01/06/2008. Outer scale in pupil plane is 2*D \n");
printf(" Max nber of frames: %d , max nber of photons: %d\n",MAXFRAMES,MAXPHOTONS);
printf(" r0=%f  nph=%d nframes=%d photo_counting_mode=%d\n", 
       r0, nph, nframes, photo_counting_mode);
printf(" iseed = %d  \n",iseed);
printf(" Image scale is lambda/(2*D) per pixel\n"); 
printf(" With lamdba=550nm, D=2m: 28.3mas/pixel; with 64x64 pixels: 1.82x1.82 arcsec\n"); 

/* *********************************************************** */ 
isize = DIM;
nn[0]=IDIM; nn[1]=IDIM;
nx1=IDIM; ny1=IDIM; idim1=IDIM;

/* JLP format: */
JLP_INQUIFMT();
 
/* object:  read object file, Fourier transform and output in objpc. */
istat = get_object(objectfile, ofname, naper, idim1);

JLP_RANDOM_INIT(&iseed);

/* Creation of the output data file (photon-counting mode) */
if(photo_counting_mode) fd1=creat(Ofname,0666);

/************************************************************ */ 
for(iframe=1; iframe<=nframes; iframe++)  {  /* Main loop starts here. */

get_psf(psf, work_frame, idim1, naper, iframe, r0);
 
/* ******************************************************************** */
/* PSF */

/* psf data are packed for fourn routine. */
   for(i = 0; i < isize; i++) {
       psfpc[i*2]=psf[i];     
       psfpc[i*2+1]=0;
       }
#ifdef FOURN1
   ioption=1; FOURN1(&psfpc[0],&nn[0],&ndim,&ioption);
#else
   fourn(psfpc-1,nn-1,2,1); 
#endif
 
/* The PSF and the object are convolved via the convolution theorem
(i.e. product of FFT's). */
   for(i = 0; i < isize; i++) {
      prfrm[i*2] = psfpc[i*2]*objpc[i*2] - psfpc[i*2+1]*objpc[i*2+1];
      prfrm[i*2+1] = psfpc[i*2]*objpc[i*2+1] + psfpc[i*2+1]*objpc[i*2];
    }

/* Inverse FFT : */
#ifdef FOURN1
    ioption=-1; FOURN1(&prfrm[0],&nn[0],&ndim,&ioption);
#else
    fourn(prfrm-1,nn-1,2,-1); 
#endif
 
/* Unpacking the real part of prfrm to prf: */ 
    for(i = 0; i < isize; i++) prf[i]=prfrm[i*2];
 
/* Shifts array so that [0,0] is at [nx1/2, ny1/2]. 
 This is convenient because the fourn routine centers the transform on
 [0,0] and in output images, we want the center of the image at 
 [nx1/2, ny1/2] (center of screen). */
   recent_real(prf,nx1,ny1,idim1);

/* Create obj_conv.fits files */
if(naper==1 && iframe == 1) 
       {strcpy(ffname,"obj_conv.fits");
	printf(" Output of %s \n",ffname);
	strcpy(comments," Object convolved with PSF -analogic- //");
        JLP_WRITEIMAG(prf,&nx1,&ny1,&idim1,ffname,comments);
       }

/* ****************************************************************** */
/* Photon clipping simulating photo-counting detection... */
   if(photo_counting_mode) {
      compute_and_write_photon_coordinates(prf, work_frame, pacfrme, isize, 
                    idim1, &iseed, nph, &nphr, fd1, iframe, quiet);
      nphotons[iframe-1] = nphr;
   }
    
 
/* Output #1's frame to FITS file */
   if(iframe==1)  {
      if(!photo_counting_mode) for(i=0; i<isize; i++)  work_frame[i]=prf[i];
      strcpy(ffname,"frame1.fits");
      printf(" Output of %s \n",ffname);
      sprintf(comments," First frame");
      JLP_WRITEIMAG(work_frame,&nx1,&ny1,&idim1,ffname,comments);
    }
}
/* end of the main loop */
 
/* Output last frame to FITS file */
  if(!photo_counting_mode) for(i=0; i<isize; i++)  work_frame[i]=prf[i];
  strcpy(ffname,"endfrm.fits");
  printf(" Output of %s \n",ffname);
  sprintf(comments," Last frame //                     ");
  JLP_WRITEIMAG(work_frame,&nx1,&ny1,&idim1,ffname,comments);
 
/* Create 1-d file containing the number of photons in each frame */
/* Only if photon clipping... */
   if(photo_counting_mode) {
	strcpy(ffname,"nphotons.fits");
	nx1=nframes; ny1=1; idim1=nframes;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Nphotons //                        ");
        JLP_WRITEIMAG(nphotons,&nx1,&ny1,&idim1,ffname,comments);
 
/* Write a negative integer as entry in output.*/
       tail[0]=(-1)*(nframes+1);  
/* Before:
       n_written=write(fd1,tail,sizeof(short));  
*/
       n_written=write(fd1,tail,1);  
       close(fd1);
    }

/* *************************************************************** */
/* End : */

  if(quiet == 0) printf("simu3: End calculation.\n");

return(0);
}
/************************************************************ 
* Generate Fourier filter: 
*************************************************************/ 
static int get_filter(float *filter, char *ffname, float fradius, int naper,
                      int idim1)
{
float w1,w2;
register int i, j;
int nx1, ny1, nn[2];
#ifdef FOURN1
int ndim = 2, ioption;
#endif
float d, sum;
char comments[81];

#ifdef DEBUG
printf(" Succesful entry in get_filter, naper = %d \n",naper);
#endif
/*  circular top hat filter centered on (idim1/2, idim1/2). */
        for(j = 0; j < idim1; j++) 
        {   w1=SQUARE(j-(idim1/2)); 
	    for(i = 0; i < idim1; i++) 
	    {
            w2=SQUARE(i-(idim1/2)); 
	    d=sqrt(w1+w2);
	    if(d<fradius)
	        filter[i + j * idim1]=1.;
	    else
	        filter[i + j * idim1]=0.;
	    }
	}
 
/* Normalization: */
	sum = 0.;
        for(i=0; i< idim1 * idim1; i++) sum = sum + filter[i];
        if(sum != 0.) 
            for(i=0; i< idim1 * idim1; i++) filter[i] = filter[i]/sum;


/* Packing arrays for Fourier transform routine*/
        for(i = 0; i < idim1 * idim1; i++)  { 
          packfil[i*2]=filter[i];
	  packfil[i*2+1]=0.;
	  }
 
/* Taking transforms of the filter.           */
  nn[0] = idim1; nn[1] = idim1;
#ifdef FOURN1
  ioption = 1; FOURN1(&packfil[0],&nn[0],&ndim,&ioption);
#else
    fourn(packfil-1,nn-1,2,1); 
#endif

/* Unpacking transformed arrays. */
      for(i = 0; i < idim1 * idim1; i++)  { 
        filhatr[i] = packfil[i*2];
        filhati[i] = packfil[i*2+1];
	}

/* ****************** Output of filter: ************************* */
if(naper==1)  {
/* Create aper. #1's fits file */
	strcpy(ffname,"filtfft.fits");
	nx1 = idim1; ny1 = idim1; 
	printf(" Output of %s \n",ffname);
	sprintf(comments," FFT of the filter -real- radius: %f //",fradius);
        JLP_WRITEIMAG(filhatr,&nx1,&ny1,&idim1,ffname,comments);
        }

return(0);
}

/********************************************************************* */
/* object :  read object file, Fourier transform and output in objpc. */
/********************************************************************* */
static int get_object(int objectfile, char *ofname, int naper, int idim1)
{
register int i;
int nx1, ny1, nn[2];
#ifdef FOURN1
int ndim = 2, ioption;
#endif
float objsum;
char comments[81];

 if(objectfile==0)      
/* if an input file is not specified creates a point source 
* in [idim1/2, idim1/2] */
       {
       for(i = 0; i < idim1; i++) object[i]=0.;
       nx1 = IDIM; ny1 = nx1; 
       object[nx1/2 + (ny1/2) * idim1] = 1.;
       }
 else 
 {
        JLP_READIMAG(object,&nx1,&ny1,&idim1,ofname,comments);
 }
   
/* Shifting center of array "object"*/
   recent_real(object,nx1,ny1,idim1);

/*
	strcpy(ofname,"object100.fits                        ");
        JLP_WRITEIMAG(object,&nx1,&ny1,&idim1,ofname,comments);
    if(i>-100){JLP_END(); exit(-1);}
*/

/* "object" is normalized, packed and transformed.  */
   objsum=0;             
   for(i = 0; i < idim1 * idim1; i++) objsum += object[i];

   for(i = 0; i < idim1 * idim1; i++)  {    
       objpc[i*2]=object[i]/objsum;
       objpc[i*2+1]=0;
       }

   nn[0] = idim1; nn[1] = idim1;
#ifdef FOURN1
  ioption=1; FOURN1(&objpc[0],&nn[0],&ndim,&ioption);
#else
   fourn(objpc-1,nn-1,2,1); 
#endif

return(0);
}

/***************************************************************
* To generate PSF 
*
* INPUT:
*  work_frame[]: work array used for saving the image of the PSF
*
***************************************************************/
static int get_psf(float *psf, float *work_frame, int idim1, int naper, 
                   int iframe, float r0)
{
register int i;
int nn[2];
#ifdef FOURN1
int ndim = 2, ioption;
#endif
int nx1, ny1;
float L0, cent_obscu, diam, norm_fourier;
char re_name[61], im_name[61], re_comments[81], im_comments[81];

nn[0] = idim1; nn[1] = idim1;
nx1 = idim1; ny1 = idim1; 
diam = 2.; cent_obscu = 0.;
/* Outer scale should be at least twice the telescope diameter (support of the
autocorrelation of the pupil) */
L0 = 2*diam; 
  
/* *********************************************************** */ 
/* Complex Gaussian random array: */
 get_gaussft(packhat, nx1, ny1, idim1);

#ifdef DEBUG
     printf("Gaussian random array real, imag %f,%f\n",packhat[0],packhat[1]);
#endif

/* Multiply with square root of Kolmogorov spectrum */
 get_phift(packhat, nx1, ny1, idim1, L0, r0);

#ifdef DEBUG
     printf("Kolmogorov random array real, imag %f,%f\n",packhat[0],packhat[1]);
#endif

/* Shifting center of array */
recent_complex(packhat, nx1, ny1, idim1);

#ifdef DEBUG
     printf("Kolmogorov array (decentred) real, imag %f,%f\n",packhat[0],packhat[1]);
#endif

/* *********************************************************** */ 
/* Inverse FFT to compute the phase term: */
#ifdef FOURN1
   ioption=-1; FOURN1(&packhat[0],&nn[0],&ndim,&ioption);
#else
   fourn(packhat-1,nn-1,2,-1);
#endif

/* WARNING!!!!! 
 FOURN1 (and fourn) divide by nx*ny when ioption=-1 (inverse FFT).
*/
norm_fourier = sqrt((double)(nx1 * ny1));
for(i = 0; i < nx1*ny1*2; i++) packhat[i] *= norm_fourier;
 
#ifdef DEBUG
     printf(" packhat: real[0],imag[0] %f,%f \n",packhat[0],packhat[1]);
#endif

/*************************************************************/
/* We only keep the real part, and compute cos(phi) and sin(phi): */
   get_phasor(packhat, nx1, ny1, idim1);

/* Multiply with pupil: */
/* The convolution is masked (pupil of a prime mirror with secondary) */
/* "packhat" here becomes the aperture function */
   diam = (float)(nx1/2);
   mask_pupil(packhat, nx1, ny1, idim1, diam, cent_obscu);

/* Output pupil: */
/* Create aper. #1's fits file */
if(naper==1 && iframe<3)
  {
  for(i = 0; i < 80; i++) re_comments[i]=' ';
  re_comments[80]='\0';
  sprintf(re_name,"aper%d_re.fits",iframe);
  sprintf(re_comments," Pupil function -real part- #%d //",iframe);

  for(i = 0; i < 80; i++) im_comments[i]=' ';
  im_comments[80]='\0';
  sprintf(im_name,"aper%d_im.fits",iframe);
  sprintf(re_comments," Pupil function -imaginary part- #%d //",iframe);

  output_pupil(packhat,nx1,ny1,idim1,re_name,im_name,re_comments,im_comments);
  }
 

/* Shifting center of array */
   recent_complex(packhat,nx1,ny1,idim1);

/* FFT to go from pupil to image plane: */
#ifdef FOURN1
   ioption=1; FOURN1(&packhat[0],&nn[0],&ndim,&ioption);
#else
   fourn(packhat-1,nn-1,2,1);
#endif


/* *********************************************************** */ 
/* Now "packhat" is the FFT of the aperture function, so the squared modulus
of "packhat" is the point spread function. */
/* From field to intensity, and normalisation to 1 */
   get_intensity(packhat,psf,nx1,ny1,idim1,1);

return(0);
}
/*************************************************************
* Clipping routine to generate simulated photo-counting "frame" 
* with Poisson noise
* Calls poidev() to generate Poisson noise,
* from Numerical recipees Chap 7, in "poidev.c".
*
* WARNING: before calling poidev, 
*          the total flux should be normalised to nphot!
*
* INPUT:
* in_array: continuous image
* nph: theoretical number of photons of the image 
* iseed: value of the seed for random generation
*
* OUTPUT:
* out_array: photon-counting simulated image
* nphr: actual number of photons of the simulated image 
* iseed: modified value of the seed for random generation
*
* JLP version 10-09-91
*************************************************************/ 
static int Poisson_clipping(float *in_array, float *out_array, int isize,
                            int *iseed, int nph, int *nphr)
{
float sum;
register int i;

/* First normalise the image (prf) to nph photons (but still analogic) */
    sum=0; 
    for(i = 0; i < isize; i++) sum += in_array[i];
    sum /= (float)nph;
    for(i = 0; i < isize; i++) in_array[i] /= sum;

#ifdef DEBUG
    printf(" # of photons before clipping: %f \n",sum);
#endif

/* Before calling poidev, the sum of prf has to be normalised to nphr
*/
   *nphr = 0;
   for(i = 0; i < isize; i++)  {
/* Poisson generator */ 
       out_array[i]=(int)poidev(in_array[i], iseed);
/* Computing total number of photons: */
       *nphr += out_array[i];
    }
#ifdef DEBUG
 printf("End of clipping: sum=%d photons\n", *nphr);
#endif

return(0);
}
/***********************************************************************
*
* INPUT:
*  prf: analog simulated image
*  fd1: pointer of output file 
*  quiet: flag set to 0 for verbose mode
*               or to 1 for quiet mode
*  iframe: loop index (from 0 to nframes)
*  idim1: X size of prf and work_frame arrays 
*  isize: idim1 * idim1
*
* 
***********************************************************************/
static int compute_and_write_photon_coordinates(float *prf, float *work_frame, 
                     char *pacfrme, int isize, int idim1, int *iseed, int nph, 
                     int *nphr, int fd1, int iframe, int quiet)
{
float frvalue;
int n_written;
register int i, j, k, n;

/* ****************************************************************** */
/* Photon clipping simulating photo-counting detection... */
    Poisson_clipping(prf, work_frame, isize, iseed, nph, nphr);

    if(quiet == 0) {
/* modulus 50 ...*/
        if((iframe % 50) == 1)
           printf("simu3/Nber of photons in frame #%d = %d \n", iframe, *nphr); 
        }
 
/* Conversion of frame information to pacfrme for writing
 to the output file.
 Writes only the position of non-zero pixels, and as many times
 as there are photons in that location */
  k=0;  
  pacfrme[0] = (-1) * iframe;  
    for(j =0; j< idim1; j++)  {
	for(i = 0; i < idim1; i++)  {
             frvalue = work_frame[i + j * idim1];
	     for(n = 0; n < frvalue; n++)  {
		 pacfrme[k*2+1]=i;
                 pacfrme[k*2+2]=j;
                 k=k+1;
                 }
             }
         }
 
/* OUTPUT photo-counting file:
* - The first value is a negative integer with the number of the frame 
* - Then the coordinates ix0 iy0 ix1 iy1 ix2 iy2 ... 
*/
/*
     n_written=write(fd1,pacfrme,((*nphr)*2+1)*sizeof(short));
*/
     n_written = write(fd1, pacfrme, ((*nphr)*2+1));

return(0);
}
