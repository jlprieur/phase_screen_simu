/*  simu2.c
 Creates simulated speckle data for a IDIMxIDIM array.
 From Jeff Morgan and Elliott (kiloho::jeff or kiloho::elliott)

 JLP
 Version 20-04-90

   The call is:
      simu2 [-n nframes] [-p nph] [-f ffname] [-o ofname] [-O Ofname]
            [-s seed] [-a] [-q]
    where the flags have the following functions:
	    -n  the program makes 'nframes' frames of simulated data
	    -p  the average number of photons per frame (default is 150) 
	    -f  the program uses the input image 'ffname' for filtering the
		frame
            -o  the program uses the input image 'ofname' as the object
		file (default is a pt source)
            -O  the output is written to 'Ofname' (default is
	        'simu2.data1')
            -s  'seed' is used as the seed for the random number
	         generator
	    -r  'fradius' of the top-hat used for filtering 
		the aperture function (linked with seeing)
            -a  aperture functions and psf's for frames 1 through 4
		outputted in images  'aper1.bdf' etc. and 
		'psf1.bdf' etc (default is no output)
            -q  frame information is not printed to the terminal 
		(default is to print)
     This program must be compiled with the libraries '-lzdump','-lnr'
     and '-lm'.

 It contains calls to an FFT routine called 'fourn',
 which is nearly a literal transcription of the routine of the same
 name found in the book "Numerical Recipes" by Press, et al.

     */
 
/*
#define DEBUG
*/

#ifdef ibm
#define FOURN1 fourn1
#else
#define FOURN1 fourn1_
#endif
#define IDIM 64
#define DIM 4096
#define MYRAN  2147483647.0
/* In previous version " random()/MYRAN " was used 
to generate pseudo-random numbers between 0 and 1 */
#define MAXPHOTONS 12000
#define MAXFRAMES 5000

#include <stdio.h>
#include <math.h>
#include <fcntl.h>         /* For file handling */
#include <jlp_ftoc.h>
 
 
/* Data storage arrays */
static float filter[DIM],object[DIM];
static float prf[DIM],psf[DIM];
static float filhatr[DIM],filhati[DIM];
static float frame[DIM],rhat[DIM];
static float nphotons[MAXFRAMES];
static int pacfrme[MAXPHOTONS*2+1];
 
/* Packed arrays for FFT routine */
static float psfpc[2*DIM], objpc[2*DIM], prfrm[2*DIM];
static float packhat[2*DIM], packfil[2*DIM], hat[2*DIM];
 
main(argc,argv)  
int argc;
char **argv;
{
/* input variables */
int nframes, nph, naper;
int iseed,quiet;
char ofname[41],Ofname[41],ffname[41];
int ffile,objectfile;
 
/* main frame loop index */
int ncount;
 
/* file access variables */
int fd1,n_written;

/* misc. variables */
int frvalue;
float nphr, sum, rq;
float fradius;
int sign1,sign2;
char *s;
int i, j, k, n, istat;
int tail[1], nn[2];

/* JLP variables */
int ndim = 2, ioption;
int nx1,ny1,idim1;
float w1, w2;
char comments[81];
 
 
 /* add a '\0' to prevent any accident... */
 comments[80]='\0'; ofname[41]='\0'; Ofname[41]='\0'; ffname[41]='\0';

    /* Set default command line parameters */
   nframes = 1;                  /* # of frames to be made = 1       */
   nph = 150;                    /* Ave. # of photons per frame = 150*/
   ffile = 0;                    /* No filter file default is a      */
				 /* circular top hat (r=8) filter.  */
   fradius = 8.;                  /* Radius of the top hat circular filter */
   strcpy(Ofname,"simu2.data1");   /* output written to simu2.data1      */
   objectfile=0;                 /* No object file,default is pt     */
                                 /* source.                          */
   quiet = 0;                    /* print out frame info             */
   iseed = 1;                    /* seed for rand. # gen set to 1    */
   naper = 0;                    /* no aperture fn or psf's outputted*/
/*   naper = 1; */                    /* aperture and psf outputted */   
 
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
				printf(" Error: Max number of frames is %d\n",
				nframes);
				}
                                --argc;
				break;
                         case 'p':
				sscanf(*++argv,"%d",&nph);
				if(nph > MAXPHOTONS*0.9)
				{ nph = MAXPHOTONS*0.9;
				printf(" Error: Max number of photons is %d\n",
				nph);
				}
                                --argc;
				break;
                         case 'r':
				sscanf(*++argv,"%f",&fradius);
                                --argc;
				break;
                         case 'f':
				strcpy(ffname,*++argv);
                                --argc;
				ffile=1;
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
				printf("simu2 [-n nframes]");
				printf(" [-p nph] [-f ffname] [-o ofname]");
				printf(" [-r filter_radius] [-O Ofname]");
				printf(" [-s seed] [-a] [-q]\n");
                                exit(-1);
				break;
                         }

/* *********************************************************** */ 
nn[0]=IDIM; nn[1]=IDIM;
/* JLP format: */
JLP_BEGIN();
JLP_INQUIFMT();
 
/* Creation of the output data file (photon-counting mode) */
fd1=creat(Ofname,0666);

/* Filter: create filthati and filthatr */
istat = get_filter(ffile,ffname,fradius,naper);

/* object :  read object file, Fourier transform and output in objpc. */
istat = get_object(objectfile,ofname,naper);

JLP_RANDOM_INIT(&iseed);

/************************************************************ */ 
for(ncount=1;ncount<=nframes;ncount++)  {  /* Main loop starts here. */

 get_psf(naper,ncount);
 
/* ******************************************************************** */
/* PSF /

/* psf are packed for fourn routine. */
   for(i=0;i<DIM;i++)  {  
       psfpc[i*2]=psf[i];     
       psfpc[i*2+1]=0;
       }
#ifdef FOURN1
   ioption=1; FOURN1(&psfpc[0],&nn[0],&ndim,&ioption);
#else
   fourn(psfpc-1,nn-1,2,1); 
#endif
 
/* The PSF and the object are convolved via the convolution theorem. */
   for(i=0;i<DIM;i++) {
      prfrm[i*2]=psfpc[i*2]*objpc[i*2]-psfpc[i*2+1]*objpc[i*2+1];
      prfrm[i*2+1]=psfpc[i*2]*objpc[i*2+1]+psfpc[i*2+1]*objpc[i*2];
    }

/* Inverse FFT : */
#ifdef FOURN1
    ioption=-1; FOURN1(&prfrm[0],&nn[0],&ndim,&ioption);
#else
    fourn(prfrm-1,nn-1,2,-1); 
#endif
 
/* Unpacking : */ 
     for(i=0;i<DIM;i++)  {
        prf[i]=prfrm[i*2];
        }
 
/* Shifting center of array again */
   for(j=0; j<IDIM/2; j++)   
   {for(i=0; i<IDIM/2; i++)   
      {w1=prf[i+IDIM*j];
       w2=prf[i+IDIM/2+IDIM*j];
       prf[i+IDIM*j]=prf[i+IDIM/2+IDIM*(j+IDIM/2)];
       prf[i+IDIM/2+IDIM*j]=prf[i+IDIM*(j+IDIM/2)];
       prf[i+IDIM/2+IDIM*(j+IDIM/2)]=w1;
       prf[i+IDIM*(j+IDIM/2)]=w2;
    }
   }    

/* Create obj_conv.bdf files */
if(naper==1 && ncount == 1) {
	strcpy(ffname,"obj_conv.bdf");
	nx1=IDIM; ny1=IDIM; idim1=IDIM;
	printf(" Output of %s \n",ffname);
	strcpy(comments," Object convolved with PSF -real- //");
        JLP_WRITEIMAG(prf,&nx1,&ny1,&idim1,ffname,comments);
    }

/* ****************************************************************** */

/* First normalise the image (prf) to nph photons (but still analogic) */
    nphr=nph;
    sum=0;  
    for(i=0; i<DIM; i++) sum=sum+prf[i];
    sum = sum/nphr;
    for(i=0; i<DIM; i++) prf[i]=prf[i]/sum;

    clipping(&nphr,iseed);

    if(quiet == 0)
        printf("# of photons in frame %d = %f \n",ncount,nphr); 
    nphotons[ncount-1]=nphr;
 
/* Conversion of frame information to pacfrme for writing
 to the output file.
 Writes only the position of non-zero pixels, and as many times
 as there are photons in that location */
    k=0;  
/* The first value is a negative integer with the number of the frame */
    pacfrme[0]=(-1)*ncount;  
    for(j=0; j<IDIM; j++)  {
	for(i=0; i<IDIM; i++)  {
             frvalue=frame[i+j*IDIM];
	     for(n=0; n<frvalue; n++)  {
		 pacfrme[k*2+1]=i;
                 pacfrme[k*2+2]=j;
                 k=k+1;
                 }
             }
         }
 
/* data written to output. */
    n_written=write(fd1,pacfrme,((int)nphr*2+1)*sizeof(int));
    
if(ncount==1)  {
/* Create #1's frame file */
	strcpy(ffname,"frame1.bdf");
	nx1=IDIM; ny1=IDIM; idim1=IDIM;
	printf(" Output of %s \n",ffname);
	sprintf(comments," First frame //                   ");
        JLP_WRITEIMAG(frame,&nx1,&ny1,&idim1,ffname,comments);
    }
 
if(ncount==nframes)  {
/* Create last frame file */
	strcpy(ffname,"endfrm.bdf");
	nx1=IDIM; ny1=IDIM; idim1=IDIM;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Last frame //                     ");
        JLP_WRITEIMAG(frame,&nx1,&ny1,&idim1,ffname,comments);
 
/* Create 1-d file containing the number of photons in each frame */
	strcpy(ffname,"nphotons.bdf");
	nx1=nframes; ny1=1; idim1=nframes;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Nphotons //                        ");
        JLP_WRITEIMAG(nphotons,&nx1,&ny1,&idim1,ffname,comments);
 
/*Write a negative integer as entry in output.*/
    tail[0]=(-1)*(nframes+1);  
    n_written=write(fd1,tail,sizeof(int));  
    }
/* end of the main loop */
}

/* *************************************************************** */
/* End : */

  if(quiet == 0) printf("simu2: End calculation.\n");

  JLP_END();
}
/* *********************************************************** */ 
/* Generate Fourier filter: */
int get_filter(ffile,ffname,fradius,naper)
char *ffname;
float fradius;
int ffile, naper;
{
float w1,w2;
register int i, j;
int nx1,ny1,idim1;
int ndim = 2, ioption, nn[2];
float r, d, sum;
char comments[81];

#ifdef DEBUG
printf(" Succesful entry in get_filter, naper = %d \n",naper);
#endif
    if(ffile == 0)  {
/*  circular top hat filter centered on (IDIM/2,IDIM/2). */
        for(j=0;j<IDIM;j++) 
        {   w1=(j-(IDIM/2))*(j-(IDIM/2)); 
	    for(i=0;i<IDIM;i++) 
	    {
            w2=(i-(IDIM/2))*(i-(IDIM/2)); 
	    d=sqrt(w1+w2);
	    if(d<fradius)
	        filter[i+j*IDIM]=1.;
	    else
	        filter[i+j*IDIM]=0.;
	    }
	}
       } 
    else 
       {idim1=IDIM;
          JLP_READIMAG(filter,&nx1,&ny1,&idim1,ffname,comments);
       }
 
/* Normalization: */
	sum = 0.;
        for(i=0; i<DIM; i++) sum = sum + filter[i];
        if(sum != 0.) {for(i=0; i<DIM; i++) filter[i] = filter[i]/sum;}

/*Shifts array so that [0,0] is at [IDIM/2,IDIM/2]. 
 This is convenient because the fourn routine centers the transform on
 [0,0] and in output images, we want the center of the image at 
 [IDIM/2,IDIM/2] (center of screen). */

   /*
   for(j=0; j<IDIM/2; j++)   
   {for(i=0; i<IDIM/2; i++)   
      {w1=filter[i+IDIM*j];
       w2=filter[i+IDIM/2+IDIM*j];
       filter[i+IDIM*j]=filter[i+IDIM/2+IDIM*(j+IDIM/2)];
       filter[i+IDIM/2+IDIM*j]=filter[i+IDIM*(j+IDIM/2)];
       filter[i+IDIM/2+IDIM*(j+IDIM/2)]=w1;
       filter[i+IDIM*(j+IDIM/2)]=w2;
    }
   }    
   */
 
/*Packing arrays for transform routine*/
        for(i=0;i<DIM;i++)  { 
          packfil[i*2]=filter[i];
	  packfil[i*2+1]=0.;
	  }
 
/* Taking transforms of the filter.           */
  nn[0]=IDIM; nn[1]=IDIM;
#ifdef FOURN1
  ioption = 1; FOURN1(&packfil[0],&nn[0],&ndim,&ioption);
#else
    fourn(packfil-1,nn-1,2,1); 
#endif

/* Unpacking transformed arrays. */
      for(i=0;i<DIM;i++)  { 
        filhatr[i]=packfil[i*2];
        filhati[i]=packfil[i*2+1];
	}

/* ****************** Output of filter: ************************* */
if(naper==1)  {
/* Create aper. #1's bdf file */
	strcpy(ffname,"filtfft.bdf");
	nx1=IDIM; ny1=IDIM; idim1=IDIM;
	printf(" Output of %s \n",ffname);
	sprintf(comments," FFT of the filter -real- radius: %f //",fradius);
        JLP_WRITEIMAG(filhatr,&nx1,&ny1,&idim1,ffname,comments);
        }

return(0);
}

/********************************************************************* */
/* object :  read object file, Fourier transform and output in objpc. */
/********************************************************************* */

int get_object(objectfile,ofname,naper)
char *ofname;
int objectfile, naper;
{
register int i, j;
int nx1,ny1,idim1;
int ndim = 2, ioption, nn[2];
float objsum, w1, w2;
char comments[81];

 if(objectfile==0)      
/* if an input file is not specified creates a point source in [IDIM/2,IDIM/2] */
       {for(i=0;i<DIM;i++) object[i]=0.;
       object[DIM/2]=1.;}
 else 
 {
	idim1=IDIM;
        JLP_READIMAG(object,&nx1,&ny1,&idim1,ofname,comments);
 }
   
/* Shifting center of array "object"*/
   for(j=0; j<IDIM/2; j++)   
   {for(i=0; i<IDIM/2; i++)   
      {w1=object[i+IDIM*j];
       w2=object[i+IDIM/2+IDIM*j];
       object[i+IDIM*j]=object[i+IDIM/2+IDIM*(j+IDIM/2)];
       object[i+IDIM/2+IDIM*j]=object[i+IDIM*(j+IDIM/2)];
       object[i+IDIM/2+IDIM*(j+IDIM/2)]=w1;
       object[i+IDIM*(j+IDIM/2)]=w2;
    }
   }    


/*
	strcpy(ofname,"object100.bdf                        ");
        JLP_WRITEIMAG(object,&nx1,&ny1,&idim1,ofname,comments);
    if(i>-100){JLP_END(); exit(-1);}
*/

/* "object" is normalized, packed and transformed.  */
   objsum=0;             
   for(i=0;i<DIM;i++)
       objsum=objsum+object[i];

   for(i=0;i<DIM;i++)  {    
       objpc[i*2]=object[i]/objsum;
       objpc[i*2+1]=0;
       }

   nn[0]=IDIM; nn[1]=IDIM;
#ifdef FOURN1
  ioption=1; FOURN1(&objpc[0],&nn[0],&ndim,&ioption);
#else
   fourn(objpc-1,nn-1,2,1); 
#endif

return(0);
}

/***************************************************************
* To generate PSF 
***************************************************************/

int get_psf(naper,ncount)
int naper, ncount;
{
int sign1,sign2;
int i, j, k, n, istat;
int nn[2];
int ndim = 2, ioption;
int nx1,ny1,idim1;
float d, w1, w2, w3, w4;
char ffname[41],comments[81];

nn[0]=IDIM; nn[1]=IDIM;
  
 /* Get complex numbers of unit magnitude and random phase. */
        for (i=0; i<DIM; i++)  { 

/* Real part */
           JLP_RANDOM(&w1);
	   sign1=(floor)(w1*2.-1.);
           sign1=sign1*2+1;
           JLP_RANDOM(&w3);
	   packhat[i*2]=w3*(float)sign1;

/* Imaginary  part */
           JLP_RANDOM(&w1);
	   sign2=(floor)(w1*2.-1.);
           sign2=sign2*2+1;
           w4=sqrt(1-w3*w3);
           packhat[i*2+1]=w4*(float)sign2;
           } 

#ifdef DEBUG
     printf(" packhat real, imag %f,%f \n",packhat[0],packhat[1]);
     printf(" sign1, sign2 %d,%d \n",sign1,sign2);
#endif

/* *********************************************************** */ 
/* Now taking transforms of the random array           */

#ifdef FOURN1
   ndim=2; ioption=1; FOURN1(&packhat[0],&nn[0],&ndim,&ioption);
#else
    fourn(packhat-1,nn-1,2,1);
#endif
 
#ifdef DEBUG
     printf(" packhat: real[0],imag[0] %f,%f \n",packhat[0],packhat[1]);
#endif

/* *********************************************************** */ 
/* Below, the transforms of the two arrays are
  multiplied together to get the transform of their convolution.  */
for(i=0;i<DIM;i++)  {
    hat[i*2]=packhat[i*2]*filhatr[i]-packhat[i*2+1]*filhati[i];
    hat[i*2+1]=packhat[i*2+1]*filhatr[i]+packhat[i*2]*filhati[i];
    }
 
/*The transf. of the conv. is inverse-transfd*/ 
#ifdef FOURN1
    ioption=-1; FOURN1(&hat[0],&nn[0],&ndim,&ioption);
#else
    fourn(hat-1,nn-1,2,-1); 
#endif
 
/* The convolution is masked (pupil of a prime mirror with secondary) */
/* "hat" here becomes the aperture function */
    for(j=0;j<IDIM;j++)  {  
       for(i=0;i<IDIM;i++)  {
          d=(j-(IDIM/2))*(j-(IDIM/2))+(i-(IDIM/2))*(i-(IDIM/2));
	  d=sqrt(d);
          if(d>(IDIM/2)) {                 
	      hat[(j * IDIM + i) * 2]=0.;  
	      hat[(j * IDIM + i) * 2 + 1]=0.;  
	      }
	  /* Hole in the center:
          else if(d<7.0)  {
	      hat[(j * IDIM + i) * 2]=0.;  
	      hat[(j * IDIM + i) * 2 + 1]=0.;  
	      }
	  */
        }
    }
 
if(naper==1 && ncount<3)  {

       for(k=0;k<DIM;k++) rhat[k]=hat[k*2];

/* Create aper. #1's bdf file */
	sprintf(ffname,"aper%d.bdf",ncount);
	nx1=IDIM; ny1=IDIM; idim1=IDIM;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Aperture function -real part- #%d //",ncount);
        JLP_WRITEIMAG(rhat,&nx1,&ny1,&idim1,ffname,comments);
    }
 
#ifdef FOURN1
   ioption=1; FOURN1(&hat[0],&nn[0],&ndim,&ioption);
#else
   fourn(hat-1,nn-1,2,1); 
#endif
/* Now "hat" is the FFT of the aperture function, so the squared modulus
of "hat" is the point spread function. */
   for(i=0;i<DIM;i++)
       psf[i]=hat[i*2]*hat[i*2]+hat[i*2+1]*hat[i*2+1];
 
/* Create psf*.bdf files */
if(naper==1 && ncount<3)  {
	sprintf(ffname,"psf%d.bdf",ncount);
	nx1=IDIM; ny1=IDIM; idim1=IDIM;
	printf(" Output of %s \n",ffname);
	sprintf(comments," PSF -modsq- #%d //",ncount);
        JLP_WRITEIMAG(psf,&nx1,&ny1,&idim1,ffname,comments);
	}
}
/*************************************************************
* Clipping routine to generate simulated photo-counting "frame" 
* simply look for the right number of photons to fit the
* total number of photons (which is random around a mean value)
* Not satisfactory... (JLP 10-09-1991)
* Old version, from Jeff Morgan
*************************************************************/ 
int clipping1(nphr,iseed)
float *nphr;
int iseed;
{
float w1, w2, fudge, sum, np;
register int i;
int k, sign1, sign2;

/* Now compute np (ie number of photons of current frame) */
    JLP_RANDOM(&w1);
    sign1=(floor)(w1*2-1);
    sign1=sign1*2+1;
    JLP_RANDOM(&w2);
    np=*nphr+(float)sign1*w2*sqrt(*nphr*2.2);
    for(i=0;i<6;i++)  {
        JLP_RANDOM(&w1);
        sign2=(floor)(w1*2-1);
        sign2=sign2*2+1;
        JLP_RANDOM(&w2);
        np=np+(float)sign2*w2*sqrt(fabs(np))*.3;
	}

/* Before calling this routine, the sum of prf has been normalised to nphr*/
    sum = *nphr;

/* Computes the value of fudge such that there are roughly 
np photons in the frame. Proceeds by adjusting "fudge" */
    fudge=1.0;
    k=0;
    while(sum >= np+.22*sqrt(np) || sum <= np-.22*sqrt(np)) {
       sum=0;
       k=k+1;
       for(i=0;i<DIM;i++)  {
           prf[i]=prf[i]*fudge;
           frame[i]=floor(prf[i]+.5);
           sum=sum+frame[i];
           }
       fudge=sqrt(np/(sum+1));
#ifdef DEBUG
       printf(" k = %d ",k);
#endif
       if(k>40)
	   fudge=sqrt(fudge);
       if(k>80)
	   fudge=sqrt(fudge);
       if(k>120)
	   fudge=sqrt(fudge);
       if(k>200)  {
	   printf("simu2/Fatal error: Clipping routine did not converge.\n");
	   *nphr = 0.;
	   exit(-1);
	   }
    }
#ifdef DEBUG
 printf(" np, sum %f %f\n",np,sum);
#endif
 *nphr = sum;
 return (0);
}
/*************************************************************
* Clipping routine to generate simulated photo-counting "frame" 
* with Poisson noise
*
* JLP version 10-09-91
*************************************************************/ 
int clipping(nphr,iseed)
float *nphr;
int iseed;
{
float sum, poidev();
int idum;
register int i;

/* Before calling this routine, the sum of prf has been normalised to nphr*/
    sum = *nphr;

       sum=0;
       for(i=0;i<DIM;i++)  {
/* Poisson generator */ 
           frame[i]=(int)poidev(prf[i],&iseed);
/* Computing total number of photons: */
           sum=sum+frame[i];
           }
 printf(" sum %f \n",sum);
 *nphr = sum;
 return (0);
}
