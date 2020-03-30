/*  simu2.c
 Creates simulated speckle data for a 64x64 array.

 JLP
 Version 20-04-90

 From Jeff Morgan and Elliott (kiloho::jeff or kiloho::elliott)
 It contains calls to an FFT routine called 'fourn',
 which is nearly a literal transcription of the routine of the same
 name found in the book "Numerical Recipes" by Press, et al.

   The call is:
      simu2 [-n nframes] [-t nttag] [-p nph] [-f ffname] [-o ofname] [-O Ofname]
            [-s seed] [-a] [-q]
    where the flags have the following functions:
	    -n  the program makes 'nframes' frames of simulated data
	    -t  the number of time-tags per frame is set to 'nttag'
		(default is 1)
	    -p  the average number of photons per frame (default is 150) 
	    -f  the program uses the input image 'ffname' for filtering the
		frame
            -o  the program uses the input image 'ofname' as the object
		file (default is a pt source)
            -O  the output is written to 'Ofname' (default is
	        'simu2.data1')
            -s  'seed' is used as the seed for the random number
	         generator
            -a  aperture functions and psf's for frames 1 through 4
		outputted in images  'aper1.bdf' etc. and 
		'psf1.bdf' etc (default is no output)
            -q  frame information is not printed to the terminal 
		(default is to print)
     This program must be compiled with the libraries '-lzdump','-lnr'
     and '-lm'.
     */
 
#include <stdio.h>
#include <math.h>
#include <fcntl.h>         /* For file handling */
 
 
/* Data storage arrays */
float real1[4096],imag1[4096];
float filter[4096],object[4096];
float prf[4096],prfi[4096],tmp[4096];
float realhat[4096],imaghat[4096];
float filhatr[4096],filhati[4096];
float psf[4096];
int sign1[4096],sign2[4096];
float frame[4096],rhat[4096];
float nphotons[5000];
int pacfrme[8000];
 
/* Packed arrays for FFT routine */
float psfpc[8192],objpc[8192],prfrm[8192];
float pack[8192],packfil[8192];
float hat[8192];
 
/* Did I really use this many indices? */
int i,j,n,n1,n2,n3,n4,n5;
int n10,n11;
int n25;
int n30,n31,n32,n33;
int n41,n42,n43,n44,n50,n51;
 
/* input variables */
int nframes,nttag,nph,naper;
int iseed,quiet;
char ofname[41],Ofname[41],ffname[41];
int ffile,objectfile;
 
/* main frame loop index */
int ncount;
 
/* file access variables */
int fd1,n_written;
int id,ipts,irtn;
 
/* misc. variables */
int frmint,sum,ssum2;
float sum2,nphr,rq;
float fudge,sign3,sign4;
float d,r,objsum,fradius;
int nn[2];
char *s;
int tail[1];

/* JLP variables */
int ndim,ioption;
int nx1,ny1,idim1;
char comments[81];
 
main(argc,argv)  
int argc;
char **argv;
{
 /* add a '\0' to prevent any accident... */
 comments[80]='\0'; ofname[41]='\0'; Ofname[41]='\0'; ffname[41]='\0';

    /* Set default command line parameters */
   nframes = 1;                  /* # of frames to be made = 1       */
   nttag = 1;                    /* number of time-tags per frame=1  */
   nph = 150;                    /* Ave. # of photons per frame = 150*/
   ffile = 0;                    /* No filter file default is a      */
				 /* circular top hat (r=8) filter.  */
   fradius = 8;                  /* Radius of the top hat circular filter */
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
                                --argc;
				break;
                         case 't':
				sscanf(*++argv,"%d",&nttag);
                                --argc;
				break;
                         case 'p':
				sscanf(*++argv,"%d",&nph);
                                --argc;
				break;
                         case 'r':
				sscanf(*++argv,"%d",&fradius);
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
				printf("simu2 [-n nframes] [-t nttag]");
				printf(" [-p nph] [-f ffname] [-o ofname]");
				printf(" [-r filter_radius] [-O Ofname]");
				printf(" [-s seed] [-a] [-q]\n");
                                exit();
				break;
                         }

/* *********************************************************** */ 
/* JLP format: */
jlp_begin_();
jlp_inquifmt_();
 
/* Creation of the output data file (photon-counting mode) */
fd1=creat(Ofname,0666);

/* *********************************************************** */ 
for(ncount=1;ncount<=nframes;ncount++)  {  /* Main loop starts here. */

    for(n=0;n<=iseed+200*ncount;n++)
        rq=random();

 /* Get complex numbers of unit magnitude and random phase. */
    if(ncount==1)  {
        for (n1=0;n1<=4095;n1++)  { 
           imag1[n1]=random()/2147483647.0;
           real1[n1]=sqrt(1-imag1[n1]*imag1[n1]);
 
	   tmp[n1]=random()/2147483647.0*1.9999999;
	   sign1[n1]=tmp[n1];
	   sign1[n1]=sign1[n1]*2-1;
 
	   tmp[n1]=random()/2147483647.0*1.9999999;
	   sign2[n1]=tmp[n1];
	   sign2[n1]=sign2[n1]*2-1;
 
           imag1[n1]=imag1[n1]*sign1[n1];
	   real1[n1]=real1[n1]*sign2[n1];
 
           } 
     }
 
/* *********************************************************** */ 
     for(n42=0;n42<=63;n42++)  {
        for(n41=64/nttag;n41<=63;n41++)  {
	    real1[n42*64+n41-64/nttag]=real1[n42*64+n41];
	    imag1[n42*64+n41-64/nttag]=imag1[n42*64+n41];
	    }
        for(n43=0;n43<=64/nttag-1;n43++)  {
	   n44=n42*64+63-n43;
           imag1[n44]=random()/2147483647.0;
           real1[n44]=sqrt(1-imag1[n44]*imag1[n44]);
 
	   tmp[n44]=random()/2147483647.0*1.9999999;
	   sign1[n44]=tmp[n44];
	   sign1[n44]=sign1[n44]*2-1;
 
	   tmp[n44]=random()/2147483647.0*1.9999999;
	   sign2[n44]=tmp[n44];
	   sign2[n44]=sign2[n44]*2-1;
 
           imag1[n44]=imag1[n44]*sign1[n44];
	   real1[n44]=real1[n44]*sign2[n44];
	   }
        }

 
/* *********************************************************** */ 
/* Filter : */
if(ncount==1)  {
    if(ffile == 0)  {
/*  circular top hat filter centered on (32,32). */
        for(j=0;j<=63;j++) 
        { for(i=0;i<=63;i++) 
	  d=sqrt((j-32)*(j-32)+(i-32)*(i-32))
	  if(d<=fradius)
	  filter[i+j*64]=1.;
	  else
	  filter[i+j*64]=0.;
	}
       } 
    else 
/*     id=kzopen(ffname,0); ipts=kzread(id,filter,4096); irtn=kzclos(id); */
       { idim1=64;
          jlp_readimag_(filter,&nx1,&ny1,&idim1,ffname,comments);
       }
 
/*Shifts array so that [0,0] is at [32,32]. 
 This is convenient because the fourn routine centers the transform on
 [0,0] and in output images, we want the center of the image at 
 [32,32] (center of screen). */

/* The first step is to swap the lower and upper half lines : */
        for(i=0;i<=2047;i++)    
            tmp[i]=filter[i+2048];
        for(i=2048;i<=4095;i++) 
            tmp[i]=filter[i-2048];
/* then swap the lower an upper columns : */
        for(i=0;i<=4095;i++)  { 
/* r determines if the current pixel is in the first half or second half : */
            r=((float)i)/64.0;               
            if(r-i/64<.5)           
  	        filter[i]=tmp[i+32];
            else 
  	        filter[i]=tmp[i-32];
            } 
 
/*Packing arrays for transform routine*/
        for(i=0;i<=4095;i++)  { 
          packfil[i*2]=filter[i];
	  packfil[i*2+1]=0.;
	  }
 
/* Taking transforms of the filter.           */
        nn[0]=64; nn[1]=64;
        fourn(packfil-1,nn-1,2,1); 
/*      ioption=1; fourn1_(&packfil[0],&nn[0],&ndim,&ioption); */

/* Unpacking transformed arrays. */
      for(i=0;i<=4095;i++)  { 
        filhatr[i]=packfil[i*2];
        filhati[i]=packfil[i*2+1];
	}
    }

/* *********************************************************** */ 
/* Now taking transforms of the random array           */

/*Packing arrays for transform routine*/
        for(i=0;i<=4095;i++)  {  
            pack[i*2]=real1[i];
	    pack[i*2+1]=imag1[i];
	    }
    nn[0]=64; nn[1]=64;
    fourn(pack-1,nn-1,2,1);
/*    ndim=2; ioption=1; fourn1_(&pack[0],&nn[0],&ndim,&ioption); */
 
/* Unpacking transformed arrays. */
    for(i=0;i<=4095;i++)  { 
       realhat[i]=pack[i*2];
       imaghat[i]=pack[i*2+1];
       }
 
/* *********************************************************** */ 
/* Below, the transforms of the two arrays are
  multiplied together to get the transform of their convolution.  */
for(i=0;i<=4095;i++)  {
    hat[i*2]=realhat[i]*filhatr[i]-imaghat[i]*filhati[i];
    hat[i*2+1]=imaghat[i]*filhatr[i]+realhat[i]*filhati[i];
    }
 
/*The transf. of the conv. is inverse-transfd*/ 
    fourn(hat-1,nn-1,2,-1); 
    /* ioption=-1; fourn1_(&hat[0],&nn[0],&ndim,&ioption); */ 
 
/* The convolution is masked (pupil of a prime mirror with secondary) */
/* "hat" here becomes the aperture function */
    for(n10=0;n10<=63;n10++)  {  
       for(n11=0;n11<=63;n11++)  {
          d=(n10-32)*(n10-32)+(n11-32)*(n11-32);
	  d=sqrt(d);
          if(d>32.0) {                 
	      hat[n10*128+n11*2]=0.;  
	      hat[n10*128+n11*2+1]=0.;
	      }
          else if(d<7.0)  {
	      hat[n10*128+n11*2]=0.;
	      hat[n10*128+n11*2+1]=0.;
	      }
        }
    }
     for(i=0;i<=4095;i++)
	 rhat[i]=hat[i*2];
 
if(naper==1)  {
    if(ncount<=3)  {
/* Create aper. #1's bdf file */
	sprintf(ffname,"aper%d.bdf",ncount);
	nx1=64; ny1=64; idim1=64;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Aperture function (real) #%d //",ncount);
        jlp_writeimag_(rhat,&nx1,&ny1,&idim1,ffname,comments);
        }
    }
 
   fourn(hat-1,nn-1,2,1); 
   /* ioption=1; fourn1_(&hat[0],&nn[0],&ndim,&ioption);  */
/* Now "hat" is the FFT of the aperture function, so the modulus
squared of "hat" is the point spread function. */
   for(i=0;i<=4095;i++)
       psf[i]=hat[i*2]*hat[i*2]+hat[i*2+1]*hat[i*2+1];
 
/* ******************************************************************** */
/* object : */

if(ncount==1)  {
    if(objectfile==0)      
/* if an input file is not specified creates a point source in [32,32] */
       {for(i=0;i<=4095;i++) object[i]=0.;
       object[2080]=1.;}
    else 
/*   id=kzopen(ofname,0); ipts=kzread(id,object,4096); irtn=kzclos(id); */
	{idim1=64;
        jlp_readimag_(object,&nx1,&ny1,&idim1,ofname,comments);
	}
   
/* Shifting center of array "object"*/
   for(i=0;i<=2047;i++)  
       tmp[i]=object[i+2048];
   for(i=2048;i<=4095;i++)
       tmp[i]=object[i-2048];
   for(i=0;i<=4095;i++)  {
       r=((float)i)/64.0;
       if(r-i/64<.5)
	  object[i]=tmp[i+32];
       else 
	  object[i]=tmp[i-32];
        } 

/* "object" is normalized, packed and transformed.  */
   objsum=0;             
   for(i=0;i<=4095;i++)
       objsum=objsum+object[i];
   for(i=0;i<=4095;i++)  {    
       objpc[i*2]=object[i]/objsum;
       objpc[i*2+1]=0;
       }
   fourn(objpc-1,nn-1,2,1); 
/*   ioption=1; fourn1_(&objpc[0],&nn[0],&ndim,&ioption); */
   }
 
/* ******************************************************************** */
/* PSF /

/* psf are packed for fourn routine. */
   for(i=0;i<=4095;i++)  {  
       psfpc[i*2]=psf[i];     
       psfpc[i*2+1]=0;
       }
   fourn(psfpc-1,nn-1,2,1); 
/*   ioption=1; fourn1_(&psfpc[0],&nn[0],&ndim,&ioption); */ 
 
/* The PSF and the object are convolved via the convolution theorem. */
   for(i=0;i<=4095;i++) {
      prfrm[i*2]=psfpc[i*2]*objpc[i*2]-psfpc[i*2+1]*objpc[i*2+1];
      prfrm[i*2+1]=psfpc[i*2]*objpc[i*2+1]+psfpc[i*2+1]*objpc[i*2];
    }

/* Inverse FFT : */
    fourn(prfrm-1,nn-1,2,-1); 
/*    ioption=1; fourn1_(&prfrm[0],&nn[0],&ndim,&ioption); */
 
/* Unpacking : */ 
     for(i=0;i<=4095;i++)  {
        prf[i]=prfrm[i*2];
        prfi[i]=prfrm[i*2+1];
        }
 
/* Shifting center of array again */
   for(i=0;i<=2047;i++)   
       tmp[i]=prf[i+2048];
   for(i=2048;i<=4095;i++)
       tmp[i]=prf[i-2048];
   for(i=0;i<=4095;i++)  {
       r=((float)i)/64.0;               
       if(r-i/64<.5)
	  prf[i]=tmp[i+32];
       else 
	  prf[i]=tmp[i-32];
        } 

/* Create psf. #1's bdf file */
if(naper==1) {
    if(ncount<=3)  {
	sprintf(ffname,"psf%d.bdf",ncount);
	nx1=64; ny1=64; idim1=64;
	printf(" Output of %s \n",ffname);
	sprintf(comments," PSF  #%d //",ncount);
        jlp_writeimag_(prf,&nx1,&ny1,&idim1,ffname,comments);
        }
    }
 
/* ****************************************************************** */
/* Clipping routine starts here. */

    nphr=nph*1.0;
    sum=0;  
    for(i=0;i<=4095;i++)
       sum=sum+prf[i];
    for(i=0;i<=4095;i++)  
        prf[i]=prf[i]/sum*nphr;
    fudge=1.0;
    sign3=floor(random()/2147483647.0*2)*2.0-1.0;
    nphr=nphr+sign3*random()/2147483647.0*sqrt(nphr*2.2);
    for(n51=0;n51<=5;n51++)  {
        sign4=floor(random()/2147483647.0*2)*2.0-1.0;
        nphr=nphr+sign4*random()/2147483647.0*sqrt(fabs(nphr))*.3;
	}
    n25=0;
    while(sum2>=nphr+.22*sqrt(nphr) || sum2<=nphr-.22*sqrt(nphr)) {
       sum2=0;
       n25=n25+1;
       for(i=0;i<=4095;i++)  {
           prf[i]=prf[i]*fudge;
           frame[i]=floor(prf[i]+.5);
           sum2=sum2+frame[i];
           }
       fudge=sqrt(nphr/(sum2+1));
       if(n25>40)
	   fudge=sqrt(fudge);
       if(n25>80)
	   fudge=sqrt(fudge);
       if(n25>120)
	   fudge=sqrt(fudge);
       if(n25>200)  {
	   printf("simu2: Clipping routine did not converge.\n");
	   exit();
	   }
    }

    if(quiet == 0)
        printf("# of photons in frame %d = %f \n",ncount,sum2); 
    nphotons[ncount-1]=sum2;
 
/* Conversion of frame information to pacfrme for writing
 to the output file.
 Writes only the position of non-zero pixels, and as many times
 as there are photons in that location */
    n33=0;  
/* The first value is a negative integer with the number of the frame */
    pacfrme[0]=(-1)*ncount;  
    for(n30=0;n30<=63;n30++)  {
	for(n31=0;n31<=63;n31++)  {
             frmint=frame[n30*64+n31];
	     for(n32=0;n32<=frmint-1;n32++)  {
		 pacfrme[n33*2+1]=n30;
                 pacfrme[n33*2+2]=n31;
                 n33=n33+1;
                 }
             }
         }
 
/* data written to output. */
    ssum2=sum2;
    n_written=write(fd1,pacfrme,(ssum2*2+1)*sizeof(int));
    
if(ncount==1)  {
/* Create #1's frame file */
	strcpy(ffname,"frame1.bdf");
	nx1=64; ny1=64; idim1=64;
	printf(" Output of %s \n",ffname);
	sprintf(comments," First frame //");
        jlp_writeimag_(frame,&nx1,&ny1,&idim1,ffname,comments);
    }
 
if(ncount==nframes)  {
/* Create last frame file */
	strcpy(ffname,"endfrm.bdf");
	nx1=64; ny1=64; idim1=64;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Last frame //");
        jlp_writeimag_(frame,&nx1,&ny1,&idim1,ffname,comments);
 
/* Create 1-d file containing the number of photons in each frame */
	strcpy(ffname,"nphotons.bdf");
	nx1=nframes; ny1=1; idim1=nframes;
	printf(" Output of %s \n",ffname);
	sprintf(comments," Nphotons //");
        jlp_writeimag_(nphotons,&nx1,&ny1,&idim1,ffname,comments);
 
/*Write a negative integer as entry in output.*/
    tail[0]=(-1)*(nframes+1);  
    n_written=write(fd1,tail,sizeof(int));  
    }
/* end of the main loop */
}

/* *************************************************************** */
/* End : */

  if(quiet == 0) printf("simu2: End calculation.\n");

  jlp_end_();
}
