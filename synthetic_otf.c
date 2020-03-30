/********************************************************************
* synthetic_otf.c
*
* To compute the synthetic Optical Transfer Function 
* and the synthetic Bispectral Transfer Function of a telescope
*
* Calls: 
* covera_mask, 
* (in jlp_cover_mask.c)
*  
* JLP
* Version 10-03-94
*******************************************************************/
#include <stdio.h>
#include <math.h>
#include <jlp_ftoc.h>

#define DEBUG

main()
{
float *mask, *modsq, *bisp;
int nx_mask, ny_mask, ir, max_nclosure, nbeta, ngamma;
int nx_modsq, ny_modsq, nx_bisp, ny_bisp; 
int pntr1, isize, iopt;
register int i;
char mask_name[61], mask_comments[81], buffer[81];
char modsq_name[61], modsq_comments[81];
char bisp_name[61], bisp_comments[81];

printf(" bisp_to_image           -- Version 08-03-94 --\n");
printf(" Computation of the Optical Transfer Function \n");
printf(" and Bispectral Transfer Function of a telescope \n");

printf(" Enter your choice =< ");
gets(buffer);sscanf(buffer,"%d",&iopt);
printf(" OK: iopt=%d\n",iopt);

/***************************************************************/
JLP_BEGIN();
JLP_INQUIFMT();

  printf(" Parameters for the bispectrum: \n");
  printf(" Radius of uv-coverage (IR) in pixels and max_nclosure:\n"); 
  gets(buffer); sscanf(buffer,"%d,%d",&ir,&max_nclosure);
#ifdef DEBUG
  printf(" OK, ir=%d max_nclosure=%d \n",ir,max_nclosure);
#endif

/* Read mask file: */
  printf(" Input mask (u-v coverage) (0 if no mask) =< "); gets(mask_name);
  if(mask_name[0] == '0')
    {
      nx_mask = 2 * (ir + 1);
      ny_mask = nx_mask;
      isize = nx_mask * ny_mask * sizeof(float);
      JLP_GVM(&mask,&isize);

/* Create a filled mask: */
      for(i = 0; i < nx_mask * ny_mask; i++) mask[i] = 1.;
    }
  else
    {
     JLP_VM_READIMAG(&pntr1,&nx_mask,&ny_mask,mask_name,mask_comments);
     JLP_FROM_MADRID(&pntr1,&mask);
    }

/* Compute u-v coverage with the mask and within a disk of radius ir: */
  covera_mask(mask,&nx_mask,&ny_mask,&ir,&max_nclosure,&nbeta,&ngamma);

/* Free memory: */
  JLP_FVM(&mask);

/* Read square modulus file */
  printf(" Output square modulus (OTF) =< "); gets(modsq_name);
  printf(" Output size of OTF  (square image)  nx =< "); 
  gets(buffer); scanf(buffer,"%d",&nx_modsq);
/* Get memory space: */
  ny_modsq = nx_modsq;
  isize = nx_modsq * ny_modsq * sizeof(float);
  JLP_GVM(&modsq,&isize);

/* Check size is correct: */
  if(nx_modsq < 2 * ir || ny_modsq < 2 * ir)
    {
    printf("Fatal error: modulus file is too small, does not fit u-v coverage\n");
    exit(-1);
    } 

/* Read bispectrum file */
  printf(" Output bispectrum =< "); gets(bisp_name);
  ny_bisp = ngamma; ny_bisp = 3;
/* Get memory space: */
  isize = nx_bisp * ny_bisp * sizeof(float);
  JLP_GVM(&bisp,&isize);

/* Now compute bispectrum slice: 
*/
   bisp_to_2D_image(modsq,&nx_modsq,&ny_modsq,
                    bisp,&nx_bisp,&nbeta,&ngamma,&iopt);

/* Output image: */
/*
  sprintf(modsq_comments,"D=%d from %s and %s\n",iopt,modsq_name,bisp_name);
*/
  JLP_WRITEIMAG(modsq,&nx_modsq,&ny_modsq,&nx_modsq,modsq_name,modsq_comments);
  JLP_WRITEIMAG(&bisp,&nx_bisp,&ny_bisp,&ny_bisp,bisp_name,bisp_comments);

/* Free memory: */
  JLP_FVM(&modsq);
  JLP_FVM(&bisp);

exit(0);
}
