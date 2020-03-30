/******************************************************************************
*
* Program object_bin.c
* To create a synthetic binary star for simulations.
*
* OUTPUT:
*    Image file : output object (nx*ny) 
*
*
* From: object_cerga.c
* (G. Schumacher and  P. Cruzalebes, version 1990 )
*
* JLP
* Version 12-02-2008
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <jlp_ftoc.h>

#define PI 3.141592653589793
#define IDIM 128
	
float tab[IDIM * IDIM];
static int  nx, ny;
void gauss();

main(argc,argv)
char *argv[];
int argc;
{
   int i,j,k,l,err,nobj,rep,idim;
   float bmax, dx, dy, ratio, eten, xc, yc;
   char name[80], buffer[80], comments[81];

   JLP_BEGIN();

   printf("object_binary : to create synthetic objects (Version 12-02-2008\n");
   printf("Maximum output size: %d %d \n",IDIM,IDIM);

/* Prompt the format (Environment variable JLP_FORMAT)*/
   JLP_INQUIFMT();

   if (argc != 4 )
   {
      printf("USAGE:  object_bin output_size dx,dy,ratio output_file\n");
      printf("\n Enter output size  nx,ny :"); gets(buffer);
      sscanf(buffer,"%d,%d",&nx,&ny);
      printf("\n Enter dx,dy,ratio :"); gets(buffer);
      sscanf(buffer,"%f,%f,%f",&dx,&dy,&ratio);
      printf(" Enter output file name: "); gets(name);
   }
   else
   {
      sscanf(argv[1],"%d,%d",&nx,&ny);
      sscanf(argv[2],"%f,%f,%f",&dx,&dy,&ratio);
      strcpy(name,argv[3]);
   }

/* Check if values are correct: */
   if(nx < 0 || nx > IDIM || ny < 0 || ny > IDIM)
      { printf(" Fatal error: wrong input size: nx = %d ny = %d \n",nx,ny);
        exit(-1);
      }  

/* Reset object to zero: */
      for(j = 0; j < ny; j++)
         for(i = 0; i < nx; i++) tab[i + IDIM * j]=0.;
		
/* Center: */
   xc = (nx/2) - dx/2.;
   yc = (ny/2) - dy/2.;
/* Diametre a mi-hauteur en pixels */ 
   eten = 4.;
/* Brillance maximum au centre */ 
   bmax = 1.;
   gauss(xc,yc,eten,bmax);

/* Center: */
   xc = (nx/2) + dx/2.;
   yc = (ny/2) + dy/2.;
/* Diametre a mi-hauteur en pixels */ 
   eten = 4.; 
/* Brillance maximum au centre */ 
   bmax = ratio;
   gauss(xc,yc,eten,bmax);

/* Write output image: */
      sprintf(comments," Object_bin: %f,%f,%f",dx,dy,ratio);
      idim = IDIM;
      JLP_WRITEIMAG(tab,&nx,&ny,&idim,name,comments);

   JLP_END();
}
/**************************************************************
*
**************************************************************/
void gauss(xc,yc,eten,bmax)
float eten, xc, yc, bmax;
{
   register int i,j;
   int x2,y2;
   float l2,r2;

   printf("\n");

   l2 = (float) ( eten * eten ) / ( (-4.) * log (2.) );

   for (j = 0; j < ny; j++) 
   {
      y2 = (j-yc)*(j-yc);
      for (i = 0; i < nx; i++)
      {
         x2 = (i-xc)*(i-xc);     
         r2 = (float) (x2 + y2);
	 tab[i + IDIM * j] += bmax * ( exp ( (double) ( r2 / l2 ) ) ); 
      }
   }
		
}		
