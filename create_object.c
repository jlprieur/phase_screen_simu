/******************************************************************************
* Program create_object.c
* To create a synthetic object for simulations.
* Possibility of selecting standard functions or setting the values 
* of individual pixels.
*
* OUTPUT:
*    Image file : output object (nx*ny) 
*
*
* Derived from an old program 
* written by G. Schumacher and P. Cruzalebes (version of 1989)
*
* JLP
* Version 29-02-2008
******************************************************************************/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <jlp_ftoc.h>

#define DEBUG

static int interactive_menu(float *tab, int nx, int ny, int idim, int *ioption,
                            float *xc, float *yc, float *xw1, float *xw2,
                            float *yw1, float *yw2, float *max_level);
static int batch_menu(FILE *fp_batch, float *tab, int nx, int ny, int idim, 
                      int *ioption, float *xc, float *yc, float *xw1, 
                      float *xw2, float *yw1, float *yw2, float *max_level,
                      int italk);
static int add_one_object(float *tab, int nx, int ny, int idim,
                          int ioption, float xc, float yc, 
                          float xw1, float xw2, float yw1, float yw2,
                          float max_level);
static int pixel(float *tab, int nx, int ny, int idim, float xc, float yc, 
                 float level);
static int rectangle(float *tab, int nx, int ny, int idim, float xc, float yc, 
                     float xwidth, float ywidth, float max_level);
static int pyramide(float *tab, int nx, int ny, int idim, float xc, float yc, 
                    float xwidth, float ywidth, float max_level);
static int ellipse(float *tab, int nx, int ny, int idim, float xc, float yc, 
                  float xwidth, float ywidth, float max_level);
static int ellip_ring(float *tab, int nx, int ny, int idim, float xc, float yc, 
                      float xwidth1, float ywidth1, float xwidth2, 
                      float ywidth2, float max_level);
static int cone(float *tab, int nx, int ny, int idim, float xc, float yc, 
                float xwidth, float ywidth, float max_level);
static int auto_circ(float *tab, int nx, int ny, int idim, float xc, float yc, 
                     float xwidth, float ywidth, float max_level);
static int truncated_cosinus(float *tab, int nx, int ny, int idim, float xc, 
                             float yc, float xwidth, float ywidth, 
                             float max_level);
static int truncated_cosinus2(float *tab, int nx, int ny, int idim, float xc, 
                              float yc, float xwidth, float ywidth, 
                              float max_level);
static int truncated_sinc2(float *tab, int nx, int ny, int idim, float xc, 
                           float yc, float xwidth, float ywidth, 
                           float max_level);
static int truncated_airy(float *tab, int nx, int ny, int idim, float xc, 
                          float yc, float xwidth, float ywidth, 
                          float max_level);
static int dimmed_ellipse(float *tab, int nx, int ny, int idim, float xc, 
                          float yc, float xwidth, float ywidth, 
                          float max_level);
static int gauss(float *tab, int nx, int ny, int idim, float xc, float yc, 
                 float xwidth, float ywidth, float max_level);
static int sinc2(float *tab, int nx, int ny, int idim, float xc, float yc, 
                 float xwidth, float ywidth, float max_level);
static int airy(float *tab, int nx, int ny, int idim, float xc, float yc, 
                float xwidth, float ywidth, float max_level);
static int moffat(float *tab, int nx, int ny, int idim, float xc, float yc, 
                  float xwidth, float ywidth, float max_level);

int main(int argc, char *argv[])
{
int i, j, nx, ny, status, interactive_mode, ival, ioption, italk = 1;
float *tab, xc, yc, xw1, xw2, yw1, yw2, max_level;
char out_fname[80], out_comments[81], batch_fname[80], buffer[80];
FILE *fp_batch;

   JLP_BEGIN_C();

   printf("create_object: to create synthetic objects (Version 29-02-2008)\n");

/* Prompt the format (Environment variable JLP_FORMAT)*/
   JLP_INQUIFMT();

if(argc == 7 && *argv[3] != ' ' && *argv[3]) argc = 4;
if(argc == 7 && *argv[2] != ' ' && *argv[2]) argc = 3;
if(argc == 7 && *argv[1] != ' ' && *argv[1]) argc = 2;
switch (argc) {
  case 2:
    strcpy(batch_fname, argv[1]);
    interactive_mode = 0;
    if(italk) printf("Batch file = %s\n",batch_fname);
    if((fp_batch = fopen(batch_fname,"r")) == NULL) {
      fprintf(stderr, "Fatal error reading batch file %s\n", batch_fname);
      exit(-1);
      }
    break;
  case 1:
    interactive_mode = 1;
    if(italk) printf("No batch file: will work in interactive mode\n");
    break;
  default:
    printf("USAGE:\n");
    printf("  if batch mode :  create_object command_file\n");
    printf("  if interactive:  create_object \n");
    exit(-1);
  }

if(interactive_mode) {
  printf("\n Enter output size  nx,ny: "); 
  scanf("%d,%d", &nx, &ny);
  printf(" Enter output file name: "); 
  scanf("%s", out_fname);
  } else {
  buffer[0] = '#';
  while(!feof(fp_batch) && buffer[0] == '#') fgets(buffer,80,fp_batch);
  if(sscanf(buffer, "size=%d,%d", &nx, &ny) != 2){
    fprintf(stderr, "Fatal error: x y  size not found!\n");
    exit(-1);
    }

  buffer[0] = '#';
  while(!feof(fp_batch) && buffer[0] == '#') fgets(buffer,80,fp_batch);
  if((ival = sscanf(buffer, "outfile=%s", out_fname)) != 1) {
    fprintf(stderr, ">%s<\n", buffer);
    fprintf(stderr, "Fatal error: output filename not found! ival=%d\n", ival);
    exit(-1);
    }
  }
if(italk) printf("OK: outfile=%s nx=%d ny=%d\n", out_fname, nx, ny);

  if((tab = (float *)malloc(nx * ny * sizeof(float))) == NULL) {
    fprintf(stderr, "create_object/Fatal error allocating array with nx=%d ny=%d\n",
             nx, ny);
    exit(-1);
    }
/* Reset frame to zero: */
   for(j = 0; j < ny; j++)
      for(i = 0; i < nx; i++) tab[i + nx * j]=0.;
		
/* Prompts desired object function: */
  status = 0;
  while(!status) {
   if(interactive_mode){ 
   status = interactive_menu(tab, nx, ny, nx, &ioption, &xc, &yc, &xw1, &xw2,
                             &yw1, &yw2, &max_level);
   } else {
   status = batch_menu(fp_batch, tab, nx, ny, nx, &ioption, &xc, &yc, &xw1, &xw2,
                             &yw1, &yw2, &max_level, italk);
   }
   if(!status) status = add_one_object(tab, nx, ny, nx, ioption, xc, yc, 
                                       xw1, xw2, yw1, yw2, max_level);
   }

/* Write output image: */
 if(interactive_mode)
    strcpy(out_comments,"From create_object: interactive mode");
 else
    sprintf(out_comments,"From create_object: batch file=%.40s", batch_fname);

 JLP_WRITEIMAG(tab,&nx,&ny,&nx,out_fname,out_comments);

JLP_END();
free(tab);

return(0);
}
/******************************************************************************
* interactive_menu()
******************************************************************************/
static int interactive_menu(float *tab, int nx, int ny, int idim, int *ioption,
                            float *xc, float *yc, float *xw1, float *xw2,
                            float *yw1, float *yw2, float *max_level)
{
int status = 0;

*ioption = 0;
*xc = *yc = *xw1 = *xw2 = *yw1 = *yw2 = *max_level = 0;

   printf("*********** MENU: ****************\n");
   printf("******* I. Pixel function: ************\n");
   printf("- Setting the value of one pixel: ................ (1)\n");
   printf("******* II. Object with bounded support: ************\n");
   printf("   rectangle with uniform brightness.............. (2)\n");
   printf("   pyramide with rectangular support \n");
   printf("      and triangular brightness profile .......... (3)\n");
   printf("   ellipse with uniform brightness ............... (4)\n");
   printf("   elliptical ring with uniform brightness ....... (5)\n");
   printf("   ellipse with (truncated) following brightness function:\n");
   printf("     linear (cone) ............................... (6)\n");
   printf("     autocorrelation of a uniform disk ........... (7)\n");
   printf("     cosinus ..................................... (8)\n");
   printf("     cosinus^2 ................................... (9)\n");
   printf("     sinc2   (i.e. [sin(PI r)/ PI r]^2 ) ......... (10)\n");
   printf("     [J1(r)/r]^2 (central lobe of Airy function).. (11)\n");
   printf("   edge-dimmed ellipse ........................... (12)\n");
   printf("****** III. Object with unbounded support: ****\n");
   printf("   ellipse with following brightness function:\n");
   printf("     Gaussian .................................... (13)\n");
   printf("     sinc2   (i.e. [sin(PI r)/ PI r]^2 ) ......... (14)\n");
   printf("     [J1(r)/r]2 (Airy function) .................. (15)\n");
   printf("     Moffat profile .............................. (16)\n");
   printf("\n");
   printf("   Exit .........................................  (0)\n");
   printf("\n");
   printf("Enter your choice : ");
   scanf("%d", ioption );
   if(*ioption == 0) return(1);

  printf("Coordinates of the center (RELATIVE TO THE CENTER OF THE FRAME!)\n");
  scanf("%f,%f", xc, yc);

  printf("\n");

*max_level = 0.;

/* Pixel setting */
  if(*ioption == 1) {
    printf("Brightness level:\n");
    scanf("%f",  max_level);

/* Elliptical ring: */
  } else if(*ioption == 5) {
    printf("Values of the major axes: X min and max, Y min and max, and brightness level:\n");
    scanf("%f,%f,%f,%f,%f", xw1, xw2, yw1, yw2, max_level);

/* Autocorrelation of a disk: only a disk is allowed: */
  } else if(*ioption == 7) {
    printf("Diameter and brightness level:\n");
    scanf("%f,%f", xw1, max_level);
    yw1 = xw1;
/* Gaussian */
  } else if(*ioption == 13) {
    printf("X and Y diameters at half width in pixels and (max) brightness level:\n");
    scanf("%f,%f,%f", xw1, yw1, max_level);
/* sinc2 and Airy */
  } else if(*ioption == 14 || *ioption == 15) {
    printf("X and Y diameters of the first dark ring in pixels and (max) brightness level:\n");
    scanf("%f,%f,%f", xw1, yw1, max_level);
/* General case: */
  } else {
    printf("X and Y width in pixels and (max) brightness level:\n");
    scanf("%f,%f,%f", xw1, yw1, max_level);
  }

/* Handling of error: */
if(*max_level == 0) {
  printf("interactive_menu/Error max level is null, please enter a new value: (0 to exit) \n");
  scanf("%f", max_level);
  if(*max_level == 0) status = -1;
  }

return(status);
}
/******************************************************************************
* batch_menu()
*
* INPUT:
* italk: 1 if talkative (i.e. verbose) 
******************************************************************************/
static int batch_menu(FILE *fp_batch, float *tab, int nx, int ny, int idim, 
                      int *ioption, float *xc, float *yc, float *xw1, 
                      float *xw2, float *yw1, float *yw2, float *max_level,
                      int italk)
{
int status = 0;
char buffer[80];

*ioption = 0;
*xc = *yc = *xw1 = *xw2 = *yw1 = *yw2 = *max_level = 0;

/* Read the option first: */
 buffer[0] = '#';
 while(!feof(fp_batch) && buffer[0] == '#') fgets(buffer,80,fp_batch);
 if(sscanf(buffer, "option=%d", ioption) != 1){
    printf("Fatal error reading batch file: option not found!\n %s\n", buffer);
    exit(-1);
    }
 if(italk) printf("ioption=%d\n", *ioption);
 if(*ioption == 0) return(1);

/* Read coordinates of the center (RELATIVE TO THE CENTER OF THE FRAME!)
*/
 buffer[0] = '#';
 while(!feof(fp_batch) && buffer[0] == '#') fgets(buffer,80,fp_batch);
 if(sscanf(buffer, "center=%f,%f", xc, yc) != 2){
    printf("Fatal error reading batch file: center coord. not found!\n %s\n", buffer);
    exit(-1);
    }

/* Read other parameters 
*  Values of the major axes: X min and max, Y min and max, and brightness level:\n");
*/
 buffer[0] = '#';
 while(!feof(fp_batch) && buffer[0] == '#') fgets(buffer,80,fp_batch);
 if(sscanf(buffer, "param=%f,%f,%f,%f,%f", xw1, xw2, yw1, yw2, max_level) != 5){
    printf("Fatal error reading batch file: 5 param. not found!\n %s\n", buffer);
    exit(-1);
    }

 if(italk) printf("xc=%f yc=%f xw1=%f xw2=%f yw1=%f yw2=%f max_level=%f\n", 
                   *xc, *yc, *xw1, *xw2, *yw1, *yw2, *max_level);

/* Error handling: */
  if(*max_level == 0) status = -1;
  if(*ioption > 1 && *xw1 == 0 && *yw1 == 0) status = -1;

return(status);
}
/************************************************************************
* Add one object to the frame
*
************************************************************************/
static int add_one_object(float *tab, int nx, int ny, int idim,
                          int ioption, float xc, float yc, 
                          float xw1, float xw2, float yw1, float yw2,
                          float max_level)
{
int status = 0;

/* Main switch according to "ioption" : */
 switch (ioption) {
   case 1 : 
     pixel(tab, nx, ny, idim, xc, yc, max_level); 
     break;
   case 2 : 
     rectangle(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 3 : 
     pyramide(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 4 : 
     ellipse(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 5 : 
     ellip_ring(tab, nx, ny, idim, xc, yc, xw1, yw1, xw2, yw2, 
                max_level); 
     break;
   case 6 : 
     cone(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 7 : 
     auto_circ(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 8 : 
     truncated_cosinus(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 9 : 
     truncated_cosinus2(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 10 : 
     truncated_sinc2(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 11 : 
     truncated_airy(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 12 : 
     dimmed_ellipse(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 13 : 
     gauss(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 14 : 
     sinc2(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 15 : 
     airy(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   case 16 : 
     moffat(tab, nx, ny, idim, xc, yc, xw1, yw1, max_level); 
     break;
   default :
      fprintf(stderr,"Bad choice\n");
      status = -1;
      break;
   }
	
return(status);
}

/******************************************************************************
* pixel: setting the value of one pixel
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the pixel (relative to the center of the frame)
* level: brightness level
******************************************************************************/
static int pixel(float *tab, int nx, int ny, int idim, float xc, float yc, 
                 float level)
{
int status = 0;
int i, j;

   i = NINT(xc + nx/2);
   j = NINT(yc + ny/2);

   if(i >= 0 && i < nx && j >= 0 && j < ny) {
     tab[i + idim * j] = level;
   } else {
     printf("\n pixel/Warning: pixel out of window!\n"); 
     status = -1;
     }

return(status);
}

/******************************************************************************
* rectangle with uniform brightness
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y size of the rectangular support 
* level: brightness level
******************************************************************************/
static int rectangle(float *tab, int nx, int ny, int idim, float xc, float yc, 
                     float xwidth, float ywidth, float max_level)
{
register int i,j;
int xmin,xmax,ymin,ymax;

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

   for (j = ymin; j <= ymax; j++)
     for (i = xmin; i <= xmax; i++) tab[i + idim * j] += max_level;

return(0);
}
/******************************************************************************
* Square pyramide with trangle profile 
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y size of the side of the rectangular support
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int pyramide(float *tab, int nx, int ny, int idim, float xc, float yc, 
                    float xwidth, float ywidth, float max_level)
{
register int i,j;
int xmin,xmax,ymin,ymax;

#ifdef DEBUG
printf("nx=%d ny=%d idim=%d xc=%.1f yc=%.1f xwidth=%.1f ywidth=%.1f, max_level=%.2f\n",
       nx, ny, idim, xc, yc, xwidth, ywidth, max_level);
#endif
  
   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

#ifdef DEBUG
printf("xmin=%d xmax=%d ymin=%d ymax=%d\n", xmin, xmax, ymin, ymax);
#endif

   for (j = ymin; j <= ymax; j++)
     for (i = xmin; i <= xmax; i++) 
         tab[i + idim * j] += max_level 
                              * ( 1. - 2. * fabs(i - xc) / xwidth )
                              * ( 1. - 2. * fabs(j - yc) / ywidth );

return(0);
}
/******************************************************************************
* ellipse: ellipse with constant brightness
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int ellipse(float *tab, int nx, int ny, int idim, float xc, float yc, 
                   float xwidth, float ywidth, float max_level)
{
register int i,j;
int xmin, xmax, ymin, ymax;
float x2, y2, xw2, yw2, r2;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

   for (j = ymin; j < ymax; j++) 
   {
      y2 = SQUARE(j-yc);
      for (i = xmin; i < xmax; i++)
      {
         x2 = SQUARE(i-xc);
         r2 = x2 / xw2 + y2 / yw2;
         if (r2 <= 1.) tab[i + idim * j] += max_level;
      }
      
   }

return(0);
}

/******************************************************************************
* ellip_ring: elliptical ring with constant brightness
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth1, ywidth1: x and y major axes of the internal ellipse
*                  (size of the rectangular support)
* xwidth2, ywidth2: x and y major axes of the external ellipse
*                  (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int ellip_ring(float *tab, int nx, int ny, int idim, float xc, float yc, 
                      float xwidth1, float ywidth1, float xwidth2, 
                      float ywidth2, float max_level)
{
register int i,j;
int xmin, xmax, ymin, ymax;
float x2, y2, xw12, yw12, xw22, yw22, r12, r22;

   xw12 = SQUARE(xwidth1 / 2.);
   yw12 = SQUARE(ywidth1 / 2.);
   xw22 = SQUARE(xwidth2 / 2.);
   yw22 = SQUARE(ywidth2 / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth2 / 2.);
   xmax = NINT(xc + xwidth2 / 2.);
   ymin = NINT(yc - ywidth2 / 2.);
   ymax = NINT(yc + ywidth2 / 2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

   for (j = ymin; j < ymax; j++) 
   {
      y2 = SQUARE(j-yc);
      for (i = xmin; i < xmax; i++)
      {
         x2 = SQUARE(i-xc);     
         r22 = (float) (x2 / xw22 + y2 / yw22);
         r12 = (float) (x2 / xw12 + y2 / yw12);
         if ( r22 <= 1. && r12 >= 1.) tab[i + idim * j] += max_level;
      }
   }

return(0);
}
/******************************************************************************
* cone: disk with linearly decressing brightness
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int cone(float *tab, int nx, int ny, int idim, float xc, float yc, 
                float xwidth, float ywidth, float max_level)
{
register int i,j;
int xmin, xmax, ymin, ymax;
float x2, y2, xw2, yw2, r2;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);


   for (j = ymin; j < ymax; j++) 
   {
      y2 = SQUARE(j-yc);
      for (i = xmin; i < xmax; i++)
      {
        x2 = SQUARE(i-xc);
        r2 = x2 / xw2 + y2 / yw2;
        if (r2 <= 1.) 
           tab[i + idim * j] += max_level * ( 1. - 2. * fabs(i - xc) / xwidth )
                                * ( 1. - 2. * fabs(j - yc) / ywidth );
      }
   }

return(0);
}
/******************************************************************************
* auto_circ: compute the result of the auto-correlation of a uniform disk
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int auto_circ(float *tab, int nx, int ny, int idim, float xc, float yc, 
                     float xwidth, float ywidth, float max_level)
{
register int i, j;
int xmin,xmax,ymin,ymax;
float x2, y2, xw2, yw2, r2, rad;

if(xwidth != ywidth) {
  printf("auto_circ/Sorry only a disk is allowed for this option!\n");
  return(-1);
  }

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

  for (j = ymin; j < ymax; j++) 
  {
      y2 = SQUARE(j-yc);
      for (i = xmin; i < xmax; i++)
      {
        x2 = SQUARE(i-xc);     
        r2 = x2 / xw2 + y2 / yw2;
        rad = sqrt(r2);
        if (r2 <= 1.) tab[i + idim * j] += max_level * (2./ PI) *
                                   (acos(rad) - rad * sqrt( 1. - SQUARE(rad)));
     }
  }

return(0);
}

/******************************************************************************
* Ellipse with truncated cosinus function decreasing down to the edges
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int truncated_cosinus(float *tab, int nx, int ny, int idim, float xc, 
                             float yc, float xwidth, float ywidth, 
                             float max_level)
{
register int i, j;
int xmin,xmax,ymin,ymax;
float x2, y2, xw2, yw2, r2, ss;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

  for (j = ymin; j < ymax; j++) 
  {
      y2 = SQUARE(j-yc);
      for (i = xmin; i < xmax; i++)
      {
        x2 = SQUARE(i-xc);     
        r2 = x2 / xw2 + y2 / yw2;
        ss = sqrt(r2) * PI / 2.;
        if (r2 <= 1.) tab[i + idim * j] += max_level * cos(ss); 
        }
   }

return(0);
}
/******************************************************************************
* Ellipse with truncated squared cosinus function decreasing down to the edges
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int truncated_cosinus2(float *tab, int nx, int ny, int idim, float xc, 
                              float yc, float xwidth, float ywidth, 
                              float max_level)
{
register int i, j;
int xmin,xmax,ymin,ymax;
float x2, y2, xw2, yw2, r2, ss;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

  for (j = ymin; j < ymax; j++) 
  {
      y2 = SQUARE(j-yc);
      for (i = xmin; i < xmax; i++)
      {
        x2 = SQUARE(i-xc);     
        r2 = x2 / xw2 + y2 / yw2;
        ss = sqrt(r2) * PI / 2.;
        if (r2 <= 1.) tab[i + idim * j] += max_level * cos(ss) * cos(ss); 
        }
   }

return(0);
}
/******************************************************************************
* Ellipse with truncated sinc2 function decreasing down to the edges
* [sin (PI x) / x ]^2
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int truncated_sinc2(float *tab, int nx, int ny, int idim, float xc, 
                           float yc, float xwidth, float ywidth, 
                           float max_level)
{
register int i, j;
int xmin,xmax,ymin,ymax;
float x2, y2, xw2, yw2, r2, ss;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

  for (j = ymin; j < ymax; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = xmin; i < xmax; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      ss = sqrt(r2) * PI ;
      if(r2 < 0.0001) tab[i + idim * j] += max_level;
      else if (r2 <= 1.) tab[i + idim * j] += max_level * SQUARE(sin(ss) / ss); 
    }
  }

return(0);
}
/******************************************************************************
* Ellipse with truncated Airy function decreasing down to the edges
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int truncated_airy(float *tab, int nx, int ny, int idim, float xc, 
                          float yc, float xwidth, float ywidth, 
                          float max_level)
{
register int i, j;
int xmin,xmax,ymin,ymax;
float x2, y2, xw2, yw2, r2, ss;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

  for (j = ymin; j < ymax; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = xmin; i < xmax; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      ss = sqrt(r2) * 3.83171 ;
      if(r2 < 0.001) {
        tab[i + idim * j] += max_level;
      } else if (r2 <= 1.) { 
         tab[i + idim * j] += max_level * 4. * SQUARE(j1(ss) / ss); 
      } /* EOF if */
    }  /* EOF loop on i */
  }  /* EOF loop on j */

return(0);
}
/******************************************************************************
* Dimmed ellipse decreasing down to the edges
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int dimmed_ellipse(float *tab, int nx, int ny, int idim, float xc, 
                          float yc, float xwidth, float ywidth, 
                          float max_level)
{
register int i, j;
int xmin,xmax,ymin,ymax;
float x2, y2, xw2, yw2, r2, dzeta, work, a1, a2;

   printf("Parameters of the dimmed disk:\n");
   printf(" Here we will use a polynomial of dzeta with:\n");
   printf(" dzeta = ( 1. - sqrt(1. - radius**2 / max_radius**2))  \n");
   printf(" output = bmax * (1. - a1 * dzeta - a2 * dzeta**2) \n\n"); 
   printf("Dimming polynomial coefficients (first, second order): a1,a2 \n");
   scanf("%f,%f",&a1,&a2);

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;
   xmin = NINT(xc - xwidth/2.);
   xmax = NINT(xc + xwidth/2.);
   ymin = NINT(yc - ywidth/2.);
   ymax = NINT(yc + ywidth/2.);
   xmin = MAXI( MINI(nx, xmin), 0);
   xmax = MAXI( MINI(nx, xmax), 0);
   ymin = MAXI( MINI(ny, ymin), 0);
   ymax = MAXI( MINI(ny, ymax), 0);

  for (j = ymin; j < ymax; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = xmin; i < xmax; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      if (r2 <= 1.) { 
         dzeta = 1. - sqrt(1. - r2); 
         work = max_level * (1. - a1 * dzeta - a2 * dzeta * dzeta);
         if(work > 0.) tab[i + idim * j] += work;
      } 
    }  /* EOF loop on i */
  }  /* EOF loop on j */

return(0);
}
/******************************************************************************
* Gauss function
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int gauss(float *tab, int nx, int ny, int idim, float xc, float yc, 
                 float xwidth, float ywidth, float max_level)
{
register int i, j;
float x2, y2, xw2, yw2, r2, ss;

   xw2 = SQUARE(xwidth / 2.);
   yw2 = SQUARE(ywidth / 2.);

   xc += nx/2;
   yc += ny/2;

  for (j = 0; j < ny; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = 0; i < nx; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      ss = - r2 * log(2.);
      tab[i + idim * j] += max_level * exp (ss); 
    }
  }
		
return(0);
}		
/******************************************************************************
* sinc2 function
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int sinc2(float *tab, int nx, int ny, int idim, float xc, float yc, 
                 float xwidth, float ywidth, float max_level)
{
register int i, j;
float x2, y2, xw2, yw2, r2, ss;

  xw2 = SQUARE(xwidth / 2.);
  yw2 = SQUARE(ywidth / 2.);

  xc += nx/2;
  yc += ny/2;

  for (j = 0; j < ny; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = 0; i < nx; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      ss = sqrt(r2) * PI ;
      if(r2 < 0.0001) tab[i + idim * j] += max_level;
      else tab[i + idim * j] += max_level * SQUARE(sin(ss) / ss); 
    }
  }

return(0);
}
/******************************************************************************
* Airy function
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int airy(float *tab, int nx, int ny, int idim, float xc, float yc, 
                float xwidth, float ywidth, float max_level)
{
register int i, j;
float x2, y2, xw2, yw2, r2, ss;

  xw2 = SQUARE(xwidth / 2.);
  yw2 = SQUARE(ywidth / 2.);

  xc += nx/2;
  yc += ny/2;

  for (j = 0; j < ny; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = 0; i < nx; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      ss = sqrt(r2) * 3.83171 ;
      if(r2 < 0.001) 
        tab[i + idim * j] += max_level;
      else  
        tab[i + idim * j] += max_level * 4. * SQUARE(j1(ss) / ss); 
    }  /* EOF loop on i */
  }  /* EOF loop on j */

return(0);
}
/******************************************************************************
* Moffat profile 
*
* INPUT:
* tab[]: array to be updated
* nx, ny: size of the frame
* idim: size of the first dimension of tab
* xc, yc: coordinates of the center (relative to the center of the frame)
* xwidth, ywidth: x and y major axes (size of the rectangular support)
* max_level: brightness level of the center (i.e. maximum)
******************************************************************************/
static int moffat(float *tab, int nx, int ny, int idim, float xc, float yc, 
                  float xwidth, float ywidth, float max_level)
{
register int i, j;
float x2, y2, xw2, yw2, r2, ss, bb, qq;
bb = 1.; 
qq = 1.5;
printf("Moffat's profile = max_level /(1 + (r b)^2)^q  (with b=%.1f q=%.1f)\n", 
        bb, qq);

  xw2 = SQUARE(xwidth / 2.);
  yw2 = SQUARE(ywidth / 2.);

  xc += nx/2;
  yc += ny/2;

  for (j = 0; j < ny; j++) 
  {
    y2 = SQUARE(j-yc);
    for (i = 0; i < nx; i++)
    {
      x2 = SQUARE(i-xc);     
      r2 = x2 / xw2 + y2 / yw2;
      ss = sqrt(r2) * bb ;
      tab[i + idim * j] += max_level / pow( 1. + SQUARE(ss), qq); 
    }  /* EOF loop on i */
  }  /* EOF loop on j */

return(0);
}
