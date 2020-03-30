/**************************************************************************
* Prototypes of routines defined in kolmo.c
*
* JLP
* Version 09/06/2006
**************************************************************************/
#ifndef _kolmo_h /* BOF sentry */
#define _kolmo_h
int get_gaussft(float *array, int nx, int ny, int idim);
int get_phift(float *array, int nx, int ny, int idim, float L0, float r0);
int get_phasor(float *array, int nx, int ny, int idim);
int recent_real(float *array, int nx, int ny, int idim);
int recent_complex(float *array, int nx, int ny, int idim);
int mask_pupil(float *array, int nx, int ny, int idim, float diam,
               float cent_obscu);
int output_pupil(float *array, int nx, int ny, int idim, char *re_name,
                 char *im_name, char *re_comments, char *im_comments);
int  get_intensity(float *array, float *intens, int nx, int ny, int idim,
                   int nphot);
int  autocor_cplx(float *array, int nx, int ny, int idim, float *ree,
                  float *ima);
int  autocor_real(float *array, int nx, int ny, int idim, float *ree);
int  structure_real(float *array, int nx, int ny, int idim, float *ree);
int  structure_cplx(float *array, int nx, int ny, int idim, float *ree);

#endif /* EOF sentry */
