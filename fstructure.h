/**************************************************************************
* Prototypes of routines defined in structure.c
*
* JLP
* Version 21/07/2008
**************************************************************************/
#ifndef _structure_h /* BOF sentry */
#define _structure_h
#define FSTRUCTURE_REAL         RENAME_(fstructure_real)
#define FSTRUCTURE_REAL_WITH_FT RENAME_(fstructure_real_with_ft)
#define FSTRUCTURE_CPLX         RENAME_(fstructure_cplx)
int  FSTRUCTURE_REAL(float *array, int *nx, int *ny, int *idim, float *ree, 
                    int *fourn_format);
int  FSTRUCTURE_REAL_WITH_FT(float *phase_array, int *nx, int *ny,
                            int *idim, float *fstruct, int *fourn_format);
int  FSTRUCTURE_CPLX(float *array, int *nx, int *ny, int *idim, float *ree);
#endif /* EOF sentry */
