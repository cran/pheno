/******************************************************************************/
/**                                                                          **/
/**                        P H E N O .  H                                    **/
/**                                                                          **/
/**   pheno.h: header file for libpheno.a, a library for diverse functions   **/
/**   implementation : pheno.c                                               **/
/**                                                                          **/
/**   Last changes: 31/10/2003                                               **/
/**                                                                          **/
/**   written by Joerg Schaber                                               **/
/**              Potsdam Institute for Climate Impact Research               **/
/**              P.O. Box 60 12 03                                           **/
/**              D-14412 Potsdam/Germany                                     **/
/**                                                                          **/
/**   This code is subject to the GNU PUBLIC LICENSE 2 or higher	         **/
/******************************************************************************/

void daylength(double *lat, int *iday, double *dl, double *delta);

void maxdaylength(double *lat, double *maxdl);

void connectivity(double tmp[],int *maxnr,int *maxnc,int row_classes[],int col_clases[]);
