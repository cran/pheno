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

void date2jul1(char *date[], int *doy, int *year);

void date2jul2(int *year, int *month, int *day, int *doy);

void daysbetween(char *date1[], char *date2[], int *ndays);

void daylength(double *lat, int *iday, double *dl, double *delta);

void jul2date1(int *doy, int *year, char *date[]);

void jul2date2(int *doy, int *year, int *day, int *month);

void leapyear(int *year, int *ly);

void maxdaylength(double *lat, double *maxdl);

void connectivity(double tmp[],int *maxnr,int *maxnc,int row_classes[],int col_clases[]);
