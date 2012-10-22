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

void Cdate2jul1(char *date[], int *doy, int *year);

void Cdate2jul2(int *year, int *month, int *day, int *doy);

void Cdaysbetween(char *date1[], char *date2[], int *ndays);

void Cdaylength(double *lat, int *iday, double *dl, double *delta);

void Cjul2date1(int *doy, int *year, char *date[]);

void Cjul2date2(int *doy, int *year, int *day, int *month);

void Cleapyear(int *year, int *ly);

void Cmaxdaylength(double *lat, double *maxdl);

void Cconnectivity(double tmp[],int *maxnr,int *maxnc,int row_classes[],int col_clases[]);
