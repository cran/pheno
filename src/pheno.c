/******************************************************************************/
/**                                                                          **/
/**                        P H E N O .  C                                    **/
/**                                                                          **/
/**   pheno.c: implementation file for libpheno.a, a library 				 **/
/** 		for diverse functions 											 **/
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
/* ANSI C implementation of some auxiliary functions
 * for the calculation of combined phenological time series
 *
 * 1. daylength(double lat, int iday, double *dl, double *delta)
 *    calculates daylength dl [h] and declination delta  [radians]
 *    on day iday for latitude lat [degrees]
 *    declination: angle between sun rays and equatorial plane for 
 *    the whole earth (-23 degrees - + 23 degrees)
 * 2. maxdaylength(double lat)
 *    calculates maximal daylength dl [h] at a certain latitude lat [degrees]
 * 3. connectivity: finds connected data sets of a matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pheno.h"

void daylength(double *lat, int *iday, double *dl, double *delta)
{
    double arg;

    *delta = -23.45*(M_PI/180)*cos(2.*M_PI/365.*((double) *iday+10));

    /* latitude is converted to rad */
    arg = -tan(*lat*M_PI/180.)*tan(*delta);

    if( arg < -1 ) *dl = 24;
    else if ( arg > 1 ) *dl = 0;
    else *dl = (24./M_PI)*acos(arg);
}

void maxdaylength(double *lat, double *maxdl)
{
    double arg,delta;

    delta = -23.45*(M_PI/180)*cos(2.*M_PI/365.*((double) 172+10));

    /* latitude is converted to rad */
    arg = -tan(*lat*M_PI/180.)*tan(delta);

    if( arg < -1 ) *maxdl = 24;
    else if ( arg > 1 ) *maxdl = 0;
    else *maxdl = (24./M_PI)*acos(arg);
}

/* auxiliary function for Connectivity
 * Searches data entries in lines and colums 
 * of the data matrix tmp[0..maxnr-1][0..maxnc-1]
 */
void con_step(double tmp[],int maxnc,int classes[],int index,int no_classes,int row_or_col_step,int check[],int *from,int *len,int class_nr) 
{
	int i;
	
	/* already checked, index: row/column to be checked */
	if(classes[index] != -1) return;
	else {
		/* search over columns for entries */
		if(row_or_col_step == 0) {
/*			puts("over cols");
*/			
			for(i=(*from);i<no_classes;i++) { 	/* maxnc */
/*				printf("tmp[%d * %d + %d = %d] = %.1f\n",maxnc,index,i,maxnc*index+i,tmp[maxnc*index+i]);
*/
				if(tmp[maxnc*index+i] != 0) {
					(*len)++;					/* col_len */
					check[(*len)-1] = i;  		/* check_col */
					if(i == (*from)) (*from)++;
				}
			}
			classes[index] = class_nr;			/* row_classes */
		}
		/* search over rows for entries */
		else {
/*			puts("over rows");
*/			
			for(i=(*from);i<no_classes;i++) {	/* maxnr */
/*				printf("tmp[%d * %d + %d = %d] = %.1f\n",maxnc,i,index,maxnc*i+index,tmp[maxnc*index+i]);
*/
				if(tmp[maxnc*i+index] != 0) {
					(*len)++;					/* row_len */
					check[(*len)-1] = i;		/* check_row */
					if(i == (*from)) (*from)++;
				}
			}
			classes[index] = class_nr;			/* col_classes */
		}
	}
}

/* Auxiliary function for Connectivity
 * Investigates number of observations per class,
 * and finds all rows/columns connected to the current row
 */
void exhaust_class(double tmp[],int row_classes[],int col_classes[],int check_row[],int check_col[],int *from_row,int *from_col,int no_rows,int no_cols,int class_nr)
{
	int i;
	int row_step = 0;
	int col_step = 1;

	int finished;
	int row_len=1;
	int col_len=0;

	finished = row_len;

	while(finished != 0) {

/*		puts("");
		printf("check_row :[");
		for(i=0;i<row_len;i++) printf("%d,",check_row[i]);
		puts("]");
		printf("check_col :[");
		for(i=0;i<col_len;i++) printf("%d,",check_col[i]);
		puts("]");
*/
		/* do first a row step, searching over columns */
		if(row_len != 0) {
			con_step(tmp,no_cols,row_classes,check_row[0],no_cols,row_step,check_col,from_col,&col_len,class_nr);
			/* eliminate the row checked */
			for(i=0;i<row_len;i++) check_row[i]=check_row[i+1];
			row_len--;
		}

		/* then do a column step, searching over rows */
		if(col_len != 0) {
			con_step(tmp,no_cols,col_classes,check_col[0],no_rows,col_step,check_row,from_row,&row_len,class_nr);
			/* eliminate the colmun checked */
			for(i=0;i<col_len;i++) check_col[i]=check_col[i+1];
			col_len--;
		}

		finished = row_len + col_len;
/*		puts("Rowclas");
		for(i=0;i<no_rows;i++) printf("row %d: %d\n",i,row_classes[i]);
		puts("Colclas");
		for(i=0;i<no_cols;i++) printf("col %d: %d\n",i,col_classes[i]);
*/
	}
}


/* auxiliary function for Connectivity
 * new connected set, no_rows number of connected rows 
 */
int pick_row(int row_classes[], int no_rows)
{
	int i;
	int picked = -1; /* not yes classified */
	
	for(i=0;i<no_rows;i++) {
		if(row_classes[i] == -1) {
			picked = i;
			break;
		}
	}
	return picked;
}
	
/* Finds connected data sets of a matrix tmp[0..maxnr-1][0..maxnc-1].
 * tmp is saved rowwise.
 * Non-entries in the matrix are considered 0
 * Returns two vectors: 
 * row_classes[0..maxnr-1] : Class number of the respective rows
 * col_classes[0..maxnc-1] : Class number of the respective cols
 */
void connectivity(double tmp[],int *maxnr,int *maxnc,int rowclasses[],int colclasses[])
{
	int *check_row;
	int *check_col;

	int class_nr = 0;			/* number of connected sets */
	int from_row = 1;			/* starting row to search, first row is already picked */
	int from_col = 0;			/* starting column to search */
	int picked,i;

	check_row=(int *) calloc((*maxnr)*(*maxnc),sizeof(int)); /* checked rows for each set */
	check_col=(int *) calloc((*maxnr)*(*maxnc),sizeof(int)); /* checked colums for each set */
	
/*	puts("TMP");
	for(i=0;i<(*maxnr)*(*maxnc);i++) printf("%.0f ",tmp[i]);
	puts("");
*/
	for(i=0;i<*maxnr;i++) rowclasses[i]=-1;
	for(i=0;i<*maxnc;i++) colclasses[i]=-1;

	/* new connected block */
	picked = pick_row(rowclasses,*maxnr);

	/* as long as there are rows not yet classified */
	while(picked != -1) {
		class_nr++;
		check_row[0] = picked;
		/* find all rows/column connected to the current row */
		exhaust_class(tmp,rowclasses,colclasses,check_row,check_col,&from_row,&from_col,*maxnr,*maxnc,class_nr);
		picked = pick_row(rowclasses,*maxnr);
	}

/*	printf("connectedSets found %d connected sets\n",class_nr);
*/
	free(check_row);
	free(check_col);
}
