/******************************************************************************/
/**                                                                          **/
/**                        P H E N O .  C                                    **/
/**                                                                          **/
/**   pheno.c: implementation file for R pheno package C functions	     **/
/**                                                                          **/
/**   Last changes: 2/2/2006                                                 **/
/**                                                                          **/
/**   written by Joerg Schaber                                               **/
/**				Max Planck Institute for Molecular Genetics  **/
/**                                                                          **/
/**   This code is subject to the GNU PUBLIC LICENSE 2 or higher	     **/
/******************************************************************************/
/* ANSI C implementation of some auxiliary functions
 *
 * 1. date2jul1(int *doy, int *year, char *date[])
 *     wandet ein Datumsstring der Form DD.MM.YYYY in eine
 *     julianischen Datum um.
 *
 * 2. date2jul2(int *doy, int *year, int *month, int *day)
 *      wandelt ein Datum in DOY um, Schaltjahre werden mit
 *      beruecksichtigt
 *
 * 3. daysbetween(char *date1[], char *date2[], int *ndays)
 *    liefert die Anzahl Tage ziwschen datew 1und date2 als int,
 *    wobei Datum im Format DD.MM.YYYY angegeben sein muss.
 *
 * 4. daylength(double *lat, int *iday, double *dl, double *delta)
 *    calculates daylength dl [h] and declination delta  [radians]
 *    on day iday for latitude lat [degrees]
 *    declination: angle between sun rays and equatorial plane for 
 *    the whole earth (-23 degrees - + 23 degrees)
 *
 * 5. jul2date1(int *doy, int *year, char *date[]);
 *     wandelt ein julianischen Datum (Tag im Jahr, Jahr)
 *     in ein Datum als String der Form DD.MM.YYYY um.
 *
 * 6. jul2date2(int *doy, int *year, int *day, int *month)
 *     wandelt ein julianischen Datum (Tag im Jahr, Jahr)
 *     in ein Datum als int monat, Tag um
 *
 * 7. leapyear(int *year, int *ly)
 *     ly = 1 = leapyear
 *     ly = 0 = kein leapyear
 *
 * 8. maxdaylength(double *lat, double *maxdl)
 *    calculates maximal daylength dl [h] at a certain latitude lat [degrees]
 * 
 * 9. connectivity: finds connected data sets of a matrix

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pheno.h"

/* auxiliary function for integer conversion */
void myitoa(int n, char *s)
{
    int len=0,sign,i=0,j,c;

    if((sign = n) < 0) n = -n; /* Vorzeichen merken */

    /* Ziffern von recht her generieren */
    do {

        *(s+i) = n % 10 + '0';  /* naechste Ziffer */
        len++;                  /* Laenge mitzaehlen*/
        i++;

    } while ((n /= 10 ) > 0);

    if (sign < 0)
    {
        *(s+i) = '-';
        len++;
    }

    /*  jetzt umdrehen  */
    for(i=0,j=len-1; i < j; i++, j--)
    {
        c = *(s+i);
        *(s+i) = *(s+j);
        *(s+j) = c;
    }

    /*  abschliessen  */
    *(s+len) = '\0';
}

void date2jul1(char *date[], int *doy, int *year)
{
    int day, month, ly;
    int month_end[13]={0,31,59,90,120,151,181,212,243,273,304,334,365};

    day   = atoi(*date);
    month = atoi((*date)+3);
    *year = atoi((*date)+6);

	leapyear(year,&ly);
	
    if(ly && month > 2)
    {
        *doy = month_end[month-1]+1+day;
    }
    else *doy = month_end[month-1]+day;
}
 
void date2jul2(int *year, int *month, int *day, int *doy)
{
    int ly,month_end[13]={0,31,59,90,120,151,181,212,243,273,304,334,365};

	leapyear(year,&ly);

    if(ly && *month > 2)
    {
        *doy = month_end[(*month)-1]+1+(*day);
    }
    else *doy = month_end[(*month)-1]+(*day);

}

void daysbetween(char *date1[], char *date2[], int *ndays)
{
    int ly,doy1,doy2,year1,year2,i,tmp;

    date2jul1(date1,&doy1,&year1);
    date2jul1(date2,&doy2,&year2);

    if(year1==year2) *ndays = abs(doy1-doy2);
    else if( year1 > year2 )
    { 
		leapyear(&year2,&ly);

        if(ly) tmp = 366 - doy2;
        else tmp = 365 - doy2;

        for(i=year2+1;i<year1;i++)
        {
			leapyear(&i,&ly);

            if(ly) tmp = tmp + 366;
            else tmp = tmp + 365;
        }

        *ndays = (tmp + doy1);
    }
    else
    {
		leapyear(&year1,&ly);

        if(ly) tmp = 366 - doy1;
        else tmp = 365 - doy1;

        for(i=year1+1;i<year2;i++)
        {
			leapyear(&i,&ly);

            if(ly) tmp = tmp + 366;
            else tmp = tmp + 365;
        }

        *ndays = (tmp + doy2);
    }
}

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

void jul2date1(int *doy, int *year, char *date[])
{
    char day[3],monthc[3],yearc[5];
    char *day_ptr, *monthc_ptr, *yearc_ptr;
    int ly,month = 1;
    int month_end[13]={0,31,59,90,120,151,181,212,243,273,304,334,365};

	leapyear(year,&ly);

    if( ly )
    {
        while( month_end[month] < *doy )
        {
           month++;
           month_end[month]++;
        }
    }
    else
    {
        while( month_end[month] < *doy )
        {
           month++;
        }
    }

    day_ptr    = &day[0];
    monthc_ptr = &monthc[0];
    yearc_ptr  = &yearc[0];

    myitoa(*doy - month_end[month-1],day_ptr);
    myitoa(month,monthc_ptr);
    myitoa(*year,yearc_ptr);

    if ( *doy - month_end[month-1] < 10 )
    {
        date[0][0] = '0';
        date[0][1] = *day_ptr;
    }
    else
    { 
        date[0][0] = *day_ptr;
        date[0][1] = *(day_ptr+1);
    }

    date[0][2]  = '.';

    if ( month < 10 )
    {
        date[0][3] = '0';
        date[0][4] = *monthc_ptr;
    }
    else
    {
        date[0][3] = *monthc_ptr;
        date[0][4] = *(monthc_ptr+1);
    }  

    date[0][5] =  '.';
    date[0][6] = *yearc_ptr;
    date[0][7] = *(yearc_ptr+1);
    date[0][8] = *(yearc_ptr+2);
    date[0][9] = *(yearc_ptr+3);
    date[0][10] = '\0';
}

void jul2date2(int *doy, int *year, int *day, int *month)
{
    int ly, monthc = 1;
    int month_end[13]={0,31,59,90,120,151,181,212,243,273,304,334,365};

	leapyear(year,&ly);
	
    if( ly )
    {
        while( month_end[monthc] < *doy )
        {
           monthc++;
           month_end[monthc]++;
        }
    }
    else
    {
        while( month_end[monthc] < *doy )
        {
           monthc++;
        }
    }

    *day = *doy - month_end[monthc-1];
    *month = monthc;
}

void leapyear(int *year, int *ly)
{
    if( (*year)%400==0 || ( (*year)%100!=0 && (*year)%4==0 ))
    {
        *ly = 1;
    }
    else
    {
        *ly = 0;
    }
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
