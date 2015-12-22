/*	This file contains some useful small functions frequently used in 
	C code for matlab.
	Haibin Ling, hbling AT umiacs.umd.edu
	07/31/2004
*/

#ifndef COMMON_MATLAB
#define COMMON_MATLAB

#include "mex.h"
#include "matrix.h"

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"


/* common functions */

#define MAT_GET(pMat,x,y,nRow)		pMat[(x)*(nRow)+(y)]
#define MAT_SET(pMat,x,y,val,nRow)	pMat[(x)*(nRow)+(y)]=(val)


#define	XY2ID(x,y,nRow)		((x)*(nRow)+(y))
#define	IDX2X(idx,nRow)		((int)((idx)/(nRow)))
#define	IDX2Y(idx,nRow)		((idx)%(nRow))


#define SAFE_FREE(p)		if(p){free(p);p=0;}

/* some debugging macros */
#define WRONG1(b)			if(b){printf("wrong\t");}
#define WRONG2(b,s)			if(b){printf(s);}
#define WRONG3(b,s,v)		if(b){printf(s); return (v);}


/* simple math marcros */
#define MIN(a,b)	((a)<(b)?(a):(b))
#define MAX(a,b)	((a)>(b)?(a):(b))
#define ROUND(x)	((int)((x)+.5))


/* #define PI		3,14159265	*/



/*------------------------------------------------------------------------------------
 compute (high dimension) distance matrix between two point sets which are represented 
 by matrix.
 Matlab code as:	D	= sqrt(dist2(X1,X2));
*/
bool myDist(double* D, double* X1,int m1, double* X2, int m2, int dim)
{
	int i,j,k;
	double	v,t;
	if(!D||!X1||!X2)	return false;

	for(i=0;i<m1;i++){
		for(j=0;j<m2;j++){
			v = 0;
			for(k=0;k<dim;k++)	{	
				t	= MAT_GET(X1,k,i,m1)-MAT_GET(X2,k,j,m2);
				v  += t*t;
			}
			MAT_SET(D,j,i,sqrt(v),m1);
		}
	}
	return true;
}




#endif
