#include "mex.h"

bool bDebug=false;

/*-----------------------------------------------------------------------
function HC=hist_cost_2(sc1,sc2)
	
	Suppose each COLUMN of sc1,sc2 is a chape contex feature at a given 
point.

	[nsamp1,nbins]=size(sc1);
	[nsamp2,nbins]=size(sc2);
	
	sc1n	= sc1./ repmat(sum(sc1,2)+eps,[1 nbins]);
	sc2n	= sc2./ repmat(sum(sc2,2)+eps,[1 nbins]);
	tmp1	= repmat(permute(sc1n,[1 3 2]),[1 nsamp2 1]);
	tmp2	= repmat(permute(sc2n',[3 2 1]),[nsamp1 1 1]);
	HC		= 0.5*sum(((tmp1-tmp2).^2)./(tmp1+tmp2+eps),3);
	
return;

  
 ------------------------------------------------------------------------*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	mxArray	*pMxHC;						/* input and output */
	double	*pSC1,   *pSC2,   *pHC;

	double	*sc1,*sc2,*pCur;
	int		nSamp1,nSamp2,nBin1,nBin2;
	double  coff;  // for weight
	int		r,c,k;
	double	dis,tmp;


    /* Analyse input data */
    
	pSC1	= mxGetPr(prhs[0]);
    nSamp1	= mxGetN(prhs[0]);
    nBin1	= mxGetM(prhs[0]);
    
    pSC2	= mxGetPr(prhs[1]);
    nSamp2	= mxGetN(prhs[1]);
    nBin2	= mxGetM(prhs[1]);


    /* Create output matrices and initilize */
    pMxHC	= mxCreateDoubleMatrix(nSamp1,nSamp2,mxREAL);
    pHC		= mxGetPr(pMxHC);



	/* compute the distance between the r-th point in shape 1
								and the c-th point in shape 2 */
	pCur	= pHC;
	for(c=0;c<nSamp2;c++)
	{
		sc2	= pSC2+c*nBin2;
		for(r=0;r<nSamp1;r++)
		{
			sc1	= pSC1+r*nBin1;
			dis	= 0;
			for(k=0;k<nBin1;k++)
			{
				coff = (k+1 < nBin1-k) ? k+1 : nBin1-k;

				tmp	= sc1[k]-sc2[k];
				if(tmp>=0)  // coff = [1:12 12:1]
					dis += tmp/(coff);
				else
					dis -= tmp/(coff);
			}

			*pCur	= dis;
			pCur++;
		}
	}


    /* Return */
    plhs[0] = pMxHC;
}

// definition of weight:
// the weight is symmetric, i.e., w1=w24, w2=w23, w3=w22...,w12=w13
// thus w(k) = 1 / (min{k+1, nBin1+1-k-1} + 1)