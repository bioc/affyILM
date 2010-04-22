#include "ilm.h"

extern int nAG,nCT;
extern int matl;
extern int ix[NBS],iy[NBS];
extern double Xty[MDIM];
extern double logIntens[MAX_CHIPDIM][MAX_CHIPDIM];
extern int freq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
extern int parfreq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];


/*Vector generation Xt*logIntens: */
void multXty(int *xvalueAG, int *yvalueAG, int *xvalueCT, int *yvalueCT)
{
	int x,y;
	int i,rAG,rCT;
	double logI;

	for(rAG=0;rAG<nAG;rAG++){
		x=xvalueAG[rAG];
		y=yvalueAG[rAG];
		logI=logIntens[x][y];
		Xty[0]+=logI;
		for (i=1;i<9;i++){
			Xty[i]+=logIntens[x+ix[i]][y+iy[i]]*logI;
		}
	}

	for (rCT=0;rCT<nCT;rCT++) {
		x=xvalueCT[rCT];
		y=yvalueCT[rCT];
		logI=logIntens[x][y];
		Xty[9]+=logI;
		for (i=10;i<18;i++) {
			Xty[i]+=logIntens[x+ix[i-9]][y+iy[i-9]]*logI;
		}
	}

	for (rAG=0;rAG<nAG;rAG++) {
		x=xvalueAG[rAG];
		y=yvalueAG[rAG];
		logI=logIntens[x][y];
		for(i=18;i<34;i++) {
			Xty[i] += freq[x][y][i-18]*logI;
		}
		for(i=34;i<matl;i++) {
			Xty[i] += parfreq[x][y][i-34]*0.1*logI;
		}
	}

	for (rCT=0;rCT<nCT;rCT++) {
		x=xvalueCT[rCT];
		y=yvalueCT[rCT];
		logI=logIntens[x][y];
		for(i=18;i<34;i++) {
			Xty[i] += freq[x][y][i-18]*logI;
		}
		for(i=34;i<matl;i++) {
			Xty[i] += parfreq[x][y][i-34]*0.1*logI;
		}
	}
}



