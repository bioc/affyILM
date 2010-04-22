#include "ilm.h"

extern int maxl;
extern int nAG,nCT;
extern double thval;

extern int WeirdY[MAX_CHIPDIM][MAX_CHIPDIM];
extern double logIntens[MAX_CHIPDIM][MAX_CHIPDIM];
extern char seq[MAX_CHIPDIM][MAX_CHIPDIM][26];

/*Counts the number of sequences where middle nucleotide of MM sequence is purine (nAG) or pyrimidine (nCT) and where absolute intensity < threshold */
int points(double thresval) 
{
	int wieviele;
	int x,y;

	wieviele=0;
	nAG=nCT=0;
	//loops only over MM's
	for (x=4;x<maxl-4;x++) for (y=4;y<maxl-4;y+=2) {
		if (seq[x][y][0]!='X' &&(y%2==0) && (WeirdY[x][y]==0))
		if((logIntens[x][y-1]<log(thval)) &&(logIntens[x][y]<log(thval))) {
			if (seq[x][y][12]=='A' || seq[x][y][12]=='G') 
				nAG++;
			else {nCT++;}
			wieviele++;
		}
	}
	return wieviele;
}


