#include "ilm.h"

extern int matl;
extern int nAG,nCT;
extern int ix[NBS],iy[NBS];
extern int freq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
extern int parfreq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
extern double logIntens[MAX_CHIPDIM][MAX_CHIPDIM];


/*Basically matrix multiplication of matrices A and B which are not explicitely defined:
gsl_matrix M (50 x 50-dimensional) is the product of matrix multiplication of two matrices,
namely X and its transpose Xt, i.e. Xt*X=M. X has the dimension N x 58 where N is the number of
features on the array. 58 unique "properties" or "characteristics" are attributed to each sequence, that
is one column/property and one row per sequence. (For the choice of these properties we refer to
KM Kroll et al, Algorithms of Molecular Biology 2009,4:15)

The restrict the computational effort to a minimum, the matrix X is not explicitely defined
but only the matrix product as it is known beforehand which matrix element holds what
information. Also, due to the symmetry of M only half of the elements are generated, the
other half is "mirrored".
*/

void multAB(gsl_matrix *M, int *xvalueAG, int *yvalueAG, int *xvalueCT, int *yvalueCT)
{
	int x,y;
	int i,j,rAG,rCT;
	double mata,matb;

/* Computation of the first 18 elements, first for purines (A and G) then for pyrimidines (C and T).These elements take the neighboring intensities into account. Each feature has 8 neighbors plus the feature itself=9 times 2 for purines/pyrimidines=18. */
	rAG=rCT=0;
	//  mat[i][j] with 0<= i,j<9
	while (rAG<nAG) {
		x=xvalueAG[rAG];
		y=yvalueAG[rAG];
		gsl_matrix_set(M,0,0,rAG+1);
		for (j=1;j<9;j++) {
			matb=logIntens[x+ix[j]][y+iy[j]];
			gsl_matrix_set(M,0,j,gsl_matrix_get(M,0,j)+matb);
			for (i=1;i<=j;i++) {  
				// only half matrix elements
				mata=logIntens[x+ix[i]][y+iy[i]];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb) );
			}
		}
		rAG++;
	}
	for(j=0;j<9;j++)for(i=0;i<j;i++) {
		gsl_matrix_set(M,j,i,gsl_matrix_get(M,i,j));
		// complete matrix elements
	}

	//  mat[i][j] with 9<= i,j<18
	while (rCT<nCT) {
		x=xvalueCT[rCT];
		y=yvalueCT[rCT];
		gsl_matrix_set(M,9,9,rCT+1);
		for (j=10;j<18;j++) {
			matb=logIntens[x+ix[j-9]][y+iy[j-9]]; 
			gsl_matrix_set(M,9,j,gsl_matrix_get(M,9,j)+matb);
			for (i=10;i<=j;i++) {  
				// only half matrix elements
				mata=logIntens[x+ix[i-9]][y+iy[i-9]];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb) );
			}
		}
		rCT++;
	}
	for(j=9;j<18;j++) for(i=9;i<j;i++) {
		gsl_matrix_set(M,j,i,gsl_matrix_get(M,i,j));
		// complete matrix elements
	}

/*Computation of the next 16 entries reflect the influence of the frequency of particular pairs. There
are 16 pairs (DNA/RNA)*/
	//  mat[i][j] with 18<= j<34, i <=j
	for (rAG=0;rAG<nAG;rAG++) {
		x=xvalueAG[rAG];
		y=yvalueAG[rAG];
		for(j=18;j<34;j++) {
			matb=freq[x][y][j-18];
			gsl_matrix_set(M,0,j,gsl_matrix_get(M,0,j)+matb);
			for (i=1;i<9;i++) {
				mata=logIntens[x+ix[i]][y+iy[i]]; 
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
			for (i=18;i<=j;i++) {
				mata=freq[x][y][i-18];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));;
			}
		}
	}
	for (rCT=0;rCT<nCT;rCT++) {
		x=xvalueCT[rCT];
		y=yvalueCT[rCT];
		for(j=18;j<34;j++) {
			matb=freq[x][y][j-18];
			gsl_matrix_set(M,9,j,gsl_matrix_get(M,9,j)+matb);
			for (i=10;i<18;i++) {
				mata=logIntens[x+ix[i-9]][y+iy[i-9]];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));;
			}
			for (i=18;i<=j;i++) {
				mata=freq[x][y][i-18];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
		}
	}
	for(j=18;j<34;j++) for(i=0;i<j;i++) {
		gsl_matrix_set(M,j,i,gsl_matrix_get(M,i,j));
	}

/*Computation of the last 16 elements, i.e. the weighted pair frequency */
	//  mat[i][j] with j>=34
	for(rAG=0;rAG<nAG;rAG++) {
		x=xvalueAG[rAG];
		y=yvalueAG[rAG];
		for(j=34;j<matl;j++) {
			matb=parfreq[x][y][j-34]*0.1;
			gsl_matrix_set(M,0,j,gsl_matrix_get(M,0,j)+matb);
			for (i=1;i<9;i++) {
				mata=logIntens[x+ix[i]][y+iy[i]];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
			for (i=18;i<34;i++) {
				mata=freq[x][y][i-18];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
			for (i=34;i<=j;i++) {
				mata=parfreq[x][y][i-34]*0.1;
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
		}
	}
	for (rCT=0;rCT<nCT;rCT++) {
		x=xvalueCT[rCT];
		y=yvalueCT[rCT];
		for(j=34;j<matl;j++) {
			matb=parfreq[x][y][j-34]*0.1;
			gsl_matrix_set(M,9,j,gsl_matrix_get(M,9,j)+matb);
			for (i=10;i<18;i++) {
				mata=logIntens[x+ix[i-9]][y+iy[i-9]];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
			for (i=18;i<34;i++) {
				mata=freq[x][y][i-18];
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}
			for (i=34;i<=j;i++) {
				mata=parfreq[x][y][i-34]*0.1;
				gsl_matrix_set(M,i,j,gsl_matrix_get(M,i,j)+(mata*matb));
			}

		}
	}
	for(j=34;j<matl;j++)for(i=0;i<j;i++) {
		gsl_matrix_set(M,j,i,gsl_matrix_get(M,i,j));
	}
}



