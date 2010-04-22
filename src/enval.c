#include "ilm.h"

extern int matl;
extern double Xty[MDIM];

/* Calculation of the residual value of linear least square estimation */
double enval(int wieviele, double ysq, gsl_vector *p, gsl_matrix *mat2)
{
	int i,j,k;
	double res,res1,res2;

	res=res1=res2=0.0;
	for(i=0;i<matl;i++) for(j=0;j<matl;j++){
		res1 += gsl_vector_get(p,i)*gsl_matrix_get(mat2,i,j)*gsl_vector_get(p,j);
	}
	for(k=0;k<matl;k++){
		res2 += Xty[k]*gsl_vector_get(p,k);
	}
	res=res1-2*res2+ysq;
	res /= wieviele;

	return res;
}

