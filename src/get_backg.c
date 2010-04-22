#include "ilm.h"

extern int maxl;
extern int freq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
extern int parfreq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
extern double logIntens[MAX_CHIPDIM][MAX_CHIPDIM];
extern double eta[MAX_CHIPDIM][MAX_CHIPDIM];
extern char seq[MAX_CHIPDIM][MAX_CHIPDIM][PROBELENGTH+1];

/* Background calculation at MM positions (even y-values)
 * There are in total 50 parameters p[1] to p[50]:
 *
 * A) Local dependence: Each MM has 8 neighbors. Parameters are categorized according to 12th nucleotide, ie either (A||G) or (C||T)
 * 	1) 12th nucleotide A||G: p[1]-p[9]
 * 	where p[1] is a constant associated to MM under investigation and
 * 	p[2]-p[9] are associated to neighboring sequences
 * 	2) 12th nucleotide C||T: p[10]-p[18]
 * 	with p[10] as constant and
 * 	- p[11]-p[18] neighbors
 *
 * B) Frequency of nucleotide pairs: p[19]-p[33] (--> in analogy to Sugimoto free energies)
 * 
 * C) Weighted pair frequency: p[34]-p[50]
*/

void get_backg(gsl_vector *p, double *BackI, double thresval)
{
	int i;
	int x,y;
	int celindx;
	double backval;
	double Ixp1y,Ixm1y,Ixp1ym1,Ixm1ym1,Ixym1,Ixyp1,Ixm1yp1,Ixp1yp1;

	for(y=2;y<maxl;y+=2) for(x=1;x<maxl-1;x++) {
/* In C index starts at 0, hence celindx=ncols*y+x. If celindx=ncols*y+x+1, then min(celindx)=1 like in R */
		celindx=maxl*y+x;
		backval=0.0;
		if(seq[x][y][0]!='X') {
			//centered MM: p[1] and p[10] constants
			//Direct neighbors 
			Ixp1y   = logIntens[x+1][y  ];//right - p[2],p[11]
			Ixm1y   = logIntens[x-1][y  ];//left - p[3],p[12]

			//Bottom neighbors 
			Ixp1ym1 = logIntens[x+1][y-1];//right bottom - p[4],p[13]
			Ixm1ym1 = logIntens[x-1][y-1];//left bottom - p[5],p[14]
			Ixym1   = logIntens[x  ][y-1];//centered=cor. PM - p[6],p[15]
		
			//Top neighbors
			Ixp1yp1 = logIntens[x+1][y+1];//right top - p[7],p[16] 
			Ixm1yp1 = logIntens[x-1][y+1];//left top - p[8],p[17]
			Ixyp1	= logIntens[x  ][y+1];//centered -p[9],p[18]

/* Set threshold for max intensity allowed for neighboring probes to same threshold as set to probes used to calculate parameters*/
			if(Ixp1y   > log(thresval)) Ixp1y   = log(thresval);
			if(Ixm1y   > log(thresval)) Ixm1y   = log(thresval);

			if(Ixp1ym1 > log(thresval)) Ixp1ym1 = log(thresval);
			if(Ixm1ym1 > log(thresval)) Ixm1ym1 = log(thresval);
			if(Ixym1   > log(thresval)) Ixym1   = log(thresval);
			
			if(Ixp1yp1 > log(thresval)) Ixp1yp1   = log(thresval);
			if(Ixm1yp1 > log(thresval)) Ixm1yp1   = log(thresval);
			if(Ixyp1   > log(thresval)) Ixyp1     = log(thresval);

			if(seq[x][y][12]=='A' || seq[x][y][12]=='G') {
				backval=gsl_vector_get(p,0);
				backval+=gsl_vector_get(p,1)*Ixp1y;
				backval+=gsl_vector_get(p,2)*Ixm1y;

				backval+=gsl_vector_get(p,3)*Ixp1ym1;
				backval+=gsl_vector_get(p,4)*Ixm1ym1;
				backval+=gsl_vector_get(p,5)*Ixym1;
				
				backval+=gsl_vector_get(p,6)*Ixp1yp1;
				backval+=gsl_vector_get(p,7)*Ixm1yp1;
				backval+=gsl_vector_get(p,8)*Ixyp1;
			}
			else if(seq[x][y][12]=='C' || seq[x][y][12]=='T') {
  				backval=gsl_vector_get(p,9);
				backval+=gsl_vector_get(p,10)*Ixp1y;
				backval+=gsl_vector_get(p,11)*Ixm1y;

				backval+=gsl_vector_get(p,12)*Ixp1ym1;
				backval+=gsl_vector_get(p,13)*Ixm1ym1;
				backval+=gsl_vector_get(p,14)*Ixym1;

				backval+=gsl_vector_get(p,15)*Ixp1yp1;
				backval+=gsl_vector_get(p,16)*Ixm1yp1;
				backval+=gsl_vector_get(p,17)*Ixyp1;
			}

			for(i=0;i<16;i++) {
				backval+=gsl_vector_get(p,18+i)*freq[x][y][i];
			}
			for(i=0;i<16;i++) {
				backval+=gsl_vector_get(p,34+i)*parfreq[x][y][i]*0.1;
			}
		}
		eta[x][y]   = backval;//at MM position
		eta[x][y-1] = backval;//copy value to PM position
		BackI[celindx]=exp(eta[x][y]);// MM
		BackI[celindx-maxl]=exp(eta[x][y-1]);//PM
	}
}

