/* Calculation of background intensities and estimation of concentrations in picoMolar*/
#include "ilm.h"

int maxl,matl;
int nAG,nCT;
int ix[NBS]={0,1,-1,1,-1,0,1,-1,0};
int iy[NBS]={0,0,0,-1,-1,-1,1,1,1};

double thval,satLim;

int freq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
int parfreq[MAX_CHIPDIM][MAX_CHIPDIM][PAIRS];
int WeirdY[MAX_CHIPDIM][MAX_CHIPDIM];

char seq[MAX_CHIPDIM][MAX_CHIPDIM][PROBELENGTH+1];

double logIntens[MAX_CHIPDIM][MAX_CHIPDIM];
double Xty[MDIM];
double eta[MAX_CHIPDIM][MAX_CHIPDIM];

/* Comments on Delta G DNA/RNA and RNA/RNA free energy calculation */
/* convention (see line ) 
 * C=0; G=1; A=2; T=3; 
 * pairs numbered from 00 to 16  
 * Hybridization free energies (nucs on DNA strand) at T=37C in kcal/mol
 * Notation: -DeltaG given here!!! 
 *
 *
 * pair no	|00   01   02   03   04   05
 * pair	  	|CC   GC   AC   TC   CG   GG
 * Sugi DNA/RNA	|2.9  2.7  1.1  1.3  1.7  2.1
 * Xia RNA/RNA	|3.26 3.42 2.24 2.35 2.36 3.26
 *
 * 06   07   08   09   10   11
 * AG   TG   CA   GA   AA   TA
 * 0.9  0.9  1.6  1.5  0.2  0.6
 * 2.08 2.11 2.11 2.35 0.93 1.33
 *
 * 12   13   14   15
 * CT   GT   AT   TT
 * 1.8  2.1  0.9  1.0 
 * 2.08 2.24 1.10 0.93
 *
 * References:
 * DNA/RNA hyb.energies: Sugimoto et al. Biochemistry 34,11211 (1995)
 * RNA/RNA: Xia et al. Biochemistry 37, 14719 (1998)
 *   */

double dgdr[MAX_CHIPDIM][MAX_CHIPDIM],dgrr[MAX_CHIPDIM][MAX_CHIPDIM];
double dgsugi37[PAIRS]={2.9,2.7,1.1,1.3,1.7,2.1,
	0.9,0.9,1.6,1.5,0.2,0.6,
	1.8,2.1,0.9,1.0};
double dgxia37[PAIRS]={3.26, 3.42, 2.24, 2.35, 2.36, 3.26,
	2.08, 2.11, 2.11, 2.35, 0.93, 1.33,
	2.08, 2.24, 1.10, 0.93};
double dgdr_init=-3.1;
double dgrr_init=-4.09;
double dgrr_termAU=-0.45;
double dgrr_sym=-0.43;


/* ************************************ BEGIN ************************************** */


void delete_gsl_objects(gsl_matrix *Mat, gsl_matrix *Mat2, gsl_matrix *V, gsl_vector *work, gsl_vector *s)
{
	if(Mat != NULL)
        	gsl_matrix_free(Mat);
        if(Mat2 != NULL)
                gsl_matrix_free(Mat2);
	if(V != NULL)
		gsl_matrix_free(V);
	if(work != NULL)
		gsl_vector_free(work);
	if(s != NULL)
		gsl_vector_free(s);
       return;
}


/*Switch 13th nucleotide to produce correct MM sequence. Usually MM has even y value (y%2==0). Exceptions with odd y-values (y%2!=0). Vector weirdY indicated whether 'normal' or not 'normal' case.  */
void switchNuc(int x, int y, int *weirdY, int i)
{		
	if(((y%2==0) && (weirdY[i]==0)) || ((y%2!=0) && (weirdY[i]==1) )){
		switch(seq[x][y][12]){
			case 'C':
				seq[x][y][12]='G';
				break;
			case 'G':
				seq[x][y][12]='C';
				break;
			case 'A':
				seq[x][y][12]='T';
				break;
			case 'T':
				seq[x][y][12]='A';
				break;
		}
	}
}


/*Calculation of binding free energies btw two strands, i.e. DGDR (DNA/RNA) and DGRR (RNA/RNA) based on the parameters determined by Sugimoto and Xia (see header)*/
void calcDGs(int x, int y, int *weirdY, int i, int letter[])
{
	int j;
	
	if(((y%2!=0) && (weirdY[i]==0)) || ((y%2==0) && (weirdY[i]==1)))
		dgdr[x][y]=dgdr_init;
		dgrr[x][y]=dgrr_init;

	if(seq[x][y][0] == 'A' && seq[x][y][24] == 'T') {
		dgrr[x][y]+=2*dgrr_termAU;
	} else if (seq[x][y][0] == 'A' || seq[x][y][24] == 'T') {
		dgrr[x][y]+=dgrr_termAU;
	}
			
	for(j=0;j<24;j++) {
		freq[x][y][letter[j]+4*letter[j+1]]++;
		parfreq[x][y][letter[j]+4*letter[j+1]]+=(j-11)*(j-12);
		dgdr[x][y]+=dgsugi37[letter[j]+4*letter[j+1]];
		dgrr[x][y]+=dgxia37[letter[j]+4*letter[j+1]];
	}
}



void get_concs(double *concs, double sat)
{
	int x,y;
	int celindx;
	double concents[MAX_CHIPDIM][MAX_CHIPDIM];
	double ImI0,alpha,concval;
	double beta=0.74;

/*only PM's */
	for(y=1;y<maxl;y+=2) for(x=1;x<maxl;x++) {
		concents[x][y]=0.0;
/* In C index starts at 0, hence celindx=ncols*y+x. If celindx=ncols*y+x+1, then min(celindx)=1 like in R */
		celindx=maxl*y+x;
		concval=0.0;

		if(seq[x][y][0] != 'X') {
			alpha=1/(1+exp((-dgrr[x][y]+46.5)*(-0.686)));
			ImI0=exp(logIntens[x][y])-exp(eta[x][y]);
			concval=pow(10,12)*exp(-beta* (dgdr[x][y]+(log(alpha)/beta) ) )* (ImI0)/(sat-ImI0); 
		}
		if(concval >= 0.0){
			concents[x][y] = concval; 
			concs[celindx]=concents[x][y];
		}
	}
}


/*Pseudo main function*/
void GetConcsILM(int *Chipdim, int *NIntens, int *zaehler, int *Nfiles, double *threshold, double *Intens, char **Seq, double *Params, double *BackI, double *DGDR, double *DGRR, double *concs,  int *weirdY, double *aLim)
{
	int x,y;
	int i,j;
	int letter[25];

	int countAG,countCT;
	int npoints;
	int *pindaxAG,*pindayAG,*pindaxCT,*pindayCT;

	double logI,en;
	double ysq=0.0;

	gsl_matrix *Mat; 
	gsl_matrix *V,*Mat2; //output matrix V not transposed! SVD: M=USVt
	gsl_vector *s,*work;

	gsl_error_handler_t *old_handler;

	maxl=*Chipdim;
	matl=MDIM;
	thval=*threshold;
	satLim=*aLim;

	gsl_vector_view params = gsl_vector_view_array(Params,matl);

/* Turn off default error handler*/
	old_handler = gsl_set_error_handler_off();
/* gsl matrix/vector allocation and initialization*/
	Mat = gsl_matrix_calloc(matl,matl);	
	Mat2 = gsl_matrix_calloc(matl,matl);
	V = gsl_matrix_calloc(matl,matl);
	work = gsl_vector_calloc(matl);
	s = gsl_vector_calloc(matl);	
/* Restore default error handler*/
	gsl_set_error_handler(old_handler);

/* Check return values*/
	if(Mat==NULL || Mat2==NULL || V==NULL || work==NULL || s==NULL){
		delete_gsl_objects(Mat,Mat2,V,work,s);
		error("Not enough memory");
	}	


/* Initialize arrays*/
	for(x=0;x<maxl;x++)for(y=0;y<maxl;y++) {
		logIntens[x][y]=0.0;
		WeirdY[x][y]=0;
		for(j=0;j<16;j++) {
			freq[x][y][j]=parfreq[x][y][j]=0;
		}
	}
	for(j=0;j<matl;j++)for (i=0;i<matl;i++) {
		Xty[i] = 0.0; 
	}


/* Index offset! In R index begins with i=1, in C it's i=0! Hence all
 * indices are offset by 1, ie
 * instead of x=(i-1)-Ky we have x=i-yK*/
	i=0;	
	while(i<*NIntens) {
		BackI[i]=concs[i]=0.0;

		y=i/maxl;
		x=i-y*maxl;
		
		logIntens[x][y]=log(Intens[i]);
		WeirdY[x][y]=weirdY[i];
		sscanf(Seq[i],"%s",seq[x][y]);

/*Sequences of all un-annotated probes are just 'XXXXXXXXXXXXXXXXXXXXXXXX' - see line 64 in affyilm.R*/
		if(seq[x][y][1]!='X'){
/* Switch 13th nucleotide to produce correct MM sequence. */
			switchNuc(x,y,weirdY,i);
		}
		
		for(j=0;j<25;j++) {
			if	(seq[x][y][j]=='C') letter[j]=0;
			else if	(seq[x][y][j]=='G') letter[j]=1;
			else if	(seq[x][y][j]=='A') letter[j]=2;
			else  			    letter[j]=3;
		}

		if(seq[x][y][1] != 'X'){
/*Calculate free energy DeltaG for hybridization btw DNA/RNA strands (DGDR) and RNA/RNA (DGRR) (only) for PM's*/
			calcDGs(x,y,weirdY,i,letter);
		}
		DGDR[i]=dgdr[x][y];
		DGRR[i]=dgrr[x][y];
		i++;
	}
	
	npoints=points(thval);
	Rprintf("%d features used for parameter estimation\n",npoints);


/* Allocate memory for MM x and y-values of points fulfilling conditions*/	
	pindaxAG = (int *) R_alloc(nAG,sizeof(int));//middle A or G
	pindayAG = (int *) R_alloc(nAG,sizeof(int));
	pindaxCT = (int *) R_alloc(nCT,sizeof(int));//middle C or T
	pindayCT = (int *) R_alloc(nCT,sizeof(int));

	countAG=countCT=0;
	//loops only over MM's
	for (x=4;x<maxl-4;x++) for (y=4;y<maxl-4;y+=2) {
		if ((seq[x][y][0]!='X') && (y%2==0) && (WeirdY[x][y]==0)) 
		if((logIntens[x][y-1]<log(thval)) && (logIntens[x][y]<log(thval))) {
			logI=logIntens[x][y];
			ysq+=logI*logI;
			if (seq[x][y][12]=='A' || seq[x][y][12]=='G') {
				*(pindaxAG+countAG)=x;
				*(pindayAG+countAG)=y;
				countAG++;
			} else {
				*(pindaxCT+countCT)=x;
				*(pindayCT+countCT)=y;
				countCT++;
			}
		}
	}

/* Matrix multiplication*/
	multAB(Mat,pindaxAG,pindayAG,pindaxCT,pindayCT);
	multXty(pindaxAG,pindayAG,pindaxCT,pindayCT);

	gsl_vector_const_view xTy = gsl_vector_const_view_array(Xty,matl);

/* Copy of Mat - needed for en-value (residual) calculation*/
	gsl_matrix_memcpy(Mat2,Mat);

/* SVD (Singular Value Decomposition) */
	//Mat becomes U
	gsl_linalg_SV_decomp(Mat,V,s,work);
//	Rprintf("SVD done \n");

/* solve Mp=Xty for p with M=USVt */
	gsl_linalg_SV_solve(Mat,V,s,&xTy.vector,&params.vector);

/* Calculation of residual en=sum[Imm(x,y)-eta(x,y)]^2=min */
	en=enval(npoints, ysq, &params.vector, Mat2);
	Rprintf("Residual value is %lf\n",en);

/* Calculation of ln(background intensities) (eta) for each feature */
	get_backg(&params.vector,BackI,thval);
	Rprintf("Background intensities calculated \n");

/* Calculation of concentrations */
	Rprintf("Concentrations calculated\n");
	get_concs(concs,satLim);
/*free memory*/
	delete_gsl_objects(Mat,Mat2,V,work,s);
	Rprintf("Done \n");
	Rprintf("------------------ \n");

}









