/* Header files*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <R.h>

#define MAX_CHIPDIM 1200
#define MDIM 50
#define PROBELENGTH 25
#define PAIRS 16
#define NBS 9



/* Function declarations */
void GetConcsILM(int *Chipdim, int *NIntens, int *zaehler, int *Nfiles, double *threshold, double *Intens, char **Seq, double *Params, double *BackI, double *DGDR, double *DGRR, double *concs,  int *weirdY, double *aLim);

int points(double thresval);

double enval(int wieviele, double ysq, gsl_vector *p, gsl_matrix *mat2);

void multAB(gsl_matrix *M, int *xvalueAG, int *yvalueAG, int *xvalueCT, int *yvalueCT);

void multXty(int *xvalueAG, int *yvalueAG, int *xvalueCT, int *yvalueCT); 

void get_backg(gsl_vector *p, double *BackI, double thresval);

void get_concs(double *concs, double satLim);




