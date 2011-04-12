#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "ilm.h"

static R_NativePrimitiveArgType GetConcsILM_t[14]={INTSXP,INTSXP,INTSXP,INTSXP,REALSXP,REALSXP,STRSXP,REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,REALSXP};

static const R_CMethodDef cMethods[] = {
	{"GetConcsILM", (DL_FUNC) &GetConcsILM, 14, GetConcsILM_t},  
	{NULL, NULL, 0}
};


void R_init_affyILM(DllInfo *info)
{
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}

void R_unload_affyILM(DllInfo *info) {
	/* Here could be code to release resources. */
}
