#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

#include "robustbase.h"


#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _t)/sizeof(name ## _t[0]), name ##_t}


static R_NativePrimitiveArgType Qn0_t[] = {
    REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType Sn0_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType mc_C_t[] = {
    REALSXP, INTSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType wgt_himed_i_t[] = {
    REALSXP, INTSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType wgt_himed_t[] = {
    REALSXP, INTSXP, REALSXP, REALSXP
};


static R_NativePrimitiveArgType R_lmrob_S_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    /* rrhoc */ REALSXP, INTSXP, REALSXP,
    /* best_r */ INTSXP, INTSXP, INTSXP,
    /* K_s */ INTSXP, INTSXP, REALSXP,
    /* converged */ LGLSXP, INTSXP
};

static R_NativePrimitiveArgType R_lmrob_MM_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP,
    /* beta_initial */ REALSXP, REALSXP,
    /* beta_m */ REALSXP, REALSXP,
    /* max_it */ INTSXP, REALSXP, INTSXP,
    /* loss */ REALSXP, REALSXP, LGLSXP, INTSXP
};

static R_NativePrimitiveArgType R_psifun_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_chifun_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_wgtfun_t[] = {
    REALSXP, REALSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_find_D_scale_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    /* c */ REALSXP, INTSXP, INTSXP, REALSXP,
    /* max_k */ INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_calc_fitted_t[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, 
    INTSXP, INTSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(Qn0),
    CDEF(Sn0),
    CDEF(mc_C),
    CDEF(wgt_himed_i),
    CDEF(wgt_himed),
    CDEF(R_lmrob_S),
    CDEF(R_lmrob_MM),
    CDEF(R_psifun),
    CDEF(R_chifun),
    CDEF(R_wgtfun),
    CDEF(R_find_D_scale),
    CDEF(R_calc_fitted),
    {NULL, NULL, 0}
};

/* static R_CallMethodDef CallEntries[] = { */
/*     {NULL, NULL, 0} */
/* }; */


static R_FortranMethodDef FortEntries[] = {
    {"rffastmcd", (DL_FUNC) &F77_SUB(rffastmcd), 48},/* ./rffastmcd.f */
    {"rfltsreg",  (DL_FUNC) &F77_SUB(rfltsreg), 42}, /* ./rfltsreg.f */
    {NULL, NULL, 0}
};

void R_init_robustbase(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL/*CallEntries*/,
		       FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
