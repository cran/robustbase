
#include <R_ext/Rdynload.h>
#include "robustbase.h"


#define CDEF(name)  {#name, (DL_FUNC) &name, sizeof(name ## _typ)/sizeof(name ## _typ[0]), name ##_typ}
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static R_NativePrimitiveArgType Qn0_typ[] = { // ./qn_sn.c
    REALSXP, INTSXP, REALSXP, INTSXP, REALSXP
};

static R_NativePrimitiveArgType Sn0_typ[] = {
    REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};

static R_NativePrimitiveArgType mc_C_typ[] = {
    REALSXP, INTSXP, REALSXP, INTSXP, REALSXP, LGLSXP
};

static R_NativePrimitiveArgType wgt_himed_i_typ[] = {
    REALSXP, INTSXP, INTSXP, REALSXP
};
static R_NativePrimitiveArgType wgt_himed_typ[] = {
    REALSXP, INTSXP, REALSXP, REALSXP
};


static R_NativePrimitiveArgType R_lmrob_S_typ[] = {
    REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    /* rrhoc */ REALSXP, INTSXP, REALSXP,
    /* best_r */ INTSXP, INTSXP, INTSXP,
    /* K_s */ INTSXP, INTSXP, INTSXP,
    /* rel_tol*/ REALSXP, REALSXP, REALSXP, REALSXP,
    /* converged */ LGLSXP, INTSXP, INTSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_lmrob_MM_typ[] = {
    REALSXP, REALSXP, INTSXP, INTSXP,
    /* beta_initial */ REALSXP, REALSXP,
    /* beta_m */ REALSXP, REALSXP,
    /* max_it */ INTSXP, REALSXP, INTSXP,
    /* loss */ REALSXP, REALSXP, LGLSXP, INTSXP
};

static R_NativePrimitiveArgType R_find_D_scale_typ[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
    /* c */ REALSXP, INTSXP, INTSXP, REALSXP,
    /* max_k */ INTSXP, LGLSXP
};

static R_NativePrimitiveArgType R_calc_fitted_typ[] = {
    REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP,
    INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_lmrob_M_S_typ[] = {
    REALSXP, REALSXP, REALSXP, REALSXP,
    /* nn */ INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
    /* scale */ REALSXP, REALSXP, REALSXP,
    /* rho_c */ REALSXP, INTSXP, REALSXP,
    /* K_m_s */ INTSXP, INTSXP,
    /* rel_tol */ REALSXP, REALSXP, REALSXP, REALSXP,
    /* converged */ LGLSXP, INTSXP,
    /* orthogonalize */ LGLSXP, LGLSXP, LGLSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_subsample_typ[] = {
    REALSXP, REALSXP, INTSXP, INTSXP,
    REALSXP, INTSXP, INTSXP, INTSXP,
    REALSXP, REALSXP, INTSXP,
    REALSXP, REALSXP, INTSXP, INTSXP,
    INTSXP, LGLSXP, INTSXP, INTSXP, REALSXP,
    LGLSXP
};

static const R_CMethodDef CEntries[]  = {
    CDEF(Qn0),
    CDEF(Sn0),
    CDEF(mc_C),
    CDEF(wgt_himed_i),
    CDEF(wgt_himed),
    CDEF(R_lmrob_S),
    CDEF(R_lmrob_MM),
    CDEF(R_find_D_scale),
    CDEF(R_calc_fitted),
    CDEF(R_lmrob_M_S),
    CDEF(R_subsample),
    {NULL, NULL, 0}
};

static R_CallMethodDef CallEntries[] = {
    CALLDEF(R_rho_inf, 2), // -> lmrob.c
    CALLDEF(R_psifun, 4),
    CALLDEF(R_chifun, 4),
    CALLDEF(R_wgtfun, 3),
    CALLDEF(R_wgt_flex, 3), // -> rob-utils.c
    CALLDEF(R_rowMedians, 5),// -> rowMedians.c [Biobase also has rowQ for quantiles]
    {NULL, NULL, 0}
};


static R_FortranMethodDef FortEntries[] = {
    {"rffastmcd", (DL_FUNC) &F77_SUB(rffastmcd), 49},/* ./rffastmcd.f */
    {"rfltsreg",  (DL_FUNC) &F77_SUB(rfltsreg), 41}, /* ./rfltsreg.f */
    {"rllarsbi",  (DL_FUNC) &F77_SUB(rllarsbi), 18}, /* ./rllarsbi.f */
    {NULL, NULL, 0}
};

void R_init_robustbase(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, FortEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);

    R_RegisterCCallable("robustbase", "R_psifun", (DL_FUNC) &R_psifun);
    R_RegisterCCallable("robustbase", "R_chifun", (DL_FUNC) &R_chifun);
    R_RegisterCCallable("robustbase", "R_wgtfun", (DL_FUNC) &R_wgtfun);
    R_RegisterCCallable("robustbase", "rho",  (DL_FUNC) &rho);
    R_RegisterCCallable("robustbase", "psi",  (DL_FUNC) &psi);
    R_RegisterCCallable("robustbase", "psip", (DL_FUNC) &psip);
    R_RegisterCCallable("robustbase", "psi2", (DL_FUNC) &psi2);
    R_RegisterCCallable("robustbase", "wgt",  (DL_FUNC) &wgt);
    R_RegisterCCallable("robustbase", "rho_inf",  (DL_FUNC) &rho_inf);
    R_RegisterCCallable("robustbase", "normcnst", (DL_FUNC) &normcnst);
}
