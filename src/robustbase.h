/* External and interal  API  of  C and Fortran routines in robustbase */

/* C code which includes this, typically includes <R.h> */

/* ./qn_sn.c : */
#define Sint int

void Qn0(double *x, Sint *n, double *res);
void Sn0(double *x, Sint *n, Sint *is_sorted, double *res, double *a2);
/*
 * void Qn    (double *x, Sint *n, Sint *finite_corr, double *res);
 * void Sn    (double *x, Sint *n, Sint *finite_corr, double *res);
*/
void wgt_himed_i(double *x, Sint *n,  Sint *iw, double *res);
void wgt_himed  (double *x, Sint *n, double *w, double *res);


/* ./lmrob.c : */
void R_lmrob_S(double *X, double *y, int *n, int *P,
	       int *nRes, double *scale, double *beta_s,
	       double *C, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, double *rel_tol,
	       int* converged, int *trace_lev);

void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it,
		double *rho_c, double *loss, double *rel_tol,
		int *converged, int *trace_lev);

int F77_NAME(rffastmcd)(
    double *dat, int *n, int *nvar, int *nhalff, int *krep,
    double *initcov, double *initmean, int *inbest, double *det,
    int *weight, int *fit, double *coeff, int *kount, double *adcov,
    int *iseed, int *temp, int *index1, int *index2, double *nmahad,
    double *ndist, double *am, double *am2,
    double *slutn, double *med, double *mad, double *sd,
    double *means, double *bmeans, double *w, double *fv1,
    double *fv2, double *rec, double *sscp1, double *cova1,
    double *corr1, double *cinv1, double *cova2,
    double *cinv2, double *z__, double *cstock, double *mstock,
    double *c1stock, double *m1stock, double *dath,
    double *cutoff, double *chimed);

int F77_NAME(rfltsreg)(
    double *dat, int *n, int *nvar,
    int *nhalff, int *krep, int *inbest, double *objfct,
    int *intercept, int *intadjust, int *nvad, double *datt,
    int *iseed, double *weights, int *temp, int *index1, int *index2,
    double *aw2, double *aw, double *residu, double *y,
    double *nmahad, double *ndist,
    double *am, double *am2, double *slutn, int *jmiss,
    double *xmed, double *xmad, double *a, double *da,
    double *h__, double *hvec, double *c__,
    double *cstock, double *mstock, double *c1stock, double *m1stock,
    double *dath, double *sd, double *means, double *bmeans);
