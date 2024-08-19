/* -*- mode: c; kept-new-versions: 40; kept-old-versions: 40 -*-
 * Indentation (etc) style:  C-c . gnu */

/* file lmrob.c
 * was	roblm/src/roblm.c - version 0.6	 by Matias Salibian-Barreras

 * Includes the stable correct asymptotic variance estimators
 * of Croux, Dhaene, Hoorelbeke
 * Includes the fast-s algorithm
 */

/* Robust MM regression estimates *
 * ------------------------------ */

/* comment code *
 *
 * adapt other sampler <<<<<<<<<<  R's random number generator !!!!

 * replace abort for too many singular resamples by
 * returning the number of singular ones
 */

/* MM:
   -  Done:  fast_s[_large_n]() both had FIXED seed (= 37),
	     and effectively discarded the seed_rand argument below
   -  Done:  drop 'register' : today's compilers do optimize well!
   -  Done:  using Calloc() / Free() instead of malloc()/free()
*/

/* Calls [BLAS 2]   DGEMV("N", ..) == dgemv("N", ..),   which computes

        y := alpha * A*x  +  beta* y,

   where here, always, alpha = -1 ,  beta = 1 ('dmone' or 'done', respectively),  i.e.,  we do

        y := -A*x + y = y - A*x

   Now, call DGEMV =~= F77_CALL(dgemv)  only via macro  get_Res_Xb()
*/

/* kollerma:
   Added alternative psi functions callable via psifun, chifun and
   wgtfun. ipsi is used to distinguish between the different types:
   0: huber,    1: biweight="bisquare",  2: GaussWeight="welsh",
   3: optimal,  4: hampel,
   5: GGW="ggw" (Generalized Gauss Weight),
   6: LQQ="lqq"=lin.quadr.qu../ piecewise linear psi'()

   The ipsi argument works for the S-estimator as well as for the MM-estimator.

   - Added implementation of M-S algorithm.

   - Modified subsampling behaviour: avoiding singular resamples by using
     customized LU decomposition.

   - Replaced C style matrices with Fortran style matrices, with as little
     copying as possible.

   - Using LAPACK's DGELS instead of local lu() decomposition.

   - Code clean up: removed all subroutines that were unused.
*/

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
# define FCONE
#endif

#include <Rmath.h>
#include <R_ext/Applic.h>
#include <float.h>

#include "robustbase.h"
//-> <R.h>, <Rinternals.h>  -> XLENGTH, R_xlen_t


/* these  will also move to "lmrob.h" ---
 *  but first make many of these 'static' <<< FIXME!
 */
void fast_s_large_n(double *X, double *y, double s_y,
		    int n, int p, int nRes, int *max_it_scale, double *res,
		    int groups, int n_group,
		    int K, int *max_k,
		    double rel_tol, double inv_tol, double scale_tol, double zero_tol,
		    int *converged,
		    int best_r, double bb, const double rrhoc[], int ipsi,
		    double *bbeta, double *sscale, int trace_lev, int mts, Rboolean ss);

void fast_s(double *X, double *y, double s_y,
	    int n, int p, int nResample, int *max_it_scale, double *res,
	    int K, int *max_k,
	    double rel_tol, double inv_tol, double scale_tol, double zero_tol,
	    int *converged,
	    int best_r, double bb, const double rrhoc[], int ipsi,
	    double *bbeta, double *sscale, int trace_lev, int mts, Rboolean ss);

Rboolean rwls(const double X[], const double y[], int n, int p,
	 double *estimate, double *i_estimate,
	 double *resid, double *loss,
	 double scale, double epsilon,
	 int *max_it, const double rho_c[], const int ipsi, int trace_lev);

double norm2 (const double x[], int n);
double norm (const double x[], int n);
double norm1(const double x[], int n);
double mean_abs(const double y[], int n) { return norm1(y, n) / n ; }

double norm_diff2(const double x[], const double y[], int n);
double norm_diff (const double x[], const double y[], int n);
double norm1_diff(const double x[], const double y[], int n);

/* moved to robustbase.h
 *
 * double normcnst(const double c[], int ipsi);
 * double rho_inf (const double c[], int ipsi);

 * double rho(double x, const double c[], int ipsi);
 * double psi(double x, const double c[], int ipsi);
 * double psip(double x, const double c[], int ipsi);// psi'
 * double psi2(double x, const double c[], int ipsi);// psi''
 * double wgt(double x, const double c[], int ipsi);
 */

double rho_huber(double x, const double c[]);
double psi_huber(double x, const double c[]);
double psip_huber(double x, const double c[]);
double psi2_huber(double x, const double c[]);
double wgt_huber(double x, const double c[]);

double rho_biwgt(double x, const double c[]);
double psi_biwgt(double x, const double c[]);
double psip_biwgt(double x, const double c[]);
double psi2_biwgt(double x, const double c[]);
double wgt_biwgt(double x, const double c[]);

double rho_gwgt(double x, const double c[]);
double psi_gwgt(double x, const double c[]);
double psip_gwgt(double x, const double c[]);
double wgt_gwgt(double x, const double c[]);

double rho_opt(double x, const double c[]);
double psi_opt(double x, const double c[]);
double psip_opt(double x, const double c[]);
double wgt_opt(double x, const double c[]);

double rho_hmpl(double x, const double c[]);
double psi_hmpl(double x, const double c[]);
double psip_hmpl(double x, const double c[]);
double psi2_hmpl(double x, const double c[]);
double wgt_hmpl(double x, const double c[]);

double rho_ggw(double x, const double c[]);
void psi_ggw_vec(double *x, int n, void *k);
double psi_ggw(double x, const double c[]);
double psip_ggw(double x, const double c[]);
double wgt_ggw(double x, const double c[]);

double rho_lqq(double x, const double c[]);
double psi_lqq(double x, const double c[]);
double psip_lqq(double x, const double c[]);
double psi2_lqq(double x, const double c[]);
double wgt_lqq(double x, const double c[]);

double sum_rho_sc(const double r[], double scale, int n, int p,
		  const double c[], int ipsi);
int refine_fast_s(const double X[], double *wx, const double y[], double s_y, double *wy,
		  double *weights, int n, int p, double *res,
		  double *work, int lwork,
		  const double beta_cand[], double *beta_j,
		  Rboolean *conv, int kk, double rel_tol, double zero_tol,
		  int trace_lev,
		  double b, const double rrhoc[], int ipsi, double initial_scale,
		  double *beta_ref, double *scale);

void m_s_subsample(double *X1, double *y, int n, int p1, int p2,
		   int nResample, int max_it_scale,
		   double rel_tol, double inv_tol, double scale_tol, double zero_tol, double bb,
		   const double rrhoc[], int ipsi, double *sscale, int trace_lev,
		   double *b1, double *b2, double *t1, double *t2,
		   double *y_tilde, double *res, double *x1, double *x2,
		   int *NIT, int *K, int *KODE, double *SIGMA, double BET0,
		   double *SC1, double *SC2, double *SC3, double *SC4, int mts, Rboolean ss);

Rboolean m_s_descent(double *X1, double *X2, double *y,
		 int n, int p1, int p2, int K_m_s, int max_k, int max_it_scale,
		 double rel_tol, double scale_tol, double bb, const double rrhoc[], int ipsi,
		 double *sscale, int trace_lev,
		 double *b1, double *b2, double *t1, double *t2,
		 double *y_tilde, double *res, double *res2, double *x1, double *x2,
		 int *NIT, int *K, int *KODE, double *SIGMA, double BET0,
		 double *SC1, double *SC2, double *SC3, double *SC4);

int subsample(const double x[], const double y[], int n, int m,
	      double *beta, int *ind_space, int *idc, int *idr,
	      double *lu, double *v, int *p,
	      double *Dr, double *Dc, int rowequ, int colequ,
	      Rboolean sample, int mts, Rboolean ss, double tol_inv, Rboolean solve);

int fast_s_with_memory(double *X, double *y, double s_y, double *res,
		       int nn, int pp, int nRes, int *max_it_scale, int K, int *max_k,
		       double rel_tol, double inv_tol, double scale_tol, double zero_tol,
		       int trace_lev,
		       int best_r, double bb, const double rrhoc[], int ipsi,
		       int mts, Rboolean ss,
		       // output ==>
		       double **best_betas, double *best_scales);

/* for "tracing" only : */
void disp_mat(const double **a, int n, int m);
void disp_vec(const double a[], int n);
void disp_veci(const int a[], int n);

double kthplace(double *, int, int);

int find_max(const double a[], int n);

double find_scale(const double r[], double b, const double rrhoc[], int ipsi,
		  double initial_scale, int n, int p,
		  int* iter, // input: max_iter,  output: #{iterations used}
		  double scale_tol,
		  Rboolean trace);

double MAD(const double a[], int n, double center, double *tmp, double *tmp2);

void zero_mat(double **a, int n, int m);

// used in many F77_CALL()s :
static const double dmone = -1., done = 1.;
static const int one = 1;

#define get_Res_Xb(_Mp,_Np,_A_,_x_,_Y_)					\
    F77_CALL(dgemv)("N", &_Mp, &_Np, &dmone, _A_, &_Mp, _x_, &one, &done, \
                    /*--> */ _Y_, &one FCONE)

// additionally using info (is *modified), trace_lev:
#define INIT_WLS(_X_, _y_, _n_, _p_)				\
    int lwork = -1;						\
    double work0;						\
								\
    /* Determine optimal block size for work array */		\
    F77_CALL(dgels)("N", &_n_, &_p_, &one, _X_, &_n_, _y_,	\
		    &_n_, &work0, &lwork, &info FCONE);		\
    if (info) {							\
	warning(_("DGELS could not determine optimal block size, using minimum")); \
	lwork = 2*_p_;						\
    } else							\
	lwork = (int)work0;					\
								\
    if (trace_lev >= 4)						\
	Rprintf(" Optimal block size for DGELS: %d\n", lwork);	\
								\
    /* allocate */						\
    double *work = (double *) R_alloc(lwork, sizeof(double)),	\
	*weights = (double *) R_alloc(_n_,   sizeof(double));

/* Solve weighted LS, called in a loop, from rwls(), fast_refine_s() & m_s_descent():
   _x_ is a "work array" only,  input(_X, _y_, _wts_);  output _beta_ (_y_ is *modified*)
     also uses (trace_level, info, one, work, lwork) */
#define FIT_WLS(_X_, _wts_, _x_, _y_, _n_, _p_, _beta_)		\
    /* add weights to _y_ and _x_ */				\
    for (int j=0; j<_n_; j++) {					\
	double wtmp = sqrt(_wts_[j]);				\
	_y_[j] *= wtmp;						\
	for (int k=0; k<_p_; k++)				\
	    _x_[_n_*k+j] = _X_[_n_*k+j] * wtmp;			\
    }								\
    /* solve weighted least squares problem */			\
    F77_CALL(dgels)("N", &_n_, &_p_, &one, _x_, &_n_, _y_,	\
		    &_n_, work, &lwork, &info FCONE);		\
    if (info) {							\
	if (info < 0) {						\
	    error(_("DGELS: illegal %i-th argument"), -info);	\
	} else {						\
	    if (trace_lev >= 4) {				\
		Rprintf(" Robustness weights in failing step: "); \
		disp_vec(_wts_, _n_);				\
	    }							\
	    error(_("DGELS: weighted design matrix not of full rank (column %d).\n" \
	            "Use control parameter 'trace.lev = 4' to get diagnostic output"), info); \
	}							\
    }								\
    COPY(_y_, _beta_, _p_)

// *modifying* 'info':
#define SETUP_EQUILIBRATION(_n_, _p_, _X_, _large_n_)           \
    /* equilibration of matrix _X_                          */  \
    /* solve (Dr X Dc) b = Dr y with beta = Dc b instead of */  \
    /*            X beta = y                                */  \
    /* see Demmel (1997) APPLIED NUMERICAL LINEAR ALGEBRA   */  \
    /*     Section 2.5.2 Equilibration                      */  \
    double *Dr, *Dc, *Xe, rowcnd, colcnd, amax;			\
    int rowequ = 0 , colequ = 0;				\
    Dr =        (double *) R_alloc(_n_,     sizeof(double));    \
    Dc =        (double *) R_alloc(_p_,     sizeof(double));    \
    Xe =        (double *) R_alloc(_n_*_p_, sizeof(double));    \
    COPY(_X_, Xe, _n_*_p_);                                     \
    F77_CALL(dgeequ)(&_n_, &_p_, Xe, &_n_, Dr, Dc, &rowcnd,	\
    		     &colcnd, &amax, &info);                    \
    if (info) {                                                 \
	if (info < 0) {                                         \
	    error(_("DGEEQU: illegal %i-th argument"), - info); \
	} else if (info > _n_) {                                \
	    if (_large_n_) {                                    \
	        error(_("Fast S large n strategy failed. Use control parameter 'fast.s.large.n = Inf'.")); \
	    } else {						\
                error(_("DGEEQU: column %i of the design matrix is exactly zero."), info - _n_); \
	    }                                                   \
	} else if (_p_ > 1) { /* p = 1 & x[1] = 0 is "ok" */ 	\
	/* FIXME: replace dgeequ by our own version */          \
	/* that does not treat this as error */                 \
	    warning(_(" Skipping design matrix equilibration (DGEEQU): row %i is exactly zero."), info); \
	}                                                       \
    } else {							\
        /* scale _X_ */                                         \
        char equed = '?'; /* init to prevent warning */		\
	F77_CALL(dlaqge)(&_n_, &_p_, Xe, &_n_, Dr, Dc, &rowcnd,	\
			 &colcnd, &amax, &equed FCONE);		\
        rowequ = equed == 'B' || equed == 'R';                  \
	colequ = equed == 'B' || equed == 'C';                  \
    }

#define SETUP_SUBSAMPLE(_n_, _p_, _X_, _large_n_)		\
    /* (Pointers to) Arrays - to be allocated */		\
    int *ind_space = (int *)    R_alloc(_n_,     sizeof(int)),		\
	*idc =       (int *)    R_alloc(_n_,     sizeof(int)),		\
	*idr =       (int *)    R_alloc(_p_,     sizeof(int)),		\
	*pivot =     (int *)    R_alloc(_p_-1,   sizeof(int));		\
    double *lu =     (double *) R_alloc(_p_*_p_, sizeof(double)),	\
	*v =         (double *) R_alloc(_p_,     sizeof(double));	\
    SETUP_EQUILIBRATION(_n_, _p_, _X_, _large_n_);

#define COPY(from, to, len) Memcpy(to, from, len)
/* In theory BLAS should be fast, but this seems slightly slower,
 * particularly for non-optimized BLAS :*/
/* #define COPY(FROM, TO, _p_) \ */
/*     F77_CALL(dcopy)(&_p_, FROM, &one, TO, &one);  */

// NB: INFI is now combined with s_y (or SIGMA) in order to be equivariant (!) :
#define INFI 1e+20

/* Called from R's  lmrob.S() in ../R/lmrob.MM.R,
 * help() in ../man/lmrob.S.Rd, this function computes an S-regression estimator
             ~~~~~~~~~~~~~~~~~
 */
void R_lmrob_S(double *X,
	       double *y, // y[] on input;  residuals[] on output
	       int *n, int *P,
	       int *nRes, // = nResample ( = 500, by default)
	       double *scale, double *beta_s,
	       double *rrhoc, int *iipsi, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, int *max_it_scale,
	       double *rel_tol, double *inv_tol,
	       double *scale_tol, // <- new, was hardwired to EPS_SCALE := 1e-10
	       double *zero_tol, // <- very new, was hardwired to EPS_SCALE := 1e-10
	       int *converged, int *trace_lev, int *mts,
	       int *ss, // subsampling
	       int *cutoff // defining "large n" <==> using fast_s_large_n()
    )
{
    /* best_r = 't' of Salibian-Barrera_Yohai(2006),
     *	      = no. of best candidates to be iterated further ("refined")
     *        =  2, by default
     */
    if (*nRes > 0) {
	double *res = (double *) R_alloc(*n, sizeof(double)); // residuals
	double s_y = mean_abs(y, *n);
	if (*n > *cutoff) {
	    if(*trace_lev > 0)
		Rprintf("lmrob_S(n = %d, nRes = %d): fast_s_large_n():\n", *n, *nRes);
	    fast_s_large_n(X, y, s_y, *n, *P, *nRes, max_it_scale, res,
			   *Groups, *N_group,
			   *K_s, max_k, *rel_tol, *inv_tol, *scale_tol, *zero_tol, converged,
			   *best_r, *bb, rrhoc, *iipsi, beta_s, scale, *trace_lev, *mts, (Rboolean)*ss);
	} else {
	    if(*trace_lev > 0)
		Rprintf("lmrob_S(n = %d, nRes = %d): fast_s() [non-large n]:\n", *n, *nRes);
	    fast_s(X, y, s_y, *n, *P, *nRes, max_it_scale, res,
		   *K_s, max_k, *rel_tol, *inv_tol, *scale_tol, *zero_tol, converged,
		   *best_r, *bb, rrhoc, *iipsi, beta_s, scale, *trace_lev, *mts, *ss);
	}
	COPY(res, y, *n); // return the 'residuals' in 'y'
    } else { // nRes[] <= 0   <==>   'only.scale = TRUE'
	if(*trace_lev > 0)
	    Rprintf("lmrob_S(nRes = 0, n = %d): --> find_scale(*, scale=%g) only:",
		    *n, *scale);
	*scale = find_scale(y, *bb, rrhoc, *iipsi, *scale, *n, *P,
			    max_it_scale, *scale_tol, *trace_lev >= 3);
	if(*trace_lev > 0) Rprintf(" used %d iterations\n", *max_it_scale);
    }
}

/* Called from R, this function computes an M-S-regression estimator */
// not only called from ../R/lmrob.M.S.R,  but also ../inst/xtraR/m-s_fns.R
//			~~~~~~~~~~~~~~~~	    ~~~~~~~~~~~~~~~~~~~~~~~
void R_lmrob_M_S(double *X1, double *X2, double *y, double *res,
		 int *nn, int *pp1, int *pp2, int *nRes, int *max_it_scale,
		 double *scale, double *b1, double *b2,
		 double *rho_c, int *ipsi, double *bb,
		 int *K_m_s, int *max_k,
		 double *rel_tol, double *inv_tol, double *scale_tol,
		 double *zero_tol,
		 int *converged,
		 int *trace_lev,
		 int *orthogonalize, int *subsample, int *descent,
		 int *mts, int *ss)
{
    /* Initialize (some of the) memory here,
     * so that we have to do it only once */
    int i, n = *nn, p1 = *pp1, p2 = *pp2;

    if(*trace_lev > 0) Rprintf(
	"lmrob_M_S(n = %d, nRes = %d, (p1,p2)=(%d,%d), (orth,subs,desc)=(%d,%d,%d))\n",
	n, *nRes, p1, p2,
	*orthogonalize, *subsample, *descent);
    double
    *t1 =      (double *) R_alloc(n,  sizeof(double)), /* size n needed for rllarsbi */
    *t2 =      (double *) R_alloc(p2, sizeof(double)),
    *ot1 =     (double *) R_alloc(p1, sizeof(double)),
    *oT2 =     (double *) R_alloc(p2*p1, sizeof(double)),
    *y_work =  (double *) R_alloc(n,  sizeof(double)),
    *y_tilde = (double *) R_alloc(n,  sizeof(double)),
    *x1 =      (double *) R_alloc(n*p1, sizeof(double)),
    *x2 =      (double *) R_alloc(n*p2, sizeof(double));
    COPY(y, y_work, n);
    COPY(X2, x2, n*p2);

    /* Variables required for rllarsbi() := L1 / least absolute residuals - estimate */
    int NIT=0, K=0, KODE=0;
    double SIGMA = 0.,
	*SC1 = (double *) R_alloc(n,  sizeof(double)),
	*SC2 = (double *) R_alloc(p1, sizeof(double)),
	*SC3 = (double *) R_alloc(p1, sizeof(double)),
	*SC4 = (double *) R_alloc(p1, sizeof(double));
    const double BET0 = 0.773372647623; /* = pnorm(0.75) */

    /* STEP 1: Orthgonalize X2 and y from X1 */
    if (*orthogonalize) {
	COPY(X1, x1, n*p1);
	F77_CALL(rllarsbi)(x1, y_work, &n, &p1, &n, &n, rel_tol,
			   &NIT, &K, &KODE, &SIGMA, t1, y_tilde, SC1, SC2,
			   SC3, SC4, &BET0);
	COPY(t1, ot1, p1);
	for (i=0; i < p2; i++) {
	    COPY(X1, x1, n*p1);
	    double *ptr = X2+i*n; COPY(ptr, y_work, n);
	    F77_CALL(rllarsbi)(x1, y_work, &n, &p1, &n, &n, rel_tol,
			       &NIT, &K, &KODE, &SIGMA, t1, x2+i*n, SC1, SC2,
			       SC3, SC4, &BET0);
	    ptr = oT2+i*p1; COPY(t1, ptr, p1);
	}
	COPY(y_tilde, y_work, n);
	/* compare with Maronna & Yohai 2000:
	 * y_work and y_tilde now contain \tilde y, ot1 -> t_1,
	 * x2 -> \tilde x2, oT2 -> T_2 */
	if(*trace_lev >= 2) Rprintf("  orthogonalized: SIGMA=%g\n", SIGMA);
    } else {
	SIGMA = mean_abs(y, n);
	if(*trace_lev >= 2) Rprintf("  *no* orthog., SIGMA=mean(|y_i|)= %g\n", SIGMA);
    }

    /* STEP 2: Subsample */
    if (*subsample) {
	m_s_subsample(X1, y_work, n, p1, p2, *nRes, *max_it_scale,
		      *rel_tol, *inv_tol, *scale_tol, *zero_tol, *bb,
		      rho_c, *ipsi, scale, *trace_lev,
		      b1, b2, t1, t2, y_tilde, res, x1, x2,
		      &NIT, &K, &KODE, &SIGMA, BET0,
		      SC1, SC2, SC3, SC4, *mts, *ss);

	if (*scale < 0)
	    error(_("m_s_subsample() stopped prematurely (scale < 0)."));
    }

    /* STEP 3: Transform back */
    if (*orthogonalize) {
	/* t1 = ot1 + b1 - oT2 %*% b2 */
	for(int i=0; i < p1; i++) t1[i] = ot1[i] + b1[i];
	get_Res_Xb(p1,p2, oT2, b2, /*--> */ t1);
	COPY(t1, b1, p1);
	/* restore x2 */
	COPY(X2, x2, n*p2);
    }

    /* update / calculate residuals  res := y - X1*b1 - X2*b2 */
    COPY(y, res, n);
    get_Res_Xb(n,p1, X1, b1, /*--> */ res);
    get_Res_Xb(n,p2, X2, b2, /*--> */ res);

    /* STEP 4: Descent procedure */
    if (*descent) {
	*converged = m_s_descent(
	    X1, X2, y, n, p1, p2, *K_m_s, *max_k, *max_it_scale,
	    *rel_tol, *scale_tol, *bb, rho_c, *ipsi, scale, *trace_lev,
	    b1, b2, t1, t2, y_tilde, res, y_work, x1, x2,
	    &NIT, &K, &KODE, &SIGMA, BET0, SC1, SC2, SC3, SC4);
    }
}

/* This function performs RWLS iterations starting from
 * an S-regression estimator (and associated residual scale).
 * So, in itself, this is ``just'' an M-estimator -- called from R's
 * lmrob..M..fit()  [ ../R/lmrob.MM.R ]
 * ~~~~~~~~~~~~~~~
 * NOTE: rel_tol now controls the *relative* changes in beta,
 *      instead of being hard-wired to EPS = 1e-7 and bounding the
 *	absolute || beta_1 - beta_2 ||
 */
void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it, double *rho_c, int *ipsi,
		double *loss, double *rel_tol, int *converged, int *trace_lev)
{
    /* starting from the S-estimate (beta_initial), use
     * irwls to compute the MM-estimate (beta_m)  */

    if(*trace_lev > 0)
	Rprintf("lmrob_MM(): rwls():\n");

    *converged = (int)rwls(X,y,*n,*P,beta_m, beta_initial, resid, loss,
		      *scale, *rel_tol,
		      max_it, rho_c, *ipsi, *trace_lev);
    if (!converged)
	COPY(beta_initial, beta_m, *P);
}


/* Call subsample() from R, for testing purposes only */
void R_subsample(const double x[], const double y[], int *n, int *m,
		 double *beta, int *ind_space, int *idc, int *idr,
		 double *lu, double *v, int *p,
		 double *_Dr, double *_Dc, int *_rowequ, int *_colequ,
		 int *status, int *sample, int *mts, int *ss, double *tol_inv,
		 int *solve)
{
    int info;

    /*	set the seed */
    GetRNGstate();

    SETUP_EQUILIBRATION(*n, *m, x, 0);

    *status = subsample(Xe, y, *n, *m, beta, ind_space, idc, idr, lu, v, p,
			Dr, Dc, rowequ, colequ, (Rboolean)*sample, *mts, (Rboolean)*ss,
			*tol_inv, (Rboolean)*solve);

    COPY(Dr, _Dr, *n);
    COPY(Dc, _Dc, *m);
    *_rowequ = rowequ;
    *_colequ = colequ;

    PutRNGstate();
}

//---- Psi(), Rho(), Functions-----------------------------------------------------------

SEXP R_psifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_) {
    /*
     * Calculate psi for vectorized x, scaled to get  psi'(0) = 1
     * deriv -1: rho(x)  {*not* normalized}
     * deriv  0: psi(x)  = rho'(x)
     * deriv  1: psi'(x) = rho''(x)   {we always have psip(0) == 1}
     * deriv  2: psi''(x)= rho'''(x)
     */
    int nprot = 1, ipsi = asInteger(ipsi_), deriv = asInteger(deriv_);
    if (isInteger(x_)) {
	x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++;
    }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_)) error(_("Argument '%s' must be numeric or integer"), "cc");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), *cc = REAL(c_);

// put the for() loop *inside* the switch (--> speed for llength >> 1) :
#define for_i_n_NA  for(i = 0; i < n; i++) r[i] = ISNAN(x[i]) ? x[i] :

    switch(deriv) { // our rho() is rho~(), i.e., scaled to max = 1
    case -1:
	if(is_redescender(ipsi)) {
	    double rho_Inf = rho_inf(cc, ipsi);
	    for_i_n_NA rho(x[i], cc, ipsi) * rho_Inf;
	} else { // huber, ..
	    for_i_n_NA rho(x[i], cc, ipsi);
	}
	break;

    case  0: for_i_n_NA psi (x[i], cc, ipsi); break;
    case  1: for_i_n_NA psip(x[i], cc, ipsi); break;
    case  2: for_i_n_NA psi2(x[i], cc, ipsi); break;
    default: error(_("'deriv'=%d is invalid"), deriv);
    }
    UNPROTECT(nprot);
    return res;
}

SEXP R_chifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_) {
    /*
     * Calculate chi for vectorized x, i.e. rho~(.) with rho~(inf) = 1:
     * deriv 0: chi (x)  = \rho(x) / \rho(Inf) =: \rho(x) * nc  == our rho() C-function
     * deriv 1: chi'(x)  = psi(x)  * nc
     * deriv 2: chi''(x) = psi'(x) * nc
     */
    int nprot = 1, ipsi = asInteger(ipsi_), deriv = asInteger(deriv_);
    if (isInteger(x_)) {
	x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++;
    }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_)) error(_("Argument '%s' must be numeric or integer"), "cc");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), *cc = REAL(c_);

    // our rho() is rho~() == chi(), i.e., scaled to max = 1
    double rI = (deriv > 0) ? rho_inf(cc, ipsi) : 0./* -Wall */;

    switch(deriv) {
    case 0: for_i_n_NA rho(x[i], cc, ipsi); break;
    case 1: for_i_n_NA psi (x[i], cc, ipsi) / rI; break;
    case 2: for_i_n_NA psip(x[i], cc, ipsi) / rI; break;
    case 3: for_i_n_NA psi2(x[i], cc, ipsi) / rI; break;
    default: error(_("'deriv'=%d is invalid"), deriv);
    }
    UNPROTECT(nprot);
    return res;
}

SEXP R_wgtfun(SEXP x_, SEXP c_, SEXP ipsi_) {
    /*
     * Calculate wgt(x) = psi(x)/x for vectorized x
     */
    int nprot = 1, ipsi = asInteger(ipsi_);
    if (isInteger(x_)) {
	x_ = PROTECT(coerceVector(x_, REALSXP)); nprot++;
    }
    if (!isReal(x_)) error(_("Argument '%s' must be numeric or integer"), "x");
    if (!isReal(c_)) error(_("Argument '%s' must be numeric or integer"), "cc");
    R_xlen_t i, n = XLENGTH(x_);
    SEXP res = PROTECT(allocVector(REALSXP, n)); // the result
    double *x = REAL(x_), *r = REAL(res), *cc = REAL(c_);

    for_i_n_NA wgt(x[i], cc, ipsi);

    UNPROTECT(nprot);
    return res;
}

#undef for_i_n_NA

SEXP R_rho_inf(SEXP cc, SEXP ipsi) {
    if (!isReal(cc)) error(_("Argument 'cc' must be numeric"));
    if (!isInteger(ipsi)) error(_("Argument 'ipsi' must be integer"));
    return ScalarReal(rho_inf(REAL(cc), INTEGER(ipsi)[0]));
}

double rho_inf(const double k[], int ipsi) {
    /*
     * Compute  \rho(\infty) for psi functions
     * (Note that our C function rho() is "rho~" and has rho(Inf) = 1)
     */
    double c = k[0];

    switch(ipsi) {
    default: error(_("rho_inf(): ipsi=%d not implemented."), ipsi);
    case 0: return(R_PosInf); // huber
    case 1: return(c*c/6.); // biweight
    case 2: return(c*c); // GaussWeight / "Welsh"
    case 3: return(3.25*c*c); // Optimal
    case 4: return(0.5*k[0]*(k[1]+k[2]-k[0])); // Hampel
    case 5: // GGW (Generalized Gauss Weight)
	switch((int)c) {
	default:
	case 0: return(k[4]); break; // k[4] == cc[5] in R -- must be correct!
	case 1: return(5.309853); break;
	case 2: return(2.804693); break;
	case 3: return(0.3748076); break;
	case 4: return(4.779906); break;
	case 5: return(2.446574); break;
	case 6: return(0.4007054); break;
	};
    case 6: // LQQ aka 'lin psip'
	return (k[2]*k[1]*(3*k[1]+2*k[0]) + (k[0]+k[1])*(k[0]+k[1])) / (6.*(k[2]-1.));
    }
} // rho_inf()

double normcnst(const double k[], int ipsi) {
    /*
     * return normalizing constant for psi functions :=  1 / \rho(\infty)
     */

    double c = k[0];

    switch(ipsi) {
    default: error(_("normcnst(): ipsi=%d not implemented."), ipsi);
    case 0: return(0.); // huber {normcnst() should never be used for that!}
    case 1: return(6./(c*c)); // biweight
    case 2: return(1./(c*c)); // GaussWeight / "Welsh"
    case 3: return(1./3.25/(c*c)); // Optimal
    case 4: return(2./(k[0]*(k[1]+k[2]-k[0]))); // Hampel
    case 5: // GGW
	switch((int)c) {
	default:
	case 0: return(1./ k[4]); break; // k[4] == cc[5] in R -- must be correct!
	case 1: return(1./5.309853); break;
	case 2: return(1./2.804693); break;
	case 3: return(1./0.3748076); break;
	case 4: return(1./4.779906); break;
	case 5: return(1./2.446574); break;
	case 6: return(1./0.4007054); break;
	};
    case 6: // LQQ aka 'lin psip'
	return((6*(k[2]-1))/(k[2]*k[1]*(3*k[1]+2*k[0])+(k[0]+k[1])*(k[0]+k[1])));
    }
} // normcnst()


double rho(double x, const double c[], int ipsi)
{
    /*
     * return the correct rho according to ipsi
     * This rho() is normalized to 1, called rho~() or chi() in other contexts
     */
    switch(ipsi) {
    default: error(_("rho(): ipsi=%d not implemented."), ipsi);
    case 0: return(rho_huber(x, c)); // huber
    case 1: return(rho_biwgt(x, c)); // biweight
    case 2: return(rho_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(rho_opt(x, c)); // Optimal
    case 4: return(rho_hmpl(x, c)); // Hampel
    case 5: return(rho_ggw(x, c)); // GGW (Generalized Gauss Weight)
    case 6: return(rho_lqq(x, c)); // LQQ := Linear-Quadratic-Quadratic
	// was LGW := "lin psip" := piecewise linear psi'()
    }
}

double psi(double x, const double c[], int ipsi)
{
    /*
     * return the correct psi according to ipsi
     * this is actually rho' and not psi
     */
    switch(ipsi) {
    default: error(_("psi(): ipsi=%d not implemented."), ipsi);
    case 0: return(psi_huber(x, c)); // huber
    case 1: return(psi_biwgt(x, c)); // biweight
    case 2: return(psi_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(psi_opt(x, c)); // Optimal
    case 4: return(psi_hmpl(x, c)); // Hampel
    case 5: return(psi_ggw(x, c)); // GGW
    case 6: return(psi_lqq(x, c)); // LQQ (piecewise linear psi')
    }
}

double psip(double x, const double c[], int ipsi)
{
    /*
     * return the correct ppsi according to ipsi
     * this is actually rho'' and not psip
     */
    switch(ipsi) {
    default: error(_("psip(): ipsi=%d not implemented."), ipsi);
    case 0: return(psip_huber(x, c)); // huber
    case 1: return(psip_biwgt(x, c)); // biweight
    case 2: return(psip_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(psip_opt(x, c)); // Optimal
    case 4: return(psip_hmpl(x, c)); // Hampel
    case 5: return(psip_ggw(x, c)); // GGW
    case 6: return(psip_lqq(x, c)); // LQQ (piecewise linear psi')
    }
}

double psi2(double x, const double c[], int ipsi)
{
    /* Compute   psi''(x) == rho'''(x)
     */
    switch(ipsi) {
    // default: error(_("psi2: ipsi=%d not implemented."), ipsi);
    case 0: return(psi2_huber(x, c)); // huber
    case 1: return(psi2_biwgt(x, c)); // biweight
    case 4: return(psi2_hmpl(x, c)); // Hampel
    case 6: return(psi2_lqq(x, c)); // LQQ (piecewise linear psi')

    default: error(_("psi2(): ipsi=%d not implemented."), ipsi);
/*
    case 2: return(psi2_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(psi2_opt(x, c)); // Optimal
    case 5: return(psi2_ggw(x, c)); // GGW
*/
    }
}

double wgt(double x, const double c[], int ipsi)
{
    /*
     * return the correct wgt according to ipsi
     * wgt: rho'(x) / x
     */
    switch(ipsi) {
    default:
    case 0: return(wgt_huber(x, c)); // huber
    case 1: return(wgt_biwgt(x, c)); // biweight
    case 2: return(wgt_gwgt(x, c)); // GaussWeight / "Welsh"
    case 3: return(wgt_opt(x, c)); // Optimal
    case 4: return(wgt_hmpl(x, c)); // Hampel
    case 5: return(wgt_ggw(x, c)); // GGW
    case 6: return(wgt_lqq(x, c)); // LQQ (piecewise linear psi')
    }
}

//---  Huber's rho / psi / ...
//---  -------

/* Huber's rho():  contrary to all the redescenders below,
   this can NOT be scaled to rho(Inf)=1 : */
double rho_huber(double x, const double c[])
{
    return (fabs(x) <= c[0]) ? x*x*0.5 : c[0]*(fabs(x) - c[0]/2);
}

double psi_huber(double x, const double c[])
{
// Huber's psi = rho'()
    return (x <= -c[0]) ? -c[0] : ((x < c[0]) ? x : c[0]);
}

double psip_huber(double x, const double c[])
{
// psi' = rho'' : Second derivative of Huber's loss function
    return (fabs(x) >= c[0]) ? 0. : 1.;
}

double psi2_huber(double x, const double c[])
{
// psi'' = rho''' : Third derivative of Huber's loss function
    return 0;
// FIXME? return NaN when  |x| == c ?? -- then also for psi2_hmpl()
}

double wgt_huber(double x, const double c[])
{
/*
 * Weights for Huber's loss function w(x) = psi(x)/x
 */
    return (fabs(x) >= c[0]) ? c[0]/fabs(x) : 1.;
}

//--- Biweight = Bisquare = Tukey's Biweight ...
//--- --------------------------------------

double rho_biwgt(double x, const double c[])
{
/*
 * Tukey's bisquare loss function  == R's  tukeyChi()
 */
    if (fabs(x) > (*c))
	return(1.);
    else {
	double t = x / (*c);
	t *= t; /* = t^2 */
	return( t*(3. + t*(-3. + t)) );
    }
}

double psi_biwgt(double x, const double c[])
{
/*
 * First derivative of Tukey's bisquare loss function
 */
    if (fabs(x) > (*c))
	return(0.);
    else {
	double a = x / (*c),
	    u = 1. - a*a;
	return( x * u * u );
    }
}

double psip_biwgt(double x, const double c[])
{
/*
 * Second derivative of Tukey's bisquare loss function
 */
    if (fabs(x) > (*c))
	return(0.);
    else {
	x /= *c;
	double x2 = x*x;
	return( (1. - x2) * (1 - 5 * x2));
    }
}

double psi2_biwgt(double x, const double c[])
{
/** 3rd derivative of Tukey's bisquare loss function rho()
 *= 2nd derivative of psi() :
 */
    if (fabs(x) >= c[0]) // psi''()  *is* discontinuous at x = c[0]: use "middle" value there:
	return (fabs(x) == c[0]) ? 4*x/c[0] : 0.;
    else {
	x /= c[0];
	double x2 = x*x;
	return 4*x/c[0] * (5 * x2 - 3.);
    }
}

double wgt_biwgt(double x, const double c[])
{
/*
 * Weights for Tukey's bisquare loss function
 */
    if( fabs(x) > *c )
	return(0.);
    else {
	double a = x / (*c);
	a = (1. - a)*(1. + a);
	return( a * a );
    }
}

//---------- gwgt == Gauss Weight Loss function =: "Welsh" --------------------

double rho_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Loss function
     */
    double ac = x / (*c);
    return(-expm1(-(ac*ac)/2));
}

/*========= TODO --- these limits could be lowered for the (normal) case where we have sub-normal numbers,
            ----     i.e. 0 < x < .Machine$double.xmin   min_x = 2^-1074 = 2^-(1022+52) ~= 4e-324
*/

// Largest x  such that  exp(-x) does not underflow :
static double MIN_Exp = -708.4; // ~ = M_LN2 * DBL_MIN_EXP = -log(2) * 1022 = -708.3964 */

// Largest x  such that  exp(-x^2/2) does not underflow :
static double MAX_Ex2 = 37.7; // ~ = sqrt(- 2. * M_LN2 * DBL_MIN_EXP);
/* max {x | exp(-x^2/2) < .Machine$double.xmin } =
 * min {x |  x^2 > -2*log(2)* .Machine$double.min.exp } =
 *  = sqrt(-2*log(2)* .Machine$double.min.exp) = {IEEE double}
 *  = sqrt(log(2) * 2044) = 37.64031 */

double psi_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Psi()
     */
    double a = x / (*c);
    if(fabs(a) > MAX_Ex2)
	return 0.;
    else
	return x*exp(-(a*a)/2);
}

double psip_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Psi'()
     */
    x /= (*c);
    if(fabs(x) > MAX_Ex2)
	return 0.;
    else {
	double ac = x*x;
	return exp(-ac/2) * (1. - ac);
    }
}

double wgt_gwgt(double x, const double c[])
{
    /*
     * Gauss Weight Loss function
     */
    double a = x / (*c);
    return(exp(-(a*a)/2));
}

double rho_opt(double x, const double c[])
{
    /*
     * Optimal psi Function, thank you robust package
     */
    double ac = x / (*c), // AX=S/XK
	ax = fabs(ac); // AX=ABST/XK
    if (ax > 3) // IF (AX .GT. 3.D0) THEN
	return(1); // rlRHOm2=3.25D0*XK*XK
    else if (ax > 2.) {
	const double R1 = -1.944/2., R2 = 1.728/4., R3 = -0.312/6., R4 = 0.016/8.;
	ax *= ax; // = |x/c| ^ 2
	return (ax*(R1+ ax*(R2+ ax*(R3+ ax*R4))) +1.792)/3.25;
	// rlRHOm2=XK*XK*(R1*AX**2+R2*AX**4+R3*AX**6+R4*AX**8+1.792D0)
    }
    else
	return(ac*ac/6.5); // rlRHOm2=S2/2.D0
}

double psi_opt(double x, const double c[])
{
    /*
     * Optimal psi Function, thank you robust package
     */
    double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
    double ax, ac;
    ac = x / (*c); // AX=S/XK
    ax = fabs(ac); // AX=ABST/XK
    if (ax > 3.) //    IF (AX .GT. 3.D0) THEN
	return(0); // rlPSIm2=0.D0
    else if (ax > 2.) { //  ELSE IF(AX .GT. 2.D0) THEN
	double a2 = ac*ac;
	if (ac > 0.) //     IF (AX .GT. 0.D0) THEN
	    return fmax2(0., (*c)*((((R4*a2 +R3)*a2 +R2)*a2 +R1)*ac));
	// rlPSIm2=DMAX1(0.D0,XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
	else
	    return -fabs((*c)*((((R4*a2 +R3)*a2 +R2)*a2 +R1)*ac));
	//  rlPSIm2=-DABS(XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
    }
    else
	return x;
}

double psip_opt(double x, const double c[])
{
    /*
     * psi'() for Optimal psi Function, thank you robust package
     */
    double ac = x / (*c),
	ax = fabs(ac);
    if (ax > 3.)
	return 0;
    else if (ax > 2.) {
	const double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
	ax *= ax; // = |x/c| ^ 2
	return R1 + ax*(3*R2 + ax*(5*R3 + ax * 7*R4));
    } else
	return 1;
}

double wgt_opt(double x, const double c[])
{
    /*
     * w(.) for optimal psi Function, thank you robust package
     */
    double ac = x / (*c),
	ax = fabs(ac);
    if (ax > 3.)
	return 0.;
    else if (ax > 2.) {
	const double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
	ax *= ax; // = |x/c| ^ 2
	return fmax2(0., R1+ ax*(R2 + ax*(R3 + ax*R4)));
    }
    else
	return 1.;
}


double rho_hmpl(double x, const double k[])
{
    /*
     * rho()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     *
     * This function is normalized s.t. rho(inf) = 1
     */
    double u = fabs(x),
	nc = k[0]*(k[1]+k[2]-k[0])/2;

    if (u <= k[0])
	return( x*x/2 / nc );
    else if (u <= k[1])
	return( ( u - k[0]/2 ) * k[0] / nc );
    else if (u <= k[2])
	return( ( k[1] - k[0]/2 + (u - k[1]) * (1 - ( u - k[1] ) / ( k[2] - k[1] ) / 2 )) * k[0] / nc);
    else
	return( 1 );
}

double psi_hmpl(double x, const double k[])
{
    /*
     * psi()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    // double sx = sign(x), u = fabs(x); :
    double sx, u;
    if (x < 0) { sx = -1; u = -x; } else { sx = +1; u = x; }

    if (u <= k[0])
	return( x );
    else if (u <= k[1])
	return sx * k[0];
    else if (u <= k[2])
	return sx * k[0] * (k[2] - u) / (k[2] - k[1]);
    else
	return 0.;
}

double psip_hmpl(double x, const double k[])
{
    /*
     * psi'()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    double u = fabs(x);

    if (u <= k[0])
	return( 1 );
    else if (u <= k[1])
	return( 0 );
    else if (u <= k[2])
	return( k[0] / ( k[1] - k[2] ) );
    else
	return( 0 );
}

double psi2_hmpl(double x, const double k[])
{
    /*
     * psi''()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    return 0.; // even though psi'() is already discontinuous at k[j]
}

double wgt_hmpl(double x, const double k[])
{
    /*
     * w(x) = psi(x)/x  for Hampel's redescending psi function
     * Hampel redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    double u = fabs(x);

    if (u <= k[0])
	return( 1 );
    else if (u <= k[1])
	return( k[0] / u );
    else if (u <= k[2])
	return( k[0] * ( k[2] - u ) / ( k[2] - k[1] ) / u );
    else
	return( 0 );
}


//--- GGW := Generalized Gauss-Weight    Koller and Stahel (2011)
//--- ---
// rho() & chi()  need to be calculated by numerical integration -- apart from 6 pre-stored cases
double rho_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
    if (k[0] > 0) { // for hard-coded constants --- use a *polynomial* approximation
	const double C[6][20] = { // 0: b = 1, 95% efficiency
	    {0.094164571656733, -0.168937372816728, 0.00427612218326869,
	     0.336876420549802, -0.166472338873754, 0.0436904383670537,
	     -0.00732077121233756, 0.000792550423837942, -5.08385693557726e-05,
	     1.46908724988936e-06, -0.837547853001024, 0.876392734183528,
	     -0.184600387321924, 0.0219985685280105, -0.00156403138825785,
	     6.16243137719362e-05, -7.478979895101e-07, -3.99563057938975e-08,
	     1.78125589532002e-09, -2.22317669250326e-11},
	    // 1: b = 1, 85% efficiency
	    {0.174505224068561, -0.168853188892986, 0.00579250806463694,
	     0.624193375180937, -0.419882092234336, 0.150011303015251,
	     -0.0342185249354937, 0.00504325944243195, -0.0004404209084091,
	     1.73268448820236e-05, -0.842160072154898, 1.19912623576069,
	     -0.345595777445623, 0.0566407000764478, -0.00560501531439071,
	     0.000319084704541442, -7.4279004383686e-06, -2.02063746721802e-07,
	     1.65716101809839e-08, -2.97536178313245e-10},
	    // 2: b = 1, bp 0.5
	    {1.41117142330711, -0.168853741371095, 0.0164713906344165,
	     5.04767833986545, -9.65574752971554, 9.80999125035463,
	     -6.36344090274658, 2.667031271863, -0.662324374141645,
	     0.0740982983873332, -0.84794906554363, 3.4315790970352,
	     -2.82958670601597, 1.33442885893807, -0.384812004961396,
	     0.0661359078129487, -0.00557221619221031, -5.42574872792348e-05,
	     4.92564168111658e-05, -2.80432020951381e-06},
	    // 3: b = 1.5, 95% efficiency
	    {0.104604570079252, 0.0626649856211545, -0.220058184826331,
	     0.403388189975896, -0.213020713708997, 0.102623342948069,
	     -0.0392618698058543, 0.00937878752829234, -0.00122303709506374,
	     6.70669880352453e-05, 0.632651530179424, -1.14744323908043,
	     0.981941598165897, -0.341211275272191, 0.0671272892644464,
	     -0.00826237596187364, 0.0006529134641922, -3.23468516804340e-05,
	     9.17904701930209e-07, -1.14119059405971e-08},
	    // 4: b = 1.5, 85% efficiency
	    {0.205026436642222, 0.0627464477520301, -0.308483319391091,
	     0.791480474953874, -0.585521414631968, 0.394979618040607,
	     -0.211512515412973, 0.0707208739858416, -0.0129092527174621,
	     0.000990938134086886, 0.629919019245325, -1.60049136444912,
	     1.91903069049618, -0.933285960363159, 0.256861783311473,
	     -0.0442133943831343, 0.00488402902512139, -0.000338084604725483,
	     1.33974565571893e-05, -2.32450916247553e-07},
	    // 5: b = 1.5, bp 0.5
	    {1.35010856132000, 0.0627465630782482, -0.791613168488525,
	     5.21196700244212, -9.89433796586115, 17.1277266427962,
	     -23.5364159883776, 20.1943966645350, -9.4593988142692,
	     1.86332355622445, 0.62986381140768, -4.10676399816156,
	     12.6361433997327, -15.7697199271455, 11.1373468568838,
	     -4.91933095295458, 1.39443093325178, -0.247689078940725,
	     0.0251861553415515, -0.00112130382664914}};

	double end[6] = {18.5527638190955, 13.7587939698492, 4.89447236180905,
			 11.4974874371859, 8.15075376884422, 3.17587939698492};
	int j = ((int)k[0]) - 1;
	double c;
	switch(j) {
	    // c : identical numbers to those in SET_ABC_GGW  below
	case 0: c = 1.694;     break;
	case 1: c = 1.2442567; break;
	case 2: c = 0.4375470; break;
	case 3: c = 1.063;     break;
	case 4: c = 0.7593544; break;
	case 5: c = 0.2959132; break;
	default: error(_("rho_ggw(): case (%i) not implemented."), j+1);
	}
	x = fabs(x);
	if (x <= c)
	    return(C[j][0]*x*x);
	else if (x <= 3*c)
	    return(C[j][1] +
		   x*(C[j][2] +
		      x*(C[j][3] +
			 x*(C[j][4] +
			    x*(C[j][5] +
			       x*(C[j][6] +
				  x*(C[j][7] +
				     x*(C[j][8] +
					x*(C[j][9])))))))));
	else if (x <= end[j])
	    return(C[j][10] +
		   x*(C[j][11] +
		      x*(C[j][12] +
			 x*(C[j][13] +
			    x*(C[j][14] +
			       x*(C[j][15] +
				  x*(C[j][16] +
				     x*(C[j][17] +
					x*(C[j][18]+
					   x*(C[j][19]))))))))));
	else return(1.);
    }
    else { // k[0] == 0; k[1:4] = (a, b, c, rho(Inf)) =  "general parameters"
	x = fabs(x);
	double a = 0., epsabs = R_pow(DBL_EPSILON, 0.25), result, abserr;
	int neval, ier, last, limit = 100, lenw = 4 * limit;
	int   *iwork =    (int *) R_alloc(limit, sizeof(int));
	double *work = (double *) R_alloc(lenw,  sizeof(double));

	// --> calculate integral of psi(.);  Rdqags() is from R's official API ("Writing R Extensions")
	Rdqags(psi_ggw_vec, (void *)k, &a, &x, &epsabs, &epsabs,
	       &result, &abserr, &neval, &ier,
	       &limit, &lenw, &last,
	       iwork, work);
	if (ier >= 1)
	    error(_("Error from Rdqags(psi_ggw*, k, ...): ier = %i"), ier);
	return(result/k[4]);
    }
}

void psi_ggw_vec(double *x, int n, void *k)
{
    for (int i = 0; i<n; i++) x[i] = psi_ggw(x[i], k);
}

double psi_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
#define SET_ABC_GGW(NAME)					\
    /* set a,b,c */						\
	double a, b, c;						\
	switch((int)k[0]) {					\
	    /* user specified: */				\
	case 0: a = k[1];      b = k[2]; c = k[3]; break;	\
	    /* Set of predefined cases: */			\
	case 1: a = 0.648;     b = 1.;   c = 1.694; break;	\
	case 2: a = 0.4760508; b = 1.;   c = 1.2442567; break;	\
	case 3: a = 0.1674046; b = 1.;   c = 0.4375470; break;	\
	case 4: a = 1.387;     b = 1.5;  c = 1.063; break;	\
	case 5: a = 0.8372485; b = 1.5;  c = 0.7593544; break;	\
	case 6: a = 0.2036741; b = 1.5;  c = 0.2959132; break;	\
	default: error(#NAME "_ggw: Case not implemented.");	\
	}							\
	double ax = fabs(x);

    SET_ABC_GGW(psi);
    if (ax < c) return x;
    else {
	a = -R_pow(ax-c,b)/2/a;
	return (a < MIN_Exp) ? 0. : x*exp(a);
    }
}

double psip_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
    SET_ABC_GGW(psip);
    if (ax < c) return 1.;
    else {
	double ea;
	a *= 2.;
	ea = -R_pow(ax-c,b)/a;
	return (ea < MIN_Exp) ? 0. : exp(ea) * (1 - b/a*ax*R_pow(ax-c,b-1));
    }
}

double wgt_ggw(double x, const double k[])
{
    /*
     * Gauss Weight with constant center
     */
    SET_ABC_GGW(wgt);
    return (ax < c) ? 1. : exp(-R_pow(ax-c,b)/2/a);
}
#undef SET_ABC_GGW


//--- LQQ := Linear-Quadratic-Quadratic ("lqq") --------------------------------
//--- ---    was LGW := "lin psip" := piecewise linear psi'() ------------------

// k[0:2] == (b, c, s) :
// k[0]= b = bend adjustment
// k[1]= c = cutoff of central linear part
// k[2]= s : "slope of descending": 1 - s = min_x psi'(x) =: ms

// "lin psip" := piecewise linear psi'() :
double psip_lqq (double x, const double k[])
{
    double ax = fabs(x);
    if (ax <= k[1])
	return(1.);
    else {
	double k01 = k[0] + k[1];// = b+c
	if (/*k[1] < ax && */ ax <= k01)
	    return 1. - k[2]/k[0] * (ax - k[1]);
	else {
	    double
		s5 = 1. - k[2], // = (1-s)
 		a = (k[0] * k[2] - 2 * k01)/ s5;
	    if (/* k01 < ax && */ ax < k01 + a)
		return -s5*((ax - k01)/a -1.);
	    else
		return 0.;
	}
    }
}

// piecewise linear psi'()  ==> piecewise constant  psi''():
double psi2_lqq (double x, const double k[])
{
    // double sx = sign(x), ax = fabs(x); :
    double sx, ax;
    if (x < 0) { sx = -1; ax = -x; } else { sx = +1; ax = x; }

    // k[0:2] == (b, c, s) :
    if (ax <= k[1])
	return 0.;
    else {
	double k01 = k[0] + k[1];
	if (/*k[1] < ax && */ ax <= k01)
	    return sx * (- k[2]/k[0]);
	else {
	    double
		s5 = 1. - k[2], // = (1-s)
 		a = (k[0] * k[2] - 2 * k01)/ s5;
	    if (/* k01 < ax && */ ax < k01 + a)
		return sx * (- s5 / a);
	    else
		return 0.;
	}
    }
}

double psi_lqq (double x, const double k[])
{
    double ax = fabs(x);
    if (ax <= k[1])
	return(x);
    else {
	// k[0:2] == (b, c, s) :
	double k01 = k[0] + k[1];
	if (ax <= k01)
	    return((double) (x>0 ? 1 : (x<0 ? -1 : 0)) *
		   (ax - k[2] * pow(ax - k[1], 2.) / k[0] / 2.));
	else {
	    double
		s5 = k[2] - 1., // s - 1
		s6 = -2 * k01 + k[0] * k[2]; // numerator( -a ) ==> s6/s5 = -a
	    if (/* k01 < ax && */ ax < k01 - s6 / s5)
		return((double) (x>0 ? 1 : -1) *
		       (-s6/2. - pow(s5, 2.) / s6 * (pow(ax - k01, 2.) / 2. + s6 / s5 * (ax - k01))));
	    else
		return 0.;
	}
    }
}

double rho_lqq (double x, const double k[])
{
    double ax = fabs(x), k01 = k[0] + k[1];
    if (ax <= k[1])
	return((3. * k[2] - 3.) /
	       (k[2] * k[1] * (3. * k[1] + 2. * k[0]) +
		pow(k01, 2.)) * x * x);
    else if (/* k[1] < ax && */ ax <= k01) {
	double s0 = ax - k[1];
	return((6. * k[2] - 6.) /
	       (k[2] * k[1] * (3. * k[1] + 2. * k[0]) + pow(k01, 2.)) *
	       (x * x / 2. - k[2] / k[0] * pow(s0, 3.) / 6.));
    }
    else {
	double
	    s5 = k[2] - 1.,
	    s6 = -2 * k01 + k[0] * k[2];
	if (/* k01 < ax && */ ax < k01 - s6 / s5) {
	    double s7 = ax - k01, k01_2 = pow(k01, 2.);
	    return((6. * s5) /
		   (k[2] * k[1] * (3. * k[1] + 2. * k[0]) + k01_2) *
		   (k01_2 / 2. - k[2] * k[0] * k[0] / 6. -
		    s7/2. * (s6 + s7 * (s5 + s7 * s5 * s5 / 3. / s6))));
	}
	else
	    return 1.;
    }
}

double wgt_lqq (double x, const double k[])
{
    double ax = fabs(x);
    if (ax <= k[1])
	return(1.);
    else {
	double k01 = k[0] + k[1];
	if (ax <= k01) {
	    double s0 = ax - k[1];
	    return(1. - k[2] * s0 * s0 / (2 * ax * k[0]));
	}
	else {
	    double
		s5 = k[2] - 1.,
		s6 = -2 * k01 + k[0] * k[2];
	    if (ax < k01 - s6 / s5) {
		double s7 = ax - k01;
		return(-(s6/2. + s5 * s5 / s6 * s7 * (s7/2. + s6 / s5)) / ax);
	    }
	    else
		return(0.);
	}
    }
}
/*============================================================================*/


/* This function finds the k-th place in the vector a[];
 * in the process it permutes the elements of a
 */
double kthplace(double *a, int n, int k)
{
    k--; // (0-indexing in C)

    int l=0,
	lr=n-1;
    while (l < lr) {
	double ak=a[k];
	int jnc=l,
	    j=lr;
	while (jnc <= j) {
	    while (a[jnc] < ak) jnc++;
	    while (a[j] > ak) j--;
	    if (jnc <= j) {
		double w=a[jnc]; a[jnc]=a[j]; a[j]=w;
		jnc++;
		j--;
	    }
	}
	if (j < k) l=jnc;
	if (k < jnc) lr=j;
    }
    return(a[k]);
}

/* This is from VR's bundle, MASS package  VR/MASS/src/lqs.c : */
/*
   Sampling k from 0:(n-1) without replacement:
   x[0:(k-1)] <- k distinct values from {0, .., n-1}
 */
static void sample_noreplace(int *x, int n, int k, int *ind_space)
{
#define II ind_space
    for (int i = 0; i < n; i++) II[i] = i;
    int nn=n;
    for (int i = 0; i < k; i++) { // nn == n-i
	int j = (int)(nn * unif_rand()); // in {0, .., nn-1}
	x[i] = II[j];
	II[j] = II[--nn];
    }
#undef II
}

void get_weights_rhop(const double r[], double s, int n,
		      const double rrhoc[], int ipsi, /* --> */ double *w)
{
    for(int i=0; i < n; i++) // work correctly for s == 0 (basically 0/0 = 0 in that case):
	w[i] = wgt((r[i] == 0.) ? 0. : (r[i] / s), rrhoc, ipsi);
}

/* RWLS iterations starting from i_estimate,
 * ---- the workhorse of the "lmrob_MM" algorithm, called only from R_lmrob_MM(),
 * which itself is called only from R's  lmrob..M..fit().
 * In itself,  ``just'' an M-estimator :
 */
Rboolean rwls(const double X[], const double y[], int n, int p,
	 double *estimate, double *i_estimate,
	 double *resid, double* loss,
	 double scale, double epsilon,
	 int *max_it, /* on Input:  maximal number of iterations;
			 on Output: number of iterations */
	 const double rho_c[], const int ipsi, int trace_lev)
{
    double d_beta = 0.;
    int j, iterations = 0;
    Rboolean converged = FALSE;

    double
	*wx     = (double *) R_alloc(n*p, sizeof(double)),
	*wy     = (double *) R_alloc(n,   sizeof(double)),
	*beta0  = (double *) R_alloc(p,   sizeof(double));

    int info = 1;
    INIT_WLS(wx, wy, n, p); // -> work[] etc

    COPY(i_estimate, beta0, p);
    /* calculate residuals */
    COPY(y, resid, n);
    get_Res_Xb(n,p, X, beta0, resid);

    /* main loop */
    while(!converged &&	 ++iterations < *max_it) {
	R_CheckUserInterrupt();
	/* compute weights for WLS */
	get_weights_rhop(resid, scale, n, rho_c, ipsi, weights);
	if(trace_lev >= 5) {
	    Rprintf("  it %4d: scale=%g, resid = ", iterations, scale);
	                                             disp_vec(resid,   n);
	    Rprintf("              new weights = "); disp_vec(weights, n);
	}
	/* solve weighted least squares problem */
	COPY(y, wy, n);
	FIT_WLS(X, weights, wx, wy, n, p, /* -> */ estimate);
	/* calculate residuals */
	if(trace_lev >= 5) {
	    Rprintf(" FIT_WLS() => new estimate= "); disp_vec(estimate, p);
	}
	COPY(y, resid, n);
	get_Res_Xb(n,p, X, estimate, resid);
	d_beta = norm1_diff(beta0,estimate, p);
	if(trace_lev >= 3) {
	    /* get the loss for the new estimate */
	    *loss = sum_rho_sc(resid,scale,n,0,rho_c,ipsi);
	    Rprintf("  it %4d: L(b1) = %#.12g ", iterations, *loss);
	    if(trace_lev >= 4) {
		Rprintf("\n  b1 = (");
		for(j=0; j < p; j++)
		    Rprintf("%s%.11g", (j > 0)? ", " : "", estimate[j]);
		Rprintf(");");
	    }
	    Rprintf(" ||b0 - b1||_1 = %g\n", d_beta);
	}
	/* check for convergence */
	converged = d_beta <= epsilon * fmax2(epsilon, norm1(estimate, p));
	COPY(estimate, beta0, p);
    } /* end while(!converged & iter <=...) */

    if(trace_lev > 0) {
	if(trace_lev < 3) *loss = sum_rho_sc(resid,scale,n,0,rho_c,ipsi);
	Rprintf(" rwls() used %2d it.; last ||b0 - b1||_1 = %#g, L(b1) = %.12g; %sconvergence\n",
		iterations, d_beta, *loss, (converged ? "" : "NON-"));
    }

    *max_it = iterations;

    return converged;
} /* rwls() */

/* sets the entries of a matrix to zero */
void zero_mat(double **a, int n, int m)
{
    int i,j;
    for(i=0; i < n; i++)
	for(j=0; j < m; j++)
	    a[i][j] = 0.;
}

/*
 *
 * 2004 / 5 -- Matias Salibian-Barrera & Victor Yohai
 * Department of Statistics, University of British Columbia
 * matias@stat.ubc.ca
 * Department of Mathematics, University of Buenos Aires
 * vyohai@uolsinectis.com.ar
 *
 *
 * Reference: A fast algorithm for S-regression estimates,
 * 2005, Salibian-Barrera and Yohai.
 */

/* This function implements the "large n" strategy
 */
void fast_s_large_n(double *X, double *y, double s_y,
		    int n, int p, int nRes, int *max_it_scale, double *res,
		    int groups, int n_group,
		    int K, int *max_k,
		    double rel_tol, double inv_tol, double scale_tol, double zero_tol,
		    int *converged,
		    int best_r, double bb, const double rrhoc[], int ipsi,
		    double *bbeta, double *sscale,
		    int trace_lev, int mts, Rboolean ss)
{
/* *X  = the n x p  design matrix (incl. intercept if appropriate),
 *	 in column order as sent by R)
 * *y  = the ( n ) response vector
 * n = the length of y
 * p = the number of columns in X
 * nRes  = number of re-sampling candidates to be used in each partition

 * groups = number of groups in which to split the
 *	      random subsample
 * n_group = size of each of the (groups) groups
 *	       to use in the random subsample
 * K      = number of refining steps for each candidate (typically 1 or 2)
 * *max_k = [on Input:] number of refining steps for each candidate (typically 1 or 2,
 *                      used to be hard coded to MAX_ITER_REFINE_S = 50 )
 *	    [on Output:] effectively used iterations
 * *rel_tol= convergence tolerance for iterative refinement iterations
             [used to be hard coded to EPS = 1e-7 ]
 * *converged: will become 0(FALSE)  iff at least one of the best_r iterations
 *             did not converge (in max_k steps to rel_tol precision)
 * best_r = no. of best candidates to be iterated further ("refined")
 * bb	   = right-hand side of S-equation (typically 1/2)
 * *rrhoc  = tuning constant(s) for loss function
 *	     (this should be associated with bb)
 * ipsi    = indicator for type of psi function to be used
 * *bbeta  = final estimator 'best beta'
 * *sscale = associated scale estimator (or -1 when problem)
 */
    double sc, best_sc, worst_sc,
	INF_sc = s_y * INFI; // may overflow to 'Inf'

    int sg = groups * n_group, // (integer overflow already checked in R code)
	k  = groups * best_r;
    double
	*beta_ref          = (double *) R_alloc(p, sizeof(double)),
	*final_best_scales = (double *) R_alloc(best_r, sizeof(double)),
	      *best_scales = (double *) R_alloc(   k,	sizeof(double)),
	*xsamp =     (double *) R_alloc(n_group*p, sizeof(double)),
	*ysamp =     (double *) R_alloc(n_group,   sizeof(double));
    int *indices =   (int *)    R_alloc(sg,  sizeof(int)),
	*ind_space = (int *)    R_alloc(n,   sizeof(int));

#define CALLOC_MAT(_M_, _n_, _d_)			\
    _M_ = (double **) R_Calloc(_n_, double *);		\
    for(int i=0; i < _n_; i++)				\
	_M_[i] = (double *) R_alloc(_d_, sizeof(double))

    double **final_best_betas; CALLOC_MAT(final_best_betas, best_r, p);
    double **      best_betas; CALLOC_MAT(      best_betas,    k  , p);

    /* assume that n > 2000 */

    /*	set the seed */
    GetRNGstate();

    /* get a sample of k indices */
    sample_noreplace(indices, n, sg, ind_space);
    /* FIXME: define groups using nonsingular subsampling? */
    /*        would also need to allow observations to be part */
    /*        of multiple groups at the same time */
    /* FIXME: Also look at lqs_setup(),
     * -----  and  xr[.,.] "fortran-like" matrix can be used from there!*/

/* For each (of 'groups') group : get the best_r best betas : */

#define X(_k_, _j_) X[(_j_)*n + _k_]
#define xsamp(_k_, _j_) xsamp[(_j_)*n_group + _k_]

    int S_code;
    for(int i=0; i < groups; i++) {
	/* populate "matrix" for group [i] */
	int ij = i*n_group;
	for(int j = 0; j < n_group; j++, ij++) { // ij == i*n_group + j
	    for (k = 0; k < p; k++)
		xsamp(j, k) = X(indices[ij], k);
	    ysamp[j] = y[indices[ij]];
	}
	if (trace_lev)
	    Rprintf(" Subsampling to find candidate betas in group %d:\n", i);
	S_code = fast_s_with_memory(xsamp, ysamp, s_y, res,
				    n_group, p, nRes, max_it_scale, K, max_k,
				    rel_tol, inv_tol, scale_tol, zero_tol,
				    trace_lev, best_r, bb, rrhoc, ipsi,
				    mts, ss, /* --> */
				    best_betas + i* best_r,
				    best_scales+ i* best_r);
	switch(S_code) {
	case 0: /* normal case: --> go on */ break;
	case 1: /* singular ! */
	    *sscale = -1.; /* problem */
	    if(trace_lev)
		Rprintf(" 'singularity' from fast_s_with_memory() in group %d\n", i+1);
	    goto cleanup_and_return;
	case 2: /* scale = 0 <==> exact fit in group [i] */ ;
	    *sscale = 0.; /* and the corresponding "best beta": */
	    COPY(best_betas[i], bbeta, p);
	    if(trace_lev) {
		Rprintf(" scale=0 from fast_s_with_memory() in group %d\n", i+1);
		if(trace_lev >= 3) { Rprintf(" and bbeta[]: "); disp_vec(bbeta, p); }
	    }
	    goto cleanup_and_return;
	default: error(_("fast_s_with_memory() returned invalid code %d"), S_code);
	}
    }

#undef xsamp

/* now	iterate (refine) these "best_r * groups"
 * best betas in the (xsamp,ysamp) sample
 * with K  C-steps and keep only the "best_r" best ones
 */
    /* initialize new work matrices */
    double
     *w_beta= (double *) R_alloc(p,   sizeof(double)),
  	*wx = (double *) R_alloc(n*p, sizeof(double)), // need only k here,
	*wy = (double *) R_alloc(n,   sizeof(double)); // but n in the last step
    xsamp =   (double *) R_alloc(sg*p, sizeof(double));
    ysamp =   (double *) R_alloc(sg,   sizeof(double));

#define xsamp(_k_,_j_) xsamp[(_j_)*sg + _k_]

    for (int ij = 0; ij < sg; ij++) {
	for (k = 0; k < p; k++)
	    xsamp(ij, k) = X(indices[ij],k);
	ysamp[ij] = y[indices[ij]];
    }

    int info = 1;
    INIT_WLS(wx, wy, n, p);
    Rboolean conv = FALSE;
    int pos_worst_scale = 0;
    for(int i=0; i < best_r; i++)
	final_best_scales[i] = INF_sc;
    worst_sc = INF_sc;
    /* set the matrix to zero */
    zero_mat(final_best_betas, best_r, p);
    for(int i=0; i < (best_r * groups); i++) {
	if(trace_lev >= 3) {
	    Rprintf("  Sample[%3d]: before refine_(*, conv=FALSE):\n", i);
	    if(i > 0) {
		Rprintf("   beta_cand : "); disp_vec(best_betas[i],p);
		Rprintf("   with scale %.15g\n", best_scales[i]);
	    }
	}
	int ik = refine_fast_s(xsamp, wx, ysamp, s_y, wy, weights, sg, p, res,
			   work, lwork, best_betas[i], w_beta,
			   &conv/* = FALSE*/, K, rel_tol, zero_tol, trace_lev,
			   bb, rrhoc, ipsi, best_scales[i], /* -> */ beta_ref, &sc);
	if(trace_lev >= 3) {
	    Rprintf("   after refine: beta_ref : "); disp_vec(beta_ref,p);
	    Rprintf("   with scale %.15g\n", sc);
	    if(ik < 0) Rprintf("* exact fit! %d zero residuals", -ik);
	}
	if ( sum_rho_sc(res, worst_sc, sg, p, rrhoc, ipsi) < bb ) {
	    int scale_iter = *max_it_scale;
	    /* scale will be better */
	    sc = find_scale(res, bb, rrhoc, ipsi, sc, sg, p, &scale_iter, scale_tol, trace_lev >= 3);
	    int k2 = pos_worst_scale;
	    final_best_scales[ k2 ] = sc;
	    COPY(beta_ref, final_best_betas[k2], p);
	    pos_worst_scale = find_max(final_best_scales, best_r);
	    worst_sc = final_best_scales[pos_worst_scale];
	}
    }

/* now iterate the best "best_r"
 * betas in the whole sample until convergence (max_k, rel_tol)
 */
    best_sc = INF_sc; *converged = 1;  k = 0;
    if(trace_lev)
	Rprintf(" Now refine() to convergence for %d very best ones:\n",
		best_r);

    for(int i=0; i < best_r; i++) {
	conv = TRUE;
	if(trace_lev >= 4) Rprintf("  i=%d:\n", i);
	int it_k = refine_fast_s(X, wx, y, s_y, wy, weights, n, p, res,
				work, lwork, final_best_betas[i], w_beta,
			        &conv/* = TRUE */, *max_k, rel_tol, zero_tol, trace_lev,
			        bb, rrhoc, ipsi, final_best_scales[i], /* -> */ beta_ref, &sc);
	if(trace_lev) {
	    Rprintf("  Final best[%d]: %sconvergence ", i, conv ? "" : "NON ");
	    if(it_k >= 0)
		Rprintf("(%d iter.)", it_k);
	    else
		Rprintf("(Exact fit! %d zeroes; scale=0, sc=%g <?< best_sc=%g)",
			-it_k, sc, best_sc);
	}
	if(best_sc > sc || sc == 0) {
	    if(trace_lev && best_sc > sc)
		Rprintf(": -> improved scale to %.15g", sc);
	    best_sc = sc;
	    COPY(beta_ref, bbeta, p);
	    if(best_sc == 0.) {
		if(trace_lev) Rprintf(" = 0  ==> finish\n");
		break;
	    }
	}
	if (trace_lev) Rprintf("\n");
	if (!conv && *converged) *converged = 0;
	if (it_k >= 0 && k < it_k) k = it_k;
    }
    *sscale = best_sc;
    *max_k = k;

/* Done. Now clean-up. */

  cleanup_and_return:
    PutRNGstate();

    R_Free(best_betas);
    R_Free(final_best_betas);

#undef X
#undef xsamp

} /* fast_s_large_n() */

// (called only in 1 place above)
int fast_s_with_memory(double *X, double *y, double s_y, double *res,
		       int n, int p, int nResample, int *max_it_scale,
		       int K, int *max_k,
		       double rel_tol, double inv_tol, double scale_tol, double zero_tol,
		       int trace_lev, int best_r, double bb, const double rrhoc[],
		       int ipsi, int mts, Rboolean ss,
		       // ==> result
		       double **best_betas, double *best_scales)
{
/*
 * Called from >>> fast_s_large_n() <<<, the adjustment for large "n",
 * same as fast_s, but it returns the best_r best betas,
 * and their associated scales --- unless there was an exact fit (scale = 0) case in one group.
 *
 * x       : an n x p design matrix (including intercept if appropriate)
 * y       : an n vector
 * res     : an n vector of residuals
 * nRes   = number of re-sampling candidates to be taken
 * K	   = number of refining steps for each candidate
 * best_r = number of (refined) to be retained for full iteration
 * bb	   = right-hand side of the S-equation (typically 1/2)
 * rrhoc[] = tuning constant(s) for loss function
 *	     (this should be associated with bb)
 * ipsi  = indicator for type of loss function to be used
 * RETURN
 * *best_betas	=  the best ... coefficient vectors
 * *best_scales =  their associated residual scales
 *
 * return(sing) : TRUE if have singularity
 * */

    Rboolean conv = FALSE,
	sing = FALSE; // sing = TRUE|FALSE  the final result
    double sc,
	INF_sc = s_y * INFI, /* may overflow to Inf */
	worst_sc = INF_sc;

    int info = 1; // is set by *both* :
    SETUP_SUBSAMPLE(n, p, X, 1);
    INIT_WLS(X, y, n, p);

    double
	*wx =        (double *) R_alloc(n*p, sizeof(double)),
	*wy =        (double *) R_alloc(n,   sizeof(double)),
	*w_beta    = (double *) R_alloc(p,   sizeof(double)),
	*beta_cand = (double *) R_alloc(p,   sizeof(double)),
	*beta_ref  = (double *) R_alloc(p,   sizeof(double));

    for(int i=0; i < best_r; i++)
	best_scales[i] = INF_sc;
    int pos_worst_scale = 0;

/* resampling approximation  */

    for(int i=0; i < nResample; i++) {
	R_CheckUserInterrupt();
	/* find a candidate */
	sing = (Rboolean) // 0 |--> FALSE (= success);  {1,2} |-> TRUE
	    subsample(Xe, y, n, p, beta_cand, ind_space, idc, idr, lu, v, pivot,
		      Dr, Dc, rowequ, colequ, 1, mts, ss, inv_tol, 1);
	if (sing) {
	    for (int k=0; k < best_r; k++) best_scales[i] = -1.;
	    return sing;
	}
	/* FIXME: is_ok ?? */

	/* improve the re-sampling candidate */

	/* conv = FALSE : do K refining steps */
	int ik = refine_fast_s(X, wx, y, s_y, wy, weights, n, p, res,
			       work, lwork, beta_cand, w_beta, &conv/* = FALSE*/, K,
			       rel_tol, zero_tol, trace_lev, bb, rrhoc, ipsi, -1.,
			       /* -> */ beta_ref, &sc);
	if(ik < 0) {
	    if(trace_lev) Rprintf(" * exact fit! %d zero residuals; scale = 0\n", -ik);

/// YES: FIXME ---> return beta_cand and be done - as in fast_s(..)
	}

	if ( sum_rho_sc(res, worst_sc, n, p, rrhoc, ipsi) < bb )	{
	    int scale_iter = *max_it_scale;
	    /* scale will be better */
	    sc = find_scale(res, bb, rrhoc, ipsi, sc, n, p,
			    &scale_iter, scale_tol, trace_lev >= 3);
	    int k = pos_worst_scale;
	    best_scales[ k ] = sc;
	    for(int j=0; j < p; j++)
		best_betas[k][j] = beta_ref[j];
	    pos_worst_scale = find_max(best_scales, best_r);
	    worst_sc = best_scales[pos_worst_scale];
	    if (trace_lev >= 2) {
	      Rprintf("  Sample[%3d]: found new candidate with scale %.7g in %d iter (worst sc=%.5g)\n",
		      i, sc, scale_iter, worst_sc);
	    }
	}
    } /* for(i ) */

    return sing;
} /* fast_s_with_memory() */

void fast_s(double *X, double *y, double s_y,
	    int n, int p, int nResample, int *max_it_scale, double *res,
	    int K, int *max_k,
	    double rel_tol, double inv_tol, double scale_tol, double zero_tol,
	    int *converged,
	    int best_r, double bb, const double rrhoc[], int ipsi,
	    double *bbeta, double *sscale, int trace_lev, int mts, Rboolean ss)
{
/* *X  = the n x p  design matrix (incl. intercept if appropriate),
 *	 in column order as sent by R)
 * *y  = the ( n ) response vector
 * nn =: n = the length of y
 * pp =: p = the number of columns in X
 * nRes   = number of re-sampling candidates to be taken
 * K	   = number of refining steps for each candidate
 * best_r = number of (refined) to be retained for full iteration
 * *converged: will become FALSE  iff at least one of the best_r iterations
 *	       did not converge (in max_k steps to rel_tol precision)
 * bb	   = right-hand side of the S-equation (typically 1/2)
 * *rrhoc  = tuning constant(s) for loss function
 *	     (this should be associated with bb)
 * iipsi  = indicator for type of loss function to be used
 * *bbeta  = "best beta" = final estimator
 * *sscale = associated scale estimator (or -1 when problem)
 */
    double sc, best_sc, aux;
    double INF_sc = s_y * INFI; // may overflow to 'Inf'

    if(trace_lev >= 2) Rprintf("fast_s(*, s_y=%g, n=%d, p=%d, ipsi=%d, ..)", s_y, n, p, ipsi);
    int info;
    SETUP_SUBSAMPLE(n, p, X, 0);

    // More arrays, allocated:
    double
	*wx = (double *) R_alloc(n*p, sizeof(double)),
	*wy = (double *) R_alloc(n,   sizeof(double)),
	*w_beta      = (double *) R_alloc(p, sizeof(double)),
	*beta_cand   = (double *) R_alloc(p, sizeof(double)),
	*beta_ref    = (double *) R_alloc(p, sizeof(double)),
	*best_scales = (double *) R_alloc(best_r, sizeof(double)),
	// matrix:
	**best_betas = (double **) R_Calloc(best_r, double *);
    for(int i=0; i < best_r; i++) {
	best_betas[i] = (double*) R_alloc(p, sizeof(double));
	best_scales[i] = INF_sc;
    }

    if(trace_lev >= 2) Rprintf(" before INIT_WLS():\n");
    INIT_WLS(wx, wy, n, p); // (re-setting 'info')

    /* disp_mat(x, n, p); */

    int pos_worst_scale = 0;
    Rboolean conv = FALSE;
    double worst_sc = INF_sc;
    /* srand((long)*seed_rand); */
    GetRNGstate();

/* resampling approximation  */
    if (trace_lev)
	Rprintf(" Subsampling %d times to find candidate betas:\n", nResample);

    for(int i=0; i < nResample; i++) {

	R_CheckUserInterrupt();
	/* find a candidate */
	Rboolean sing = (Rboolean) // 0 |--> FALSE (= success);  {1,2} |-> TRUE
	    subsample(Xe, y, n, p, beta_cand, ind_space, idc, idr, lu, v, pivot,
		      Dr, Dc, rowequ, colequ, 1, mts, ss, inv_tol, 1);
	if (sing) {
	    *sscale = -1.;
	    if (trace_lev >= 1)
		Rprintf("  Sample[%3d]: singular subsample() - giving up!\n", i);
	    goto cleanup_and_return;
	}
	Rboolean trace_sample = trace_lev >= 3 && (p <= 9 || trace_lev >= 5);
	if (trace_sample) {
	    Rprintf("  Sample[%3d]: idc = ", i); disp_veci(idc, p);
	    if(p <= 3 || trace_lev >= 5) {
		Rprintf("   b^[] = "); disp_vec(beta_cand, p);
	    }
	}

	/* improve the re-sampling candidate */

	/* conv = FALSE : do K refining steps */
	int ik = refine_fast_s(X, wx, y, s_y, wy, weights, n, p, res,
		      work, lwork, beta_cand, w_beta, &conv/* = FALSE*/, K,
		      rel_tol, zero_tol, trace_lev, bb, rrhoc, ipsi, -1.,
		      /* -> */ beta_ref, &sc);
	if(trace_lev >= 3) {
	    double del = norm_diff(beta_cand, beta_ref, p);
	    if(!trace_sample) Rprintf("  Sample[%3d]:", i);
	    Rprintf("   after refine_(*, conv=F):\n"
		    "   beta_ref : "); disp_vec(beta_ref,p);
	    Rprintf("   with ||beta_ref - beta_cand|| = %.12g, --> sc = %.15g\n",
		    del, sc);
	}
	if(fabs(sc) == 0.) { /* exact zero set by refine_*() */
	    if(trace_lev >= 1)
		Rprintf(" |s|=0: Have %d (too many) exact zeroes -> leaving refinement!\n", -ik);
	    *sscale = sc;
	    /* *converged = 1; -- important when used as init in lmrob() --
	     * ---------------  but for that need at least valid residuals, fitted,
	     *  and possibly a *scale* sc > 0 ??? -------------- FIXME
	     */
	    COPY(beta_cand, bbeta, p);
	    goto cleanup_and_return;
	}
	if ( sum_rho_sc(res, worst_sc, n, p, rrhoc, ipsi) < bb )	{
	    int scale_iter = *max_it_scale;
	    /* scale will be better */
	    sc = find_scale(res, bb, rrhoc, ipsi, sc, n, p,
			    &scale_iter, scale_tol, trace_lev >= 3);
	    int k = pos_worst_scale;
	    best_scales[ k ] = sc;
	    COPY(beta_ref, best_betas[k], p);
	    pos_worst_scale = find_max(best_scales, best_r);
	    worst_sc = best_scales[pos_worst_scale];
	    if (trace_lev >= 2) {
		if(trace_lev < 3) /* not yet "Sample[..]" */ Rprintf("  Sample[%3d]:", i);
		Rprintf("   found new candidate with scale %.7g in %d iter (worst sc=%.5g)\n",
			sc, scale_iter, worst_sc);
	    }
	}
    } /* for(i in 1..nResample) */

/* now look for the very best */
    if(trace_lev)
	Rprintf(" Now refine() to convergence for %d very best ones:\n", best_r);

    best_sc = INF_sc; *converged = 1;
    int k = 0;
    for(int i=0; i < best_r; i++) {
	conv = TRUE;
	if(trace_lev >= 4) Rprintf("  i=%d:\n", i);
	int it_k = refine_fast_s(X, wx, y, s_y, wy, weights, n, p, res, work, lwork,
				best_betas[i], w_beta, &conv /* = TRUE */, *max_k,
				rel_tol, zero_tol, trace_lev, bb, rrhoc, ipsi,
				best_scales[i], /* -> */ beta_ref, &aux);
	if(trace_lev) {
	    Rprintf("  Best[%d]: %sconvergence ", i, conv ? "" : "NON ");
	    if(it_k >= 0)
		Rprintf("(%d iter.)", it_k);
	    else
		Rprintf("(%d zeroes; scale=0)", -it_k);
	}
	if(aux < best_sc) {
	    best_sc = aux;
	    COPY(beta_ref, bbeta, p);
	    if(trace_lev) {
		Rprintf(": -> improved scale to %.15g", best_sc);
		if(trace_lev >= 2) { Rprintf("; bbeta= "), disp_vec(bbeta,p); }
	    }
	}
	if(trace_lev) Rprintf("\n");
	if (!conv && *converged) *converged = 0; // will stay at 0, even when next best_betas[i] is good
	if (k < it_k) k = it_k;
    }
    *sscale = best_sc;
    *max_k = k;

  cleanup_and_return:

    PutRNGstate();

    R_Free(best_betas);

    return;
} /* fast_s() */

int refine_fast_s(const double X[], double *wx,
		  const double y[], double s_y, double *wy, double *weights,
		  int n, int p, double *res,
		  double *work, int lwork, const double beta_cand[], double *beta_j,
		  Rboolean *conv, int kk, double rel_tol, double zero_tol,
		  int trace_lev, // here: only print when trace_lev >= 3
		  double b, const double rrhoc[], int ipsi, double initial_scale,
		  double *beta_ref, double *scale)
{
/*
 * X	   = matrix (n x p) of explanatory variables
 * y	   = vector ( n )   of responses
 * weights = robustness weights wt[] * y[] (length n) Output
 * res	   = residuals	y[] - x[,] * beta  (length n) Output
 * conv: FALSE means do  kk  refining steps     (and conv stays FALSE)
 *	 TRUE  means refine until convergence(rel_tol, max{k} = kk)
 *             and in this case, 'conv' *returns* TRUE if refinements converged
 * rel_tol  =
 * zero_tol = small positive, determining "zero residuals" - was hardwired to EPS_SCALE := 1e-10
 * beta_cand= candidate beta[] (of length p)	Input
 * beta_j   = internal "working" beta[] (of length p)
 * is	    = initial scale			input

 * beta_ref = resulting beta[] (of length p)	Output
 * scale    = final scale			Output

 * for FIT_WLS, DGELS:
 * wx       = matrix (n x p)
 * wy       = vector of length n  {also used for 'beta' (<==> p <= n )}
 * work     = vector of length lwork
 * lwork    = length of vector work

 * return( <number of iteration steps> *or* -z,  where z = zeroes := #{i; |R~_i| <= e_z};
                                                                     e_z := zero_tol
 */
    Rboolean trace_beta, converged = FALSE;/* Wall */
    if (trace_lev >= 3) {
	Rprintf("   refine_fast_s(s0=%g, convChk=%s): ", initial_scale, *conv ? "TRUE" : "FALSE");
	trace_beta = (p <= 6 || trace_lev >= 5);
	if(trace_beta) {
	    Rprintf("beta_cand= "), disp_vec(beta_cand,p);
	}
    }

    /* calculate residuals */
    COPY(y, res, n);
    get_Res_Xb(n,p, X, beta_cand, /*--> */ res);

    if( initial_scale < 0. )
	initial_scale = MAD(res, n, 0., wy, weights);// wy and weights used as work arrays

    double s0 = s_y * zero_tol;
    int zeroes = 0;
    for(int j=0; j < n; j++) { // ensuring "eps_zero" to be equivariant in y[] :
	if(fabs(res[j]) <= s0)
	    zeroes++;
    }
    if (trace_lev >= 4)
	Rprintf("   |{i; |R_i| <= %.4g ~= 0}| = %d zeroes (zero_tol=%.3g, s_y=%g);\n",
		s0, zeroes, zero_tol, s_y);

/* if "perfect fit", return it with a 0 assoc. scale */
    if(initial_scale <= 0. || zeroes > ((double)n)/2.) { /* <<- FIXME: should depend on 'b' ! */
	// if(zeroes > (((double)n + (double)p)/2.)) -- corresponding to MAD_p but have MAD=MAD_0 above
	COPY(beta_cand, beta_ref, p);
	if (trace_lev >= 3)
	    Rprintf("   too many zeroes -> scale=0 & quit refinement\n");
	*scale = 0.;
	return -zeroes; // (for diagnosis + emphasize special scale=0)
    }

    s0 = initial_scale; // > 0
    int info, iS;
    if(trace_lev >= 4)
	Rprintf("   %s %d refinement iterations, starting with s0=%g\n",
		*conv ? "maximally" : "doing", kk, s0);
    COPY(beta_cand, beta_j, p);
    for(iS=0; iS < kk; iS++) {
	/* one step for the scale */
	s0 = s0 * sqrt( sum_rho_sc(res, s0, n, p, rrhoc, ipsi) / b );
	/* compute weights for WLS */
	get_weights_rhop(res, s0, n, rrhoc, ipsi, weights);
        /* solve weighted least squares problem */
	COPY(y, wy, n); // wy = y[1:n] on in input, and beta[1:p] on output
	FIT_WLS(X, weights, wx, wy, n, p, /* -> */ beta_ref);
	if(*conv) { /* check for convergence */
	    double del = norm_diff(beta_j, beta_ref, p);
	    double nrmB= norm(beta_j, p);
	    converged = (del <= rel_tol * fmax2(rel_tol, nrmB));
	    if(trace_lev >= 4)
		Rprintf("     it %4d, ||b[i]||= %#.12g, ||b[i] - b[i-1]||= %#.15g --> conv=%s\n",
			iS, nrmB, del, (converged ? "TRUE" : "FALSE"));
	    if(converged)
		break;
	}
	/* calculate residuals */
	COPY(y, res, n);
	get_Res_Xb(n,p, X, beta_ref, /*--> */ res);
	COPY(beta_ref, beta_j, p);
    } /* for(iS = 0; iS < kk ) */

    if(*conv) {
	if(converged) {
	    if(trace_lev >= 3) Rprintf("refine_() converged after %d iterations\n", iS);
	}
	else { // !converged
	    *conv = FALSE;
	    warning(_("S refinements did not converge (to refine.tol=%g) in %d (= k.max) steps"),
		    rel_tol, iS);
	}
    }
    *scale = s0;
    return iS; /* number of refinement steps */
} /* refine_fast_s() */


/* Subsampling part for M-S algorithm,  i.e. called only from  R_lmrob_M_S()
 * Recreates RLFRSTML function found in robust/src/lmrobml.f of package 'robust' */
void m_s_subsample(double *X1, double *y, int n, int p1, int p2,
		   int nResample, int max_it_scale,
		   double rel_tol, double inv_tol, double scale_tol,
		   double zero_tol,
		   double bb,
		   const double rrhoc[], int ipsi, double *sscale, int trace_lev,
		   double *b1, double *b2, double *t1, double *t2,
		   double *y_tilde, double *res, double *x1, double *x2,
		   int *NIT, int *K, int *KODE, double *SIGMA, double BET0,
		   double *SC1, double *SC2, double *SC3, double *SC4,
		   int mts, Rboolean ss)
{
    int i, p = p1 + p2, info;
    double sc = /* INF_sc := */ *SIGMA * INFI;
    if(sc > DBL_MAX) sc = DBL_MAX; // cannot use 'Inf' here
    *sscale = sc;

    if (trace_lev >= 2)
	Rprintf(" Starting M-S subsampling procedure(p1=%d, p2=%d; ini.sc=%g) .. ", p1,p2, sc);

    SETUP_SUBSAMPLE(n, p2, x2, 0);

    /*	set the seed */
    GetRNGstate();

    if (trace_lev >= 2) Rprintf(" [setup Ok]\n");

    for(i=0; i < nResample; i++) {
	R_CheckUserInterrupt();
	/* STEP 1: Draw a subsample of size p2 from (X2, y) */
	Rboolean sing = (Rboolean) // 0 |--> FALSE (= success);  {1,2} |-> TRUE
	    subsample(Xe, y, n, p2, t2, ind_space, idc, idr, lu, v, pivot,
		      Dr, Dc, rowequ, colequ, /* sample= */ TRUE, mts,
		      ss, inv_tol, /*solve = */ TRUE);
	if (sing) {
	    *sscale = -1.;
	    goto cleanup_and_return;
	}
	/* calculate partial residuals */
	COPY(y, y_tilde, n);
        get_Res_Xb(n,p2, x2, t2, y_tilde);
	/* STEP 3: Obtain L1-estimate of b1 */
	COPY(X1, x1, n*p1);
	F77_CALL(rllarsbi)(x1, y_tilde, &n, &p1, &n, &n, &rel_tol,
			   NIT, K, KODE, SIGMA, t1, res, SC1, SC2,
			   SC3, SC4, &BET0);
	if (*KODE > 1) { // KODE in {0, 1} is ok
	    REprintf("m_s_subsample(): Problem in RLLARSBI (L1-regr). KODE=%d. Exiting.\n",
		     *KODE);
	    *sscale = -1.;
	    goto cleanup_and_return;
	}
	/* STEP 4: Check if candidate looks promising */
	if (sum_rho_sc(res, *sscale, n, p, rrhoc, ipsi) < bb) {
	    int scale_iter = max_it_scale;
	    /* scale will be better */
	    /* STEP 5: Solve for sc */
	    sc = find_scale(res, bb, rrhoc, ipsi, sc, n, p, &scale_iter,
			    scale_tol, trace_lev >= 4);
	    if(trace_lev >= 2)
		Rprintf("  Sample[%3d]: new candidate with sc =%#16.9g in %d iter\n",
			i, sc, scale_iter);
	    /* STEP 6: Update best fit */
	    *sscale = sc;
	    COPY(t1, b1, p1);
	    COPY(t2, b2, p2);
	    if (sc < zero_tol) {
		REprintf("\nScale too small\n"
			 "Aborting m_s_subsample()\n\n");
		*sscale = -1.;
		goto cleanup_and_return;
	    }
	}
    } /* for(i ) */

    /* STEP 7: Clean up and return */
    if (trace_lev >= 1) {
	Rprintf(" Finished M-S subsampling with scale = %.5f\n",*sscale);
#define maybe_SHOW_b1_b2			\
	if (trace_lev >= 3) {			\
	     Rprintf("  b1: "); disp_vec(b1,p1);\
	     Rprintf("  b2: "); disp_vec(b2,p2);\
	}
	maybe_SHOW_b1_b2;
    }

  cleanup_and_return:
    PutRNGstate();
} /* m_s_subsample() */

/* Descent step for M-S algorithm
 * Return value: convergence; note that convergence is *not* guaranteed
 */
Rboolean m_s_descent(double *X1, double *X2, double *y,
		 int n, int p1, int p2, int K_m_s, int max_k, int max_it_scale,
		 double rel_tol, double scale_tol, double bb, const double rrhoc[],  int ipsi,
		 double *sscale, int trace_lev,
		 double *b1, double *b2, double *t1, double *t2,
		 double *y_tilde, double *res, double *res2, double *x1, double *x2,
		 int *NIT, int *K, int *KODE, double *SIGMA,  double BET0,
		 double *SC1, double *SC2, double *SC3, double *SC4)
{
    int nnoimpr = 0, nref = 0;
    int p = p1 + p2;
    Rboolean converged = FALSE;
    double sc = *sscale;

    COPY(b1, t1, p1);
    COPY(b2, t2, p2);
    COPY(res, res2, n);

    if (trace_lev >= 2)
	Rprintf(" Starting descent procedure...\n");

    int info = 1;
    INIT_WLS(x2, y, n, p2);

    if (trace_lev >= 3) {
	Rprintf("  Scale: %.5f\n", *sscale);
	if (trace_lev >= 5) {
	    Rprintf("   res2: "); disp_vec(res2,n);
	}
    }

    /* Do descent steps until there is no improvement for   */
    /* K_m_s steps or we are converged                      */
    /* (convergence is not guaranteed)                      */
    while ( (nref++ < max_k) & (!converged) & (nnoimpr < K_m_s) ) {
	R_CheckUserInterrupt();
	/* STEP 1: update b2 (save it to t2) */
	/* y_tilde = y - x1 %*% t1 */
	COPY(y, y_tilde, n);
	COPY(X1, x1, n*p1); // FIXME(MM)?  don't need x1, use X1 (which is unchanged!)
	get_Res_Xb(n,p1, x1, t1, y_tilde);
	/* compute weights for WLS */
	get_weights_rhop(res2, sc, n, rrhoc, ipsi, weights);
	/* solve weighted least squares problem */
	FIT_WLS(X2, weights, x2, y_tilde, n, p2,  /* -> */ t2);
        /* get (intermediate) residuals */
	COPY(y, res2, n);
	get_Res_Xb(n,p2, X2, t2, res2);
	/* STEP 2: Obtain L1-estimate of b1 */
	COPY(res2, y_tilde, n);
	F77_CALL(rllarsbi)(x1, y_tilde, &n, &p1, &n, &n, &rel_tol,
			   NIT, K, KODE, SIGMA, t1, res2,
			   SC1, SC2, SC3, SC4, &BET0);
	if (*KODE > 1) {
	    error(_("m_s_descent(): Problem in RLLARSBI (RILARS). KODE=%d. Exiting."),
		  *KODE);
	}
	/* STEP 3: Compute the scale estimate */
	int scale_iter = max_it_scale;
	sc = find_scale(res2, bb, rrhoc, ipsi, sc, n, p, &scale_iter,
			scale_tol, trace_lev >= 4); // <- here only if higher trace_lev
	/* STEP 4: Check for convergence */
	/* FIXME: check convergence using scale ? */
	double del = sqrt(norm_diff2(b1, t1, p1) + norm_diff2(b2, t2, p2));
	double nrmB = sqrt(norm2(t1, p1) + norm2(t2, p2));
	converged = (del < rel_tol * fmax2(rel_tol, nrmB));
	if (trace_lev >= 3) {
	    if(converged) Rprintf(" -->> converged\n");
	    if (trace_lev >= 4) {
		Rprintf("   Ref.step %3d: #{no-improvements}=%3d; (del,dB)=(%12.7g,%12.7g)\n",
			nref, nnoimpr, del, rel_tol * fmax2(rel_tol, nrmB));
		if (trace_lev >= 5) {
		    Rprintf("    weights: "); disp_vec(weights,n);
		    Rprintf("    t2: "); disp_vec(t2,p2);
		    Rprintf("    t1: "); disp_vec(t1,p1);
		    Rprintf("    res2: "); disp_vec(res2,n);
		}
	    }
	}
	/* STEP 5: Update best fit */
	if (sc < *sscale) {
	    COPY(t1, b1, p1);
	    COPY(t2, b2, p2);
	    COPY(res2, res, n);
	    *sscale = sc;
	    if (trace_lev >= 2)
		Rprintf("  Refinement step %3d: better fit, scale: %#10.5g\n",
			nref, sc);
	    nnoimpr = 0;
	} else {
	    if (trace_lev >= 3)
		Rprintf("  Refinement step %3d: no improvement, scale: %#10.5g\n",
			nref, sc);
	    nnoimpr++;
	}
    } // while(.)

    if ( (!converged) & (nref == max_k) )
	warning(_(" M-S estimate: maximum number of refinement steps reached."));

    if (trace_lev >= 1) {
	Rprintf(" Descent procedure: %sconverged (best scale: %.5g, last step: %.5g)\n",
		converged ? "" : "not ", *sscale, sc);
	if (nnoimpr == K_m_s)
	    Rprintf("  The procedure stopped after %d steps because there was no improvement in the last %d steps.\n  To increase this number, use the control parameter 'k.m_s'.\n", nref, nnoimpr);
	else if (trace_lev >= 2)
	    Rprintf("  No improvements in %d out of %d steps.\n", nnoimpr, nref);
	maybe_SHOW_b1_b2;
    }

    return converged;
} /* m_s_descent() */

/* draw a subsample of observations and calculate a candidate           *
 * starting value for S estimates                                       *
 * uses a custom LU decomposition, which acts on the transposed design  *
 * matrix. In case of a singular subsample, the subsample is modified   *
 * until it is non-singular (for ss == TRUE (== 1)).                    *
 *                                                                      *
 * Parts of the algorithm are based on the Gaxpy version of the LU      *
 * decomposition with partial pivoting by                               *
 * Golub G. H., Van Loan C. F. (1996) - MATRIX Computations             */
int subsample(const double x[], const double y[], int n, int m,
	      double *beta, int *ind_space, int *idc, int *idr,
	      double *lu, double *v, int *pivot,
	      double *Dr, double *Dc, int rowequ, int colequ,
	      Rboolean sample, int mts, Rboolean ss, double tol_inv, Rboolean solve)
{
    /* x:         design matrix (n x m)
       y:         response vector
       n:         length of y, nrow of x
       m:         ncol of x  ( == p )
       beta:      [out] candidate parameters (length m)
       ind_space: (required in sample_noreplace, length n)
                  holds the index permutation
       idc:       (required in sample_noreplace, !! length n !!)
		  [out] index of observations used in subsample
       idr:       work array of length m
       lu:        [out] LU decomposition of subsample of xt (m x m)
                  Note: U has is not rescaled by 1 / *cf, as it should,
		        this is done R_subsample().
       v:         work array of length m
       pivot:     [out] pivoting table of LU decomposition (length m-1)
       Dr:        row equilibration (as calculated in SETUP_EQUILIBRATION)
       Dc:        column equilibration
       rowequ:    whether rows were equilibrated
       coleq:     whether cols were equilibrated
       sample:    whether to sample or not
       mts:       the number of singular samples allowed before
                  giving up (Max Try Samples)
       ss:        type of subsampling to be used:
                  0 (FALSE): simple subsampling
                  1  (TRUE): nonsingular subsampling
       tol_inv:   tolerance for declaring a matrix singular
       solve:     solve the least squares problem on the subsample?
                  (0: no, 1: yes)

       return value ('condition'):
             0: success
             1: singular (matrix xt does not contain a m dim. full rank submatrix)
             2: too many singular resamples (simple subsampling case)
*/
    int j, k, l, mu = 0, i = 0, attempt = 0;
    Rboolean sing;

#define xt(_k_, _j_) x[idr[_k_]*n+idc[_j_]]
#define U(_k_, _j_) lu[_j_*m+_k_]
#define u(_k_, _j_) lu + (_j_*m+_k_)
#define L(_k_, _j_) lu[_j_*m+_k_]
#define l(_k_, _j_) lu + (_j_*m+_k_)

Start:
    /* STEP 1: Calculate permutation of 1:n */
    if (sample) {
	sample_noreplace(ind_space, n, n, idc);
    } else // !sample --> "trivial permutation":
	for(k=0;k<n;k++) ind_space[k] = k;
    for(k=0;k<m;k++) idr[k] = k;

    /* STEP 2: Calculate LU decomposition of the first m cols of xt     *
     *         using the order in ind_space                             */
    for(j = 0; j < m; j++) {
	sing=TRUE;
	do {
	    if (i+j == n) {
		warning(_("subsample(): could not find non-singular subsample."));
		return(1);
	    }
	    idc[j] = ind_space[i+j];
	    if (j == 0) {
		for(k=j;k<m;k++) v[k] = xt(k, j);
	    } else {
		for(k=0;k<j;k++) U(k,j) = xt(k, j);
		/* z = solve(lu[0:(j-1), 0:(j-1)], xt[0:(j-1), j]) */
		F77_CALL(dtrsv)("L", "N", "U", &j, lu, &m, u(0, j), &one FCONE FCONE FCONE);
		/* Rprintf("Step %d: z = ", j);  */
		/* for(i=0; i < j; i++) Rprintf("%lf ",U(i, j));  Rprintf("\n"); */

		/* v[j:(m-1)] = xt[j:(m-1), j] - L[j:(m-1), 0:(j-1)] %*% z */
		for(k=j;k<m;k++) {
		    v[k] = xt(k, j);
		    for(l=0;l<j;l++) v[k] -= L(k, l) * U(l, j);
		}
		/* Rprintf("v = "); disp_vec(v, m); */
	    }
	    if (j < m-1) {
		/* find pivot */
		double tmpd = fabs(v[j]); mu = j;
		for(k=j+1; k<m; k++) if (tmpd < fabs(v[k])) { mu = k; tmpd = fabs(v[k]); }
                /* debug possumDiv example, see tests/subsample.R */
		/* if (j == 36) { */
		/*     Rprintf("Step %d: ", j+1); */
		/*     for(k=j;k<m;k++) Rprintf("%lf ", fabs(v[k])); */
		/*     Rprintf("\n %d %lf\n", mu+1, v[mu]); */
		/*     Rprintf("47 > 51: %x\n", fabs(v[46]) > fabs(v[50])); */
		/*     Rprintf("47 < 51: %x\n", fabs(v[46]) < fabs(v[50])); */
		/* } */
		/* continue only if pivot is large enough */
		if (tmpd >= tol_inv) {
		    pivot[j] = mu; // swap  j <--> mu :
		        tmpd =   v[j];   v[j] =   v[mu];   v[mu] = tmpd;
		    int itmp = idr[j]; idr[j] = idr[mu]; idr[mu] = itmp;
		    for(k=j+1;k<m;k++) L(k, j) = v[k] / v[j];
		    if (j > 0) {
			for(k=0;k<j;k++) {
			    tmpd = L(j, k); L(j, k) = L(mu, k); L(mu, k) = tmpd;
			}
		    }
		}
	    }
	    if (fabs(v[j]) < tol_inv) {
		if (ss == 0) {
		    attempt++;
		    if (attempt >= mts) {
			warning(_("Too many singular resamples. Aborting subsample().\n See parameter 'subsampling; in help of lmrob.config()."));
			return(2);
		    }
		    goto Start;
		}
		/* drop observation and try next one */
		i++;
	    } else {
		sing = FALSE;
		U(j, j) = v[j];
	    }
	} while(sing);
    } /* end for loop */

    /* Rprintf("lu:");    disp_vec (lu, m*m); */
    /* Rprintf("pivot:"); disp_veci(pivot, m-1); */
    /* Rprintf("idc:");   disp_veci(idc, m); */

    /* STEP 3: Solve for candidate parameters if requested */
    if (solve == 0) {
      for(k=0;k<m;k++) beta[k] = NA_REAL;
    } else {
      for(k=0;k<m;k++) beta[k] = y[idc[k]];
      /* scale y ( = beta ) */
      if (rowequ) for(k=0;k<m;k++) beta[k] *= Dr[idc[k]];
      /* solve U\tr L\tr \beta = y[subsample] */
      F77_CALL(dtrsv)("U", "T", "N", &m, lu, &m, beta, &one FCONE FCONE FCONE);
      F77_CALL(dtrsv)("L", "T", "U", &m, lu, &m, beta, &one FCONE FCONE FCONE);
      /* scale the solution */
      if (colequ) for(k=0;k<m;k++) beta[k] *= Dc[idr[k]];
      /* undo pivoting : swap beta[k] <-> beta[pivot[k]] */
      for(k=m-2; k>=0; k--) {
	  double tmp = beta[k]; beta[k] = beta[pivot[k]]; beta[pivot[k]] = tmp;
      }
    }

    return(0);

#undef Xt
#undef U
#undef u
#undef L
#undef l
}


double find_scale(const double r[], double b, const double rrhoc[], int ipsi,
		  double initial_scale, int n, int p,
		  int* iter, // input: max_iter,  output: #{iterations used}
		  double scale_tol,
		  Rboolean trace)
{
    if(initial_scale <= 0.) {
	warning(_("find_scale(*, initial_scale = %g <= 0) -> final scale = 0"), initial_scale);
	return 0.;
    }
    // else
    double scale = initial_scale;
    if(trace) Rprintf("find_scale(*, ini.scale =%#13.11g, tol=%g):\n  it | new scale\n",
		      scale, scale_tol);
    for(int it = 0; it < iter[0]; it++) {
	scale *= sqrt( sum_rho_sc(r, scale, n, p, rrhoc, ipsi) / b ) ;
	if(trace) Rprintf("  %2d | %#13.10g\n", it, scale);
	if(fabs(scale - initial_scale) <= scale_tol*initial_scale) { // converged:
	    *iter = it; return(scale);
	}
	initial_scale = scale;
    }
    warning(_("find_scale() did not converge in '%s' (= %d) iterations with tol=%g, last rel.diff=%g"),
	    "maxit.scale", /* <- name from lmrob.control() */ *iter, scale_tol,
	    (scale - initial_scale) / initial_scale);

    return(scale);
}


// As R's which.max(a),  return()ing zero-based   k in {0,1,...,n-1}
int find_max(const double a[], int n)
{
    int k = 0;
    if(n > 1) {
	double tt = a[0];
	for(int i=1; i < n; i++)
	    if(tt < a[i]) {
		tt = a[i];
		k = i;
	    }
    }
    return k;
}

double sum_rho_sc(const double r[], double scale, int n, int p, const double c[], int ipsi)
{
    double s = 0;
    for(int i=0; i < n; i++)
	s += rho(r[i]/scale, c, ipsi);
    return(s / ((double) n - p));
}

/* ||x||_2^2 */
double norm2(const double x[], int n)
{
    double s = 0.;
    s = F77_CALL(dnrm2)(&n, x, &one);
    return( s*s );
}

/* ||x||_2 */
double norm(const double x[], int n)
{
    return(F77_CALL(dnrm2)(&n, x, &one));
}

/* ||x||_1 */
double norm1(const double x[], int n)
{
    return(F77_CALL(dasum)(&n, x, &one));
}

/* ||x-y||_2^2 */
double norm_diff2(const double x[], const double y[], int n)
{
    double s = 0;
    for(int i=0; i < n; i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
    return( s );
}

/* ||x-y||_2 */
double norm_diff(const double x[], const double y[], int n)
{
    double s = 0;
    for(int i=0; i < n; i++)
	s += (x[i]-y[i])*(x[i]-y[i]);
    return( sqrt(s) );
}


/* ||x-y||_1 */
double norm1_diff(const double x[], const double y[], int n)
{
    double s = 0;
    for(int i=0; i < n; i++)
	s += fabs(x[i]-y[i]);
    return(s);
}


double median(const double x[], int n, double *aux)
{
    for(int i=0; i < n; i++) aux[i]=x[i];
    return
	(n % 2) ? /* odd */ kthplace(aux,n, n/2+1)
 	: /* even */ (kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1)) / 2. ;
}

double median_abs(const double x[], int n, double *aux)
{
    for(int i=0; i < n; i++) aux[i] = fabs(x[i]);
    return
	(n % 2) ? /* odd */ kthplace(aux,n, n/2+1)
 	: /* even */ (kthplace(aux,n, n/2) + kthplace(aux,n, n/2+1)) / 2. ;
}

double MAD(const double a[], int n, double center,
	   double *b, // -> the centered  a[] - center
	   double *tmp)
{
/* if center == 0 then do not center */
/*     if( fabs(center) > 0.) { */
	for(int i=0; i < n; i++)
	    b[i] = a[i] - center;
/*     } */
    return( median_abs(b,n,tmp) * 1.4826 );
}

void disp_vec(const double a[], int n)
{
    for(int i=0; i < n; i++) Rprintf("%g ", a[i]);
    Rprintf("\n");
}

void disp_veci(const int a[], int n)
{
    for(int i=0; i < n; i++) Rprintf("%d ", a[i]);
    Rprintf("\n");
}

void disp_mat(const double **a, int n, int m)
{
    for(int i=0; i < n; i++) {
	Rprintf("\n");
	for(int j=0; j < m; j++) Rprintf("%10.8f ",a[i][j]);
    }
    Rprintf("\n");
}

void R_find_D_scale(double *rr, double *kkappa, double *ttau, int *llength,
		    double *sscale, double *cc, int *iipsi, int *ttype, double *rel_tol,
		    int *max_k, int *converged)
{
    /* compute D_scale using iterative algorithm
       type: 1: d1
       2: d2
       3: dt1
       4: dt2
    */
    *converged = 0;
    for (int k=0; k < *max_k; k++) {
	double scale = *sscale, tsum1 = 0, tsum2 = 0;
	// calculate weights
	for(int i=0; i < *llength; i++) {
	    double a, w = wgt(rr[i] / ttau[i] / scale, cc, *iipsi);
	    switch(*ttype) {
	    case 1: // d1
		a = rr[i]/ttau[i];
		tsum1 += a*a*w;
		tsum2 += w; break;
	    case 2: // d2
		a = rr[i]/ttau[i]*w;
		tsum1 += a*a;
		tsum2 += w*w; break;
	    default:
	    case 3: // dt1
		tsum1 += rr[i]*rr[i]*w;
		tsum2 += w*ttau[i]*ttau[i]; break;
	    case 4: // dt2
		a = rr[i]*w;
		tsum1 += a*a;
		a = ttau[i]*w;
		tsum2 += a*a; break;
	    };
	}

	*sscale = sqrt(tsum1 / tsum2 / *kkappa);

	// Rprintf("\n type = %d, scale = %10.8f \n", *ttype, *sscale);

	if (fabs(scale - *sscale) < *rel_tol * fmax2(*rel_tol, scale)) {
	    *converged = 1;
	    break;
	}
    }
}

/* specialized function calc_fitted */
/* calculates fitted values from simulation output array. */
/* this is used to process simulation output in the */
/* lmrob_simulation vignette */
void R_calc_fitted(double *XX, double *bbeta, double *RR, int *nn, int *pp, int *nnrep,
		   int *nnproc, int *nnerr)
{
    unsigned long A, B, C, D, E;
    A = (unsigned long)*nnerr; B = (unsigned long)*nnproc; C = (unsigned long)*nnrep;
    D = (unsigned long)*nn; E = (unsigned long)*pp;
    // calculate fitted values over errstr, procstr and replicates
    for(unsigned long a = 0; a < A; a++) { // errstr
	for(unsigned long b = 0; b < B; b++) { // procstr
	    for(unsigned long c = 0; c < C; c++) { // replicates
		// check for NAs
		if (!ISNA(bbeta[c + /* 0*C + */  b*C*E + a*B*E*C])) {
		    for(unsigned long d = 0; d < D; d++) { // observations
			RR[d + c*D + b*C*D + a*B*C*D] = 0; // initialize result
			for(unsigned long e = 0; e < E; e++) { // predictors
			    RR[d + c*D + b*C*D + a*B*C*D] += bbeta[c + e*C + b*C*E + a*B*E*C] *
				XX[d + e*D + c*E*D + a*C*E*D];
			}
		    }
		}
	    }
	}
    }
}
