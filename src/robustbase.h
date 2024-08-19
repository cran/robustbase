/* External and interal  API  of  C and Fortran routines in robustbase */

#define STRICT_R_HEADERS

#include <R.h>
// for SEXP
#include <Rinternals.h>

/**< For internationalized messages */
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("Matrix", String)
#else
#define _(String) (String)
#define dngettext(pkg, String, StringP, N) (N > 1 ? StringP : String)
#endif



/* --------- ./qn_sn.c : -------- */
#define Sint int

void Qn0(double *x, Sint *n, double *k, Sint *len_k, double *res);

void Sn0(double *x, Sint *n, Sint *is_sorted, double *res, double *a2);
/*
 * void Qn    (double *x, Sint *n, Sint *h, Sint *finite_corr, double *res);
 * void Sn    (double *x, Sint *n,          Sint *finite_corr, double *res);
*/

/* call via .C() from R : */
void wgt_himed_i(double *x, Sint *n,  Sint *iw, double *res);
void wgt_himed  (double *x, Sint *n, double *w, double *res);

/* call from C: */
double pull(double *a, int n, int k);

double whimed_i(double *a, int *iw, int n,
		double *acand, double *a_srt, int *iw_cand);
double whimed(double *a, double *w, int n,
	      double *acand, double *a_srt, double *w_cand);

/* --------- ./mc.c -------- */

/* call via .C() from R : */
void mc_C(double *z, int *in, double *eps, int *iter, double *out, int *scale);

/* call from C: *iter is both input and output */
double mc_C_d(const double z[], int n, const double eps[], int *iter, int scale);




/* --------- ./lmrob.c --------- */

static inline
Rboolean is_redescender(int ipsi) {// a simple wrapper for readability
    // for now, fastest:
    if(ipsi == 0)
	return FALSE;
    return TRUE;
/* if have many more, maybe
    switch(ipsi) {
    default: error(_("ipsi=%d not implemented."), ipsi);
    case 0: // huber and future other non-redescenders
	return FALSE;
    case 1: case 2: case 3: case 4: case 5: case 6:
	return TRUE;
    }
*/
}

SEXP R_rho_inf(SEXP cc, SEXP ipsi);

void R_lmrob_S(double *X, double *y, int *n, int *P,
	       int *nRes, double *scale, double *beta_s,
	       double *C, int *iipsi, double *bb,
	       int *best_r, int *Groups, int *N_group,
	       int *K_s, int *max_k, int *max_it_scale,
	       double *rel_tol, double *inv_tol, double *scale_tol,
               //     ^^^^^^^^^ = refine.tol in R
	       double *zero_tol,
	       int* converged, int *trace_lev, int *mts, int *ss, int *cutoff);

void R_lmrob_M_S(double *X1, double *X2, double *y, double *res,
		 int *n, int *p1, int *p2, int *nRes, int *max_it_scale,
		 double *scale, double *b1, double *b2,
		 double *rho_c, int *ipsi, double *bb,
		 int *K_m_s, int *max_k, double *rel_tol, double *inv_tol, double *scale_tol,
		 double *zero_tol,
		 int *converged, int *trace_lev,
		 int *orthogonalize, int *subsample,
		 int *descent, int *mts, int *ss);

void R_lmrob_MM(double *X, double *y, int *n, int *P,
		double *beta_initial, double *scale,
		double *beta_m, double *resid,
		int *max_it, double *rho_c, int *ipsi,
		double *loss, double *rel_tol,
		int *converged, int *trace_lev);

void R_subsample(const double *x, const double *y, int *n, int *m,
		 double *beta, int *ind_space, int *idc, int *idr,
		 double *lu, double *v, int *p,
		 double *_Dr, double *_Dc, int *_rowequ, int *_colequ,
		 int *status, int *sample, int *mts, int *ss, double *tol_inv,
		 int *solve);

SEXP R_psifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_);
SEXP R_chifun(SEXP x_, SEXP c_, SEXP ipsi_, SEXP deriv_);
SEXP R_wgtfun(SEXP x_, SEXP c_, SEXP ipsi_);

double rho(double x, const double c[], int ipsi);
double psi(double x, const double c[], int ipsi);
double psip(double x, const double c[], int ipsi);// psi'
double psi2(double x, const double c[], int ipsi);// psi''
double wgt(double x, const double c[], int ipsi);
double rho_inf (const double c[], int ipsi); // == \rho(\infty)
double normcnst(const double c[], int ipsi); // == 1 / \rho(\infty) ==  1 / rho_inf()


void R_find_D_scale(double *rr, double *kkappa, double *ttau, int *llength,
		    double *sscale, double *cc, int *iipsi, int *ttype, double *rel_tol,
		    int *max_k, int *converged);

void R_calc_fitted(double *XX, double *bbeta, double *RR, int *nn, int *pp, int *nnrep,
		   int *nnproc, int *nnerr);

// ------- ./rob-utils.c ---------------
SEXP R_wgt_flex(SEXP x_, SEXP c_, SEXP h_);

// ------- ./rowMedians.c ---------------
SEXP R_rowMedians(SEXP x, SEXP naRm, SEXP hasNA, SEXP byRow, SEXP keepNms);

/* ------- ./rffastmcd.f ------------ */
void F77_NAME(rffastmcd)(
    double *dat, int *n, int *nvar, int *nhalff, int *krep,          //  5
    int *nmini, int *kmini,                                          //  7
    double *initcov, double *initmean, int *inbest, double *det,     // 11
    int *weight, int *fit, double *coeff, int *kount, double *adcov, // 16
    int *temp, int *index1, int *index2, int *indexx, double *nmahad,// 21
    double *ndist, double *am, double *am2,                          // 24
    double *slutn, double *med, double *mad, double *sd,             // 28
    double *means, double *bmeans, double *w, double *fv1,           // 32
    double *fv2, double *rec, double *sscp1, double *cova1,          // 36
    double *corr1, double *cinv1, double *cova2,                     // 39
    double *cinv2, double *z__, double *cstock, double *mstock,      // 43
    double *c1stock, double *m1stock, double *dath,                  // 46
    double *cutoff, double *chimed, int *i_trace);                   // 49 args

/* ------- ./rfltsreg.f ------------ */
void F77_NAME(rfltsreg)(
    double *dat, int *n, int *nvar,
    int *nhalff, int *krep, int *inbest, double *objfct,
    int *intercept, int *intadjust, int *nvad, double *datt, double *weights,
    int *temp, int *index1, int *index2, double *aw2, double *aw,
    double *residu, double *y, double *nmahad, double *ndist,
    double *am, double *am2, double *slutn,   int *jmiss,
    double *xmed, double *xmad, double *a, double *da,
    double *h__, double *hvec,
    double *c__, double *cstock, double *mstock, double *c1stock, double *m1stock,
    double *dath, double *sd, double *means, double *bmeans, int *i_trace);

/* ------- ./rllarsbi.f -------------- */
void F77_NAME(rllarsbi)(
    double *X, double *Y, int *N, int *NP, int *MDX, int *MDT,
    double *TOL, int *NIT, int *K, int *KODE, double *SIGMA, double *THETA,
    double *RS, double *SC1, double *SC2, double *SC3, double *SC4,
    const double *BET0);
