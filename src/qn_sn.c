/*
 *  Copyright (C) 2005--2007, 2021	Martin Maechler, ETH Zurich
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* This is a merge of the C version of original files  qn.f and sn.f,
 * translated by f2c (version 20010821).	       ====	====
 * and then by f2c-clean,v 1.9 2000/01/13 13:46:53
 * and further clean-edited manually by Martin Maechler.
 *
 * Further added interface functions to be called via .C() from R or S-plus
 * Note that Peter Rousseeuw has explicitely given permission to
 * use his code under the GPL for the R project.
 */

/* Original comments by the authors of the Fortran original code,
 * (merged for Qn & Sn in one file by M.M.):

   This file contains fortran functions for two new robust estimators
   of scale denoted as Qn and Sn, decribed in Rousseeuw and Croux (1993).
   These estimators have a high breakdown point and a bounded influence
   function. The implementation given here is very fast (running in
   O(n logn) time) and needs little storage space.

	Rousseeuw, P.J. and Croux, C. (1993)
	Alternatives to the Median Absolute Deviation",
	Journal of the American Statistical Association, Vol. 88, 1273-1283.

   For both estimators, implementations in the pascal language can be
   obtained from the original authors.

   This software may be used and copied freely for scientific
   and/or non-commercial purposes, provided reference is made
   to the abovementioned paper.

Note by MM: We have explicit permission from P.Rousseeuw to
licence it under the GNU Public Licence.

See also ../inst/Copyrights
*/
#include <inttypes.h>
/*        ^^^^^^^^^^ is supposedly more common and standard than
 * #include <stdint.h>
 * or #include <sys/types.h> */
/* --> int64_t ; if people don't have the above, they can forget about it.. */
/* #include "int64.h" */

#include <Rmath.h> /* -> <math.h> and much more */

/* Interface routines to be called via .C() : */
#include "robustbase.h"

/* ----------------- Further Declarations ------------------------------ */

/* sn0() and  qn0()  --- but also  mc_C()  in ./mc.c
 * -----      ----                 ------
 use

   pull(a,n,k): finds the k-th order statistic of an
		array a[] of length n (preserving a[])
*/

/*
   whimed_i(a,iw,n): finds the weighted high median of an array
		a[] of length n, with positive int weights iw[]
		(using auxiliary arrays acand[], a_srt[] & iw_cand[]
		 all of length n).
*/

/* qn0() uses (and for C API:) */

/* Main routines for C API */
double qn(double *x, int n, int h,         int finite_corr);
double sn(double *x, int n, int is_sorted, int finite_corr);
/* these have no extra factors (no consistency factor & finite_corr): */
void  qn0(const double x[], int n, const int64_t k[], int len_k, /* ==> */ double *res);

double sn0(double *x, int n, int is_sorted, double *a2);


/* ----------- Implementations -----------------------------------*/

// === called from R ( ../R/qnsn.R ) via .C() :

void Qn0(double *x, Sint *n, double *k, Sint *len_k, double *res) {
    int l_k = (int)*len_k;
    // "hack" as R / .C() have no int_64 : copy k[] to int64 ik[]:
    int64_t *ik  = (int64_t *) R_alloc(l_k, sizeof(int64_t));
    for(int i=0; i < l_k; i++) ik[i] = (int64_t) k[i];
    qn0(x, (int)*n, ik, l_k, res);
}

void Sn0(double *x, Sint *n, Sint *is_sorted, double *res, double *a2)
{
    char *vmax = vmaxget();

    *res = sn0(x, (int)*n, (int)*is_sorted, a2);
#ifdef DEBUG_Sno
    REprintf("Sn0(* -> res=%g)\n", *res);
#endif

    vmaxset(vmax);
}

void qn0(const double x[], int n, const int64_t k[], int len_k, /* ==> */ double *res)
{
/*--------------------------------------------------------------------

   Efficient algorithm for the scale estimator:

       Q*_n = { |x_i - x_j|; i<j }_(k)	[= Qn without scaling ]

		i.e. the k-th order statistic of the |x_i - x_j|

   Parameters of the function Qn :
       x  : double array containing the observations
       n  : number of observations (n >= 2)
 */

    double *y	  = (double *)R_alloc(n, sizeof(double));
    double *work  = (double *)R_alloc(n, sizeof(double));
    double *a_srt = (double *)R_alloc(n, sizeof(double));
    double *a_cand = (double *)R_alloc(n, sizeof(double));

    int *left	  = (int *)R_alloc(n, sizeof(int));
    int *right	  = (int *)R_alloc(n, sizeof(int));
    int *p	  = (int *)R_alloc(n, sizeof(int));
    int *q	  = (int *)R_alloc(n, sizeof(int));
    int *weight	  = (int *)R_alloc(n, sizeof(int));

    const int64_t
	nn2 = (int64_t) n * (n + 1) / 2, // = choose(n+1, 2)
	n2  = (int64_t) n * n,
	// k >= k_L  ==> right[] = n (and smaller k  can start with smaller right[], hence more efficiently)
	/* NB: The correct k_L seems hard to find, till only from largish simulations and "regression fit":
	 * 2021, Jan.  3 -- 19       :  k_L = 2 + ((int64_t)(n-1) * (n-2))/2;
	 * 2021, Jan. 19 -- 27, 11:13:  k_L =     ((int64_t)(n-1) * (n-2))/2;
	 * 2021, Jan. 27 -- ... : from .../Pkg-ex/robustbase/Qn-multi-k-debugging.R
                              bnd1 <- function(n) 5 - 1.75*(n %% 2) + (0.3939 - 0.0067*(n %% 2)) * n*(n-1)
	*/
	k_L = 5 - 1.75*(n % 2) + (0.3939 - 0.0067*(n % 2)) * (int64_t) n*(n-1) ;
    int h = n / 2 + 1; // really use, only to set right[] below ?
    for (int i = 0; i < n; ++i)
	y[i] = x[i];
    R_qsort(y, 1, n); /* y := sort(x) */

#ifdef DEBUG_qn
    REprintf("qn0(x, n=%2d, .. |k|=%d): \n", n, len_k);
#endif

    // k = (int64_t)h * (h - 1) / 2;
  for(int i_k=0; i_k < len_k; i_k++) { // compute k'th quantile k' := k[i_k] ------------------
    /* Following should be `long long int' : they can be of order n^2 */
    int64_t nl = nn2, nr = n2, knew = k[i_k] + nl;/* = k + (n+1 \over 2) */
#ifdef DEBUG_qn
    REprintf(" qn0(x, ..): nl,nr= %d, %d;  k[%d] = %5d :\n", nl,nr, i_k, k[i_k]);
#endif
/* L200: */
    Rboolean found = FALSE;
    double trial = R_NaReal;/* -Wall */
    int i, j;

    for (int i = 0; i < n; ++i)
	left [i] = n - i + 1;
    if(k[i_k] >= k_L) { // for these large "quantiles", need "large" right[] boundary
	for (int i = 0; i < n; ++i)
	    right[i] = n;
    } else {
	for (int i = 0; i < n; ++i)
	    right[i] = (i <= h) ? n : n - (i - h);
	/* the n - (i-h) is from the paper (Cr. + Ro. '92 "Time-efficient"); original code had 'n' */
    }

    while(!found && nr - nl > n) {
	j = 0;
	/* Truncation to float :
	   try to make sure that the same values are got later (guard bits !) */
	for (i = 1; i < n; ++i) {
	    if (left[i] <= right[i]) {
		weight[j] = right[i] - left[i] + 1;
		int jh = left[i] + weight[j] / 2;
		work[j] = (float)(y[i] - y[n - jh]);
		++j;
	    }
	}
	trial = whimed_i(work, weight, j, a_cand, a_srt, /*iw_cand*/ p);

#ifdef DEBUG_qn
	REprintf(" ..!found: whimed(");
#  ifdef DEBUG_long
	REprintf("wrk=c(");
	for(i=0; i < j; i++) REprintf("%s%g", (i>0)? ", " : "", work[i]);
	REprintf("),\n	   wgt=c(");
	for(i=0; i < j; i++) REprintf("%s%d", (i>0)? ", " : "", weight[i]);
	REprintf("), j= %3d) -> trial= %7g\n", j, trial);
#  else
	REprintf("j=%3d) -> trial= %g:", j, trial);
#  endif
#endif
	j = 0;
	for (i = n - 1; i >= 0; --i) {
	    while (j < n && ((float)(y[i] - y[n - j - 1])) < trial)
		++j;
	    p[i] = j;
	}
#ifdef DEBUG_qn
	REprintf(" f_1: j=%2d", j);
#endif
	j = n + 1;
	for (i = 0; i < n; ++i) {
	    while ((float)(y[i] - y[n - j + 1]) > trial)
		--j;
	    q[i] = j;
	}
	int64_t
	    sump = 0,
	    sumq = 0;
	for (i = 0; i < n; ++i) {
	    sump += p[i];
	    sumq += q[i] - 1;
	}
#ifdef DEBUG_qn
	REprintf(" f_2 -> j=%2d, sump|q= %lld,%lld; ", j, sump,sumq);
#endif
	if (knew <= sump) {
	    for (i = 0; i < n; ++i)
		right[i] = p[i];
	    nr = sump;
#ifdef DEBUG_qn
	    REprintf("knew <= sump =: nr , new right[]\n");
#endif
	} else if (knew > sumq) {
	    for (i = 0; i < n; ++i)
		left[i] = q[i];
	    nl = sumq;
#ifdef DEBUG_qn
	    REprintf("knew > sumq =: nl , new left[]\n");
#endif
	} else { /* sump < knew <= sumq */
	    found = TRUE;
#ifdef DEBUG_qn
	    REprintf("sump < knew <= sumq ---> FOUND\n");
#endif
	}
    } /* while */

    if (found)
	res[i_k] = trial;
    else {
#ifdef DEBUG_qn
	REprintf(".. not fnd -> new work[]");
#endif
	j = 0;
	for (i = 1; i < n; ++i) {
	    for (int jj = left[i]; jj <= right[i]; ++jj) {
		work[j] = y[i] - y[n - jj];
		j++;
	    }/* j will be = sum_{i=2}^n (right[i] - left[i] + 1)_{+}  */
	}
#ifdef DEBUG_qn
	REprintf(" of length %d; k' = knew-nl = %d\n", j, knew-nl);
#endif

	/* return pull(work, j - 1, knew - nl)	: */
	knew -= (nl + 1); /* -1: 0-indexing */

	if(knew  > j-1) { // see this happening when the quantile number k[i_k] is close to the right end!
	    knew = j-1;
#ifdef DEBUG_qn
	    REprintf("knew >= j should never happen, setting it to j-1=%d\n", knew);
#endif
	} else if(knew < 0) {
	    knew = 0;
#ifdef DEBUG_qn
	    REprintf("knew < 0 should never happen, setting it to 0\n");
#endif
	}
	rPsort(work, j, knew);
	res[i_k] = work[knew];
    }
  } // for(int i_k=0, i_k < len_k ...) { k_ = k[i_k] ; ....
  return;
} /* qn0 */


#ifdef __never__ever__
/* NB: "old version" -- with *wrong* constant and "inaccurate" finite-sample corr
 * --  equivalent to Qn.old() in ../R/qnsn.R , see also ../man/Qn.Rd
 *
 * (This is *not* called from our R code anyway)
 */
double qn(double *x, int n, int finite_corr)
{
/* Efficient algorithm for the scale estimator:

	Qn = dn * 2.2219 * {|x_i-x_j|; i<j}_(k)
*/
    double r, dn = 1./*-Wall*/;

    r = 2.2219 * qn0(x, n); /* asymptotic consistency for sigma^2 */
    /*		 === */

    if (finite_corr) {
	if (n <= 9) {
	    if	    (n == 2)	dn = .399;
	    else if (n == 3)	dn = .994;
	    else if (n == 4)	dn = .512;
	    else if (n == 5)	dn = .844;
	    else if (n == 6)	dn = .611;
	    else if (n == 7)	dn = .857;
	    else if (n == 8)	dn = .669;
	    else if (n == 9)	dn = .872;
	} else {
	    if (n % 2 == 1)
		dn = n / (n + 1.4);
	    else /* (n % 2 == 0) */
		dn = n / (n + 3.8);
	}
	return dn * r;
    }
    else return r;
} /* qn */
#endif


double sn0(double *x, int n, int is_sorted, double *a2)
{
/*
   Efficient algorithm for the scale estimator:

       S*_n = LOMED_{i} HIMED_{j} |x_i - x_j|

   which can equivalently be written as

       S*_n = LOMED_{i} LOMED_{j != i} |x_i - x_j|

   Arguments :

       x  : double array (length >= n) containing the observations
       n  : number of observations (n>=2)

       is_sorted: logical indicating if x is already sorted
       a2 : to contain	a2[i] := LOMED_{j != i} | x_i - x_j |,
	    for i=1,...,n
*/

    /* Local variables */
    double medA, medB;

    int i, diff, half, Amin, Amax, even, length;
    int leftA,leftB, nA,nB, tryA,tryB, rightA,rightB;
    int n1_2;

    if(!is_sorted)
	R_qsort(x, 1, n);

    a2[0] = x[n / 2] - x[0];
    n1_2 = (n + 1) / 2;

    /* first half for() loop : */
    for (i = 2; i <= n1_2; ++i) {
	nA = i - 1;
	nB = n - i;
	diff = nB - nA;
	leftA  = leftB	= 1;
	rightA = rightB = nB;
	Amin = diff / 2 + 1;
	Amax = diff / 2 + nA;

	while (leftA < rightA) {
	    length = rightA - leftA + 1;
	    even = 1 - length % 2;
	    half = (length - 1) / 2;
	    tryA = leftA + half;
	    tryB = leftB + half;
	    if (tryA < Amin) {
		rightB = tryB;
		leftA = tryA + even;
	    }
	    else {
		if (tryA > Amax) {
		    rightA = tryA;
		    leftB = tryB + even;
		}
		else {
		    medA = x[i - 1] - x[i - tryA + Amin - 2];
		    medB = x[tryB + i - 1] - x[i - 1];
		    if (medA >= medB) {
			rightA = tryA;
			leftB = tryB + even;
		    } else {
			rightB = tryB;
			leftA = tryA + even;
		    }
		}
	    }
	} /* while */

	if (leftA > Amax) {
	    a2[i - 1] = x[leftB + i - 1] - x[i - 1];
	} else {
	    medA = x[i - 1] - x[i - leftA + Amin - 2];
	    medB = x[leftB + i - 1] - x[i - 1];
	    a2[i - 1] = fmin2(medA,medB);
	}
    }

    /* second half for() loop : */
    for (i = n1_2 + 1; i <= n - 1; ++i) {
	nA = n - i;
	nB = i - 1;
	diff = nB - nA;
	leftA  = leftB	= 1;
	rightA = rightB = nB;
	Amin = diff / 2 + 1;
	Amax = diff / 2 + nA;

	while (leftA < rightA) {
	    length = rightA - leftA + 1;
	    even = 1 - length % 2;
	    half = (length - 1) / 2;
	    tryA = leftA + half;
	    tryB = leftB + half;
	    if (tryA < Amin) {
		rightB = tryB;
		leftA = tryA + even;
	    } else {
		if (tryA > Amax) {
		    rightA = tryA;
		    leftB = tryB + even;
		} else {
		    medA = x[i + tryA - Amin] - x[i - 1];
		    medB = x[i - 1] - x[i - tryB - 1];
		    if (medA >= medB) {
			rightA = tryA;
			leftB = tryB + even;
		    } else {
			rightB = tryB;
			leftA = tryA + even;
		    }
		}
	    }
	} /* while */

	if (leftA > Amax) {
	    a2[i - 1] = x[i - 1] - x[i - leftB - 1];
	} else {
	    medA = x[i + leftA - Amin] - x[i - 1];
	    medB = x[i - 1] - x[i - leftB - 1];
	    a2[i - 1] = fmin2(medA,medB);
	}
    }
    a2[n - 1] = x[n - 1] - x[n1_2 - 1];

    return pull(a2, n, n1_2);
} /* sn0 */


// C API version -- *not* called from our R code :
double sn(double *x, int n, int is_sorted, int finite_corr)
{
/*
   Efficient algorithm for the scale estimator:

       Sn = cn * 1.1926 * LOMED_{i} HIMED_{j} |x_i-x_j|

   which can equivalently be written as

       Sn = cn * 1.1926 * LOMED_{i} LOMED_{j != i} |x_i-x_j|*/

    double cn, r;
    double *a2 = (double *)R_alloc(n, sizeof(double));

    r = 1.1926 * /* asymptotic consistency for sigma^2 */
	sn0(x, n, is_sorted, a2);
    /*	=== */

    cn = 1.; /* n >= 10 even, or no finite_corr[ection] */
    if (finite_corr) {
	if (n <= 9) {
	    if	    (n == 2)	cn = 0.743;
	    else if (n == 3)	cn = 1.851;
	    else if (n == 4)	cn = 0.954;
	    else if (n == 5)	cn = 1.351;
	    else if (n == 6)	cn = 0.993;
	    else if (n == 7)	cn = 1.198;
	    else if (n == 8)	cn = 1.005;
	    else if (n == 9)	cn = 1.131;
	} else if (n % 2 == 1) /* n odd, >= 11 */
	    cn = n / (n - 0.9);
    }
    return cn * r;

} /* sn */


/* pull():   auxiliary routine for (Qn and) Sn
 * ======    ========  ------------------------
 */
double pull(double *a_in, int n, int k)
{
/* Finds the k-th order statistic of an array a[] of length n
 *	     --------------------
*/
    int j;
    double *a, ax;
    char* vmax = vmaxget();
    a = (double *)R_alloc(n, sizeof(double));
    /* Copy a[] and use copy since it will be re-shuffled: */
    for (j = 0; j < n; j++)
	a[j] = a_in[j];

    k--; /* 0-indexing */
    rPsort(a, n, k);
    ax = a[k];

    vmaxset(vmax);
    return ax;
} /* pull */

/* Local variables section

 * Local variables:
 * mode: c
 * kept-old-versions: 12
 * kept-new-versions: 20
 * End:
 */
