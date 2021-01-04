/*------ Definition of a template for whimed(_i) : *
 *			 --------     ~~~~~~
 * i.e., included several times from ./wgt_himed.c
 *				     ~~~~~~~~~~~~~
 */


#ifdef _d_whimed_

# define _WHIMED_	whimed
# define _WGT_TYPE_	double
# define _WGT_SUM_TYPE_ double
# undef _d_whimed_

#elif defined (_i_whimed_)

# define _WHIMED_	whimed_i
# define _WGT_TYPE_	int
# define _WGT_SUM_TYPE_ int64_t
# undef _i_whimed_

#else
# error "must define correct  whimed_  macro !"
#endif


double _WHIMED_(double *a, _WGT_TYPE_ *w, int n,
		double* a_cand, double *a_srt, _WGT_TYPE_* w_cand)
{

/*
  Algorithm to compute the weighted high median in O(n) time.

  The whimed is defined as the smallest a[j] such that the sum
  of the weights of all a[i] <= a[j] is strictly greater than
  half of the total weight.

  Arguments:

  a: double array containing the observations
  n: number of observations
  w: array of (int/double) weights of the observations.
*/
    int i;
    /* sum of weights: `int' do overflow when  n ~>= 1e5 */
    _WGT_SUM_TYPE_ wleft, wmid, wright, w_tot, wrest;
    double trial;

    w_tot = wrest = 0;
    for (i = 0; i < n; ++i)
	w_tot += w[i];

#ifdef DEBUG_whimed
    REprintf("wgt_himed(a[], w[], n)  -- on entry: n=%d, w_tot=%g\n", n, (double)w_tot);
#endif
    if(n == 0) return NA_REAL;

/* REPEAT : */
    do {
	int n2 = n/2;/* =^= n/2 +1 with 0-indexing */
	for (i = 0; i < n; ++i)
	    a_srt[i] = a[i];
	rPsort(a_srt, n, n2);
	trial = a_srt[n2];

	wleft = 0;    wmid  = 0;    wright= 0;
	for (i = 0; i < n; ++i) {
	    if (a[i] < trial)
		wleft += w[i];
	    else if (a[i] > trial)
		wright += w[i];
	    else
		wmid += w[i];
	}
	/* wleft = sum_{i; a[i]	 < trial}  w[i]
	 * wmid	 = sum_{i; a[i] == trial}  w[i] at least one 'i' since trial is one a[]!
	 * wright= sum_{i; a[i]	 > trial}  w[i]
	 */
#ifdef DEBUG_whimed
	REprintf(" trial=%-g; w(left|mid|right) = (%g,%g,%g); ", trial,
		 (double)wleft, (double)wmid, (double)wright);
#endif

	int kcand = 0;
	if (2 * (wrest + wleft) > w_tot) {
	    for (i = 0; i < n; ++i) {
		if (a[i] < trial) {
		    a_cand[kcand] = a[i];
		    w_cand[kcand] = w[i];	++kcand;
		}
	    }
	}
	else if (2 * (wrest + wleft + wmid) <= w_tot) {
	    for (i = 0; i < n; ++i) {
		if (a[i] > trial) {
		    a_cand[kcand] = a[i];
		    w_cand[kcand] = w[i];	++kcand;
		}
	    }
	    wrest += wleft + wmid;
#ifdef DEBUG_whimed
	    REprintf(" new wrest = %g; ", (double)wrest);
#endif
	}
	else {
#ifdef DEBUG_whimed
	    REprintf(" -> found! return trial\n");
#endif
	    return trial;
	    /*==========*/
	}
	n = kcand;
#ifdef DEBUG_whimed
	REprintf("  ... and try again with  n:= kcand=%d\n", n);
#endif
	for (i = 0; i < n; ++i) {
	    a[i] = a_cand[i];
	    w[i] = w_cand[i];
	}
    } while(1);

} /* _WHIMED_ */

#undef _WHIMED_
#undef _WGT_TYPE_
#undef _WGT_SUM_TYPE_
