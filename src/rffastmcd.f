cc -*- mode: fortran; kept-new-versions: 25; kept-old-versions: 20 -*-
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  rrcov : Scalable Robust Estimators with High Breakdown Point
cc
cc  This program is free software; you can redistribute it and/or modify
cc  it under the terms of the GNU General Public License as published by
cc  the Free Software Foundation; either version 2 of the License, or
cc  (at your option) any later version.
cc
cc  This program is distributed in the hope that it will be useful,
cc  but WITHOUT ANY WARRANTY; without even the implied warranty of
cc  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
cc  GNU General Public License for more details.
cc
cc  You should have received a copy of the GNU General Public License
cc  along with this program; if not, write to the Free Software
cc  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
cc
cc  I would like to thank Peter Rousseeuw and Katrien van Driessen for
cc  providing the initial code of this function.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   Computes the MCD estimator of multivariate location and scatter.
cc   This estimator is given by the subset of h observations for which
cc   the determinant of their covariance matrix is minimal. The MCD
cc   location estimate is then the mean of those h points, and the MCD
cc   scatter estimate is their covariance matrix. This value of h may be
cc   chosen by the user; its default value is roughly n/2.
cc
cc   The MCD estimator was first introduced in:
cc
cc	Rousseeuw, P.J. (1984), "Least Median of Squares Regression,"
cc	Journal of the American Statistical Association, Vol. 79,
cc	pp. 871-881. [See page 877.]
cc
cc   The MCD is a robust estimator in the sense that the estimates are
cc   not unduly influenced by outliers in the data, even if there
cc   are many outliers. Its robustness was proved in:
cc
cc	Rousseeuw, P.J. (1985), "Multivariate Estimation with High
cc	Breakdown Point," in Mathematical Statistics and Applications,
cc	edited by  W. Grossmann, G. Pflug, I. Vincze, and W. Wertz.
cc	Dordrecht: Reidel Publishing Company, pp. 283-297.
cc
cc	Rousseeuw, P.J. and Leroy, A.M. (1987), Robust Regression and
cc	Outlier Detection, Wiley-Interscience, New York. [Chapter 7]
cc
cc   The program also computes the distance of each observation
cc   from the center (location) of the data, relative to the shape
cc   (scatter) of the data:
cc
cc   * Using the classical estimates yields the Mahalanobis distance
cc     MD(i). Often, outlying points fail to have a large Mahalanobis
cc     distance because of the masking effect.
cc
cc   * Using the MCD estimates yields a robust distance RD(i).
cc     These distances allow us to easily identify the outliers.
cc
cc   For applications of robust distances in a regression context see:
cc
cc	Rousseeuw, P.J. and van Zomeren, B.C. (1990), "Unmasking
cc	Multivariate Outliers and Leverage Points," Journal of the
cc	American Statistical Association, Vol. 85, 633-639.
cc
cc   There also a diagnostic plot is given to distinguish between
cc   regular observations, vertical outliers, good leverage points,
cc   and bad leverage points.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc   The new FAST_MCD algorithm introduced here is due to
cc
cc	Rousseeuw, P.J. and Van Driessen, K. (1997), "A Fast
cc	Algorithm for the Minimum Covariance Determinant
cc	Estimator," in preparation.
cc
cc   The algorithm works as follows:
cc
cc	 The dataset contains n cases, and nvar variables are used.
cc	 When n < 2*nmini, the algorithm will analyze the dataset as a whole.
cc	 When n >= 2*nmini, the algorithm will use several subdatasets.
cc
cc	 When the dataset is analyzed as a whole, a trial
cc	 subsample of nvar+1 cases is taken, of which the mean and
cc	 covariance matrix is calculated. The h cases with smallest
cc	 relative distances are used to calculate the next mean and
cc	 covariance matrix, and this cycle is repeated k1 times.
cc	 [For small n we can consider all subsets of nvar+1 out of n, else
cc	 the algorithm draws 500 random subsets.]
cc	 Afterwards, the best 10 solutions (covariance matrices and
cc	 corresponding means) are used as starting values for the final
cc	 iterations. These iterations stop when two subsequent determinants
cc	 become equal. (At most k3 iteration steps are taken.)
cc	 The solution with smallest determinant is retained.
cc
cc	 When the dataset contains more than 2*nmini cases, the algorithm
cc	 does part of the calculations on (at most) kmini nonoverlapping
cc	 subdatasets, of (roughly) nmini cases.
cc
cc	 Stage 1: For each trial subsample in each subdataset,
cc	 k1 iterations are carried out in that subdataset.
cc	 For each subdataset, the 10 best solutions are stored.
cc
cc	 Stage 2 considers the union of the subdatasets, called the
cc	 merged set. (If n is large, the merged set is a proper subset of
cc	 the entire dataset.) In this merged set, each of the 'best
cc	 solutions' of stage 1 are used as starting values for k2
cc	 iterations. Also here, the 10 best solutions are stored.
cc
cc	 Stage 3 depends on n, the total number of cases in the
cc	 dataset. If n <= 5000, all 10 preliminary solutions are iterated
cc	 k3 times. If n > 5000, only the best preliminary
cc	 solution is iterated, and the number of iterations decreases to 1
cc	 according to n*nvar. (If n*nvar <= 100,000 we iterate k3 times,
cc	 whereas for n*nvar > 1,000,000 we take only one iteration step.)
cc
cc   An important advantage of the algorithm FAST_MCD is that it allows
cc   for exact fit situations, where more than h observations lie on
cc   a hyperplane. Then the program still yields the MCD location and
cc   scatter matrix, the latter being singular (as it should be), as
cc   well as the equation of the hyperplane.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine rffastmcd(dat,n,nvar,nhalff,krep,initcov,initmean,
     *	  inbest,det,weight,fit,coeff,kount,adcov,
     *	  iseed,
     *	  temp, index1, index2, nmahad, ndist, am, am2, slutn,
     *	  med, mad, sd, means, bmeans, w, fv1, fv2,
     *	  rec, sscp1, cova1, corr1, cinv1, cova2, cinv2, z,
     *	  cstock, mstock, c1stock, m1stock, dath,
     *	  cutoff, chimed)

cc	VT::10.10.2005 - a DATA operator was used for computing the
cc		median and the 0.975 quantile of the chisq distribution
cc		with nvar degrees of freedom. Since now we have no
cc		restriction on the number of variables, these will be
cc		passed as parameters - cutoff and chimed

cc
	implicit integer(i-n), double precision(a-h,o-z)
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  ALGORITHM PARAMETERS:
cc
cc	To change the number of subdatasets and their size, the values of
cc	kmini and nmini can be changed.
cc
	parameter (kmini=5)
	parameter (nmini=300)
cc
cc	The number of iteration steps in stages 1,2 and 3 can be changed
cc	by adapting the parameters k1, k2, and k3.
cc
	parameter (k1=2)
	parameter (k2=2)
	parameter (k3=100)

c MM: below, there is also hardcoded "10" in places related to "stock" !

cc
cc  The parameter krep represents the total number of trial subsamples
cc  to be drawn when n exceeds 2*nmini.
cc
c	parameter (krep=500)
cc
cc  The following lines need not be modified.
cc
C	parameter (nvmax1=nvmax+1)
C	parameter (nvmax2=nvmax*nvmax)
C	parameter (nvm12=nvmax1*nvmax1)
	parameter (km10=10*kmini)
	parameter (nmaxi=nmini*kmini)
cc
	integer rfncomb,rfnbreak
	integer ierr,matz,seed,tottimes,step
	integer pnsel
	integer flag(km10)
	integer mini(kmini)
	integer subdat(2,nmaxi)
	double precision mcdndex(10,2,kmini)
	integer subndex(450)
	integer replow
	integer fit
cc	double precision chi2(50)
cc	double precision chimed(50)
cc. consistency correction now happens in R code
cc.	double precision faclts(11)
	double precision pivot,rfmahad,medi2

	integer inbest(nhalff)
	integer weight(n)
	double precision coeff(kmini,nvar)
	double precision dat(n,nvar)
	double precision initcov(nvar*nvar)
	double precision adcov(nvar*nvar)
	double precision initmean(nvar)

	double precision med1,med2
	integer temp(n)
	integer index1(n)
	integer index2(n)
	double precision nmahad(n)
	double precision ndist(n)
	double precision am(n),am2(n),slutn(n)

	double precision med(nvar)
	double precision mad(nvar)
	double precision sd(nvar)
	double precision means(nvar)
	double precision bmeans(nvar)
	double precision w(nvar),fv1(nvar),fv2(nvar)

	double precision rec(nvar+1)
	double precision sscp1((nvar+1)*(nvar+1))
	double precision cova1(nvar*nvar)
	double precision corr1(nvar*nvar)
	double precision cinv1(nvar*nvar)
	double precision cova2(nvar*nvar)
	double precision cinv2(nvar*nvar)
	double precision z(nvar*nvar)


	double precision cstock(10,nvar*nvar)
	double precision mstock(10,nvar)
	double precision c1stock(km10,nvar*nvar)
	double precision m1stock(km10,nvar*nvar)
	double precision dath(nmaxi,nvar)

	double precision percen

	logical all,part,fine,final,rfodd,class

cc  Median of the chi-squared distribution:
cc	data chimed/0.454937,1.38629,2.36597,3.35670,4.35146,
cc     *  5.34812,6.34581,7.34412,8.34283,9.34182,10.34,11.34,12.34,
cc     *  13.34,14.34,15.34,16.34,17.34,18.34,19.34,20.34,21.34,22.34,
cc     *  23.34,24.34,25.34,26.34,27.34,28.34,29.34,30.34,31.34,32.34,
cc     *  33.34,34.34,35.34,36.34,37.34,38.34,39.34,40.34,41.34,42.34,
cc     *  43.34,44.34,45.34,46.34,47.33,48.33,49.33/
cc  The 0.975 quantile of the chi-squared distribution:
cc	data chi2/5.02389,7.37776,9.34840,11.1433,12.8325,
cc     *  14.4494,16.0128,17.5346,19.0228,20.4831,21.920,23.337,
cc     *  24.736,26.119,27.488,28.845,30.191,31.526,32.852,34.170,
cc     *  35.479,36.781,38.076,39.364,40.646,41.923,43.194,44.461,
cc     *  45.722,46.979,48.232,49.481,50.725,51.966,53.203,54.437,
cc     *  55.668,56.896,58.120,59.342,60.561,61.777,62.990,64.201,
cc     *  65.410,66.617,67.821,69.022,70.222,71.420/

cc. consistency correction now happens in R code
cc	 data faclts/2.6477,2.5092,2.3826,2.2662,2.1587,
cc.	*  2.0589,1.9660,1.879,1.7973,1.7203,1.6473/


C	CALL INTPR('Entering RFFASTMCD - KREP: ',-1,KREP,1)

        call rndstart
C            -------- == GetRNGstate() in C

C	20.06.2005 - substitute the parameters nmax and nvmax
	nmax = n
	nvmax = nvar

	nvmax1=nvmax+1
	nvmax2=nvmax*nvmax
	nvm12=nvmax1*nvmax1

	nrep = krep
	part=.false.
	fine=.false.
	final=.false.
	all=.true.
	kstep=k1
	medi2=0

c. These tests are superfluous, now that nmax == n,  nvmax == nvar :

c.	if(nvar.gt.nvmax) then
c. c 9400
c.	      fit= -2
c.	      kount = nvmax
c.	      goto 9999
c.	   endif
c.
c.	   if(n.gt.nmax) then
c. c 9200
c.	      fit= -1
c.	      kount = nmax
c.	      goto 9999
c.	   endif

cc
cc  From here on, the sample size n is known.
cc  Some initializations can be made. First of all, h (= the number of
cc  observations on which the MCD is based) is given by the integer variable
cc  nhalff.
cc  If nhalff equals n, the MCD is the classical covariance matrix.
cc  The logical value class indicates this situation.
cc  The variable nbreak is the breakdown point of the MCD estimator
cc  based on nhalff observations, whereas jdefaul is the optimal value of
cc  nhalff, with maximal breakdown point. The variable percen is the
cc  corresponding percentage.
cc
	percen = (1.D0*nhalff)/(1.D0*n)

	if(nvar.lt.5) then
	  eps=1.0D-12
	else
	  if(nvar.ge.5.and.nvar.le.8) then
	    eps=1.0D-14
	  else
	    eps=1.0D-16
	  endif
	endif

	jbreak=rfnbreak(nhalff,n,nvar)
	class= .false.
	if(nhalff.ge.n) then
c	 compute *only* the classical estimate
	  class= .true.
	  goto 9500
	endif

	if(nvar.eq.1) then
	  do 23, jj=1,n
 23	    ndist(jj)=dat(jj,1)
	  call rfshsort(ndist,n)
cc. consistency correction now happens in R code
cc.	  nquant=min(int(real(((nhalff*1.D0/n)-0.5D0)*40))+1,11)
cc.	  factor=faclts(nquant)
cc.	  call rfmcduni(ndist,n,nhalff,slutn,bstd,am,am2, factor,
	  call rfmcduni(ndist,n,nhalff,slutn,bstd,am,am2, 1.d0,
     *	    n-nhalff+1)
	  initmean(1)=slutn(1)
	  adcov(1)=bstd
	  initcov(1)=bstd
	  goto 9999
	endif
cc
cc  Some initializations:
cc    seed = starting value for random generator
cc    matz = auxiliary variable for the subroutine rs, indicating whether
cc	     or not eigenvectors are calculated
cc    nsel = number of variables + 1
cc    ngroup = number of subdatasets
cc    part = logical value, true if the dataset is split up
cc    fine = logical value, becomes true when the subsets are merged
cc    final = logical value, to indicate the final stage of the algorithm
cc    all = logical value, true for small n, if all (p+1)-subsets out of
cc	    n can be drawn
cc    subdat = matrix with a first row containing indices of observations
cc	       and a second row indicating the corresponding subdataset
cc
	seed=iseed
	matz=1
	nsel=nvar+1
	ngroup=1
	part=.false.
	fine=.false.
	final=.false.
	all=.true.
	do 21,i=1,nmaxi
	  subdat(1,i)=1000000
	  subdat(2,i)=1000000
 21	continue
cc
cc  Determine whether the dataset needs to be divided into subdatasets
cc  or can be treated as a whole. The subroutine rfrdraw constructs
cc  nonoverlapping subdatasets, with uniform distribution of the case numbers.
cc  For small n, the number of trial subsamples is determined.
cc
c MM(FIXME):  The following code depends crucially on  kmini == 5

	do 22 i=1,kmini
 22	   mini(i)=0
	if(n.gt.(2*nmini-1)) then
	  kstep=k1
	  part=.true.
	  ngroup=int(n/(nmini*1.D0))
	  if(n.ge.(2*nmini) .and. n.le.(3*nmini-1)) then
	    if(rfodd(n)) then
	      mini(1)=int(n/2)
	      mini(2)=int(n/2)+1
	    else
	      mini(1)=n/2
	      mini(2)=n/2
	    endif
	  else
	    if(n.ge.(3*nmini) .and. n.le.(4*nmini-1)) then
	      if(3*(n/3) .eq. n) then
		mini(1)=n/3
		mini(2)=n/3
		mini(3)=n/3
	      else
		mini(1)=int(n/3)
		mini(2)=int(n/3)+1
		if(3*(n/3) .eq. n-1) then
		  mini(3)=int(n/3)
		else
		  mini(3)=int(n/3)+1
		endif
	      endif
	    else
	      if(n.ge.(4*nmini) .and. n.le.(5*nmini-1)) then
		if(4*(n/4) .eq. n) then
		  mini(1)=n/4
		  mini(2)=n/4
		  mini(3)=n/4
		  mini(4)=n/4
		else
		  mini(1)=int(n/4)
		  mini(2)=int(n/4)+1
		  if(4*(n/4) .eq. n-1) then
		    mini(3)=int(n/4)
		    mini(4)=int(n/4)
		  else
		    if(4*(n/4) .eq. n-2) then
		      mini(3)=int(n/4)+1
		      mini(4)=int(n/4)
		    else
		      mini(3)=int(n/4)+1
		      mini(4)=int(n/4)+1
		    endif
		  endif
		endif
	      else
		mini(1)=nmini
		mini(2)=nmini
		mini(3)=nmini
		mini(4)=nmini
		mini(5)=nmini
	      endif
	    endif
	  endif
	  nhalf=int(mini(1)*percen)
	  if(ngroup.gt.kmini) ngroup=kmini
	  nrep=int((krep*1.D0)/ngroup)
	  minigr=mini(1)+mini(2)+mini(3)+mini(4)+mini(5)
	  call rfrdraw(subdat,n,seed,minigr,mini,ngroup,kmini)
	else
	  minigr=n
	  nhalf=nhalff
	  kstep=k1
	  if(n.le.replow(nsel)) then
	     nrep=rfncomb(nsel,n)
	  else
C VT::02.09.2004 - remove the hardcoded 500 for nrep
C	     nrep=500
	    nrep=krep
	    all=.false.
	  endif
	endif
	seed=iseed

cc
cc  Some more initializations:
cc    m1stock = matrix containing the means of the ngroup*10 best estimates
cc		obtained in the subdatasets.
cc    c1stock = matrix containing the covariance matrices of the ngroup*10
cc		best estimates obtained in the subdatasets.
cc    mstock = matrix containing the means of the ten best estimates
cc	       obtained after merging the subdatasets and iterating from
cc	       their best estimates.
cc    cstock = matrix containing the covariance matrices of the ten best
cc	       estimates obtained after merging the subdatasets
cc	       and iterating from their best estimates.
cc    means = mean vector
cc    bmeans = initial MCD location estimate
cc    sd = standard deviation vector
cc    nmahad = vector of mahalanobis distances
cc    ndist = vector of general (possibly robust) distances
cc    inbest = best solution vector
cc    index1 = index vector of subsample observations
cc    index2 = index vector of ordered mahalanobis distances
cc    temp  = auxiliary vector
cc    flag = vector with components indicating the occurrence of a
cc	     singular intermediate MCD estimate.
cc
	do 31 j=1,nvmax
	   do 33 k=1,10
	      mstock(k,j)=1000000.D0
	      do 35 kk=1,kmini
		 m1stock((kk-1)*10+k,j)=1000000.D0
 35	      continue
	      do 37 i=1,nvmax
		 do 39 kk=1,kmini
		    c1stock((kk-1)*10+k,(j-1)*nvmax+i)=1000000.D0
 39		 continue
		 cstock(k,(j-1)*nvmax+i)=1000000.D0
 37	      continue
 33	   continue
	   means(j)=0.D0
	   bmeans(j)=0.D0
	   sd(j)=0.D0
 31	continue

	do 41 j=1,nmax
	   nmahad(j)=0.D0
	   ndist(j)=0.D0
	   index1(j)=1000000
	   index2(j)=1000000
	   temp(j)=1000000
 41	continue
	do 43 j=1,km10
 43	   flag(j)=1


 9500	continue
cc
cc ********* Compute the classical estimates **************
cc
	call rfcovinit(sscp1,nvar+1,nvar+1)
	do 51 i=1,n
	  do 53 j=1,nvar
 53	    rec(j)=dat(i,j)
	  call rfadmit(rec,nvar,nvar+1,sscp1)
 51	continue
	call rfcovar(n,nvar,nvar+1,sscp1,cova1,means,sd)
	do 57 j=1,nvar
	  if(sd(j).eq.0.D0)	goto 5001
 57	continue

	call rfcovcopy(cova1,cinv1,nvar,nvar)
	det=1.D0
	do 58 j=1,nvar
	  pivot=cinv1((j-1)*nvar+j)
	  det=det*pivot
	  if(pivot.lt.eps)	goto 5001

	  call rfcovsweep(cinv1,nvar,j)
 58	continue
	call rfcorrel(nvar,cova1,corr1,sd)

c      if just classical estimate, we are done
	if(class)		goto 9999

	goto 5002

c  singularity '1' (exact fit := 1) :
 5001	continue
	call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
	call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
	call rfexact(kount,n,ndist, nvmax1,nvar,
     *	     sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
	call rfcovcopy(cova1,initcov,nvar,nvar)
	call rfcovcopy(means,initmean,nvar,1)
	do 56 j=1,nvar
 56	   coeff(1,j)=z(j)
	fit=1
	goto 9999


 5002	continue
cc
cc  Compute and store classical Mahalanobis distances.
cc
	do 62 j=1,n
	  do 64 i=1,nvar
 64	    rec(i)=dat(j,i)
	  nmahad(j)=rfmahad(rec,nvar,means,cinv1)
 62	continue


cc******** Compute the MCD estimates ************** ------------------------------

cc     Main loop: inspects the subsamples.
cc     Every time the sscp of the subsample is placed in sscp1,
cc     its covariance matrix in cova1, and its inverse in cinv1 .
cc     The minimum covariance determinant matrix is placed in cova2,
cc     and its inverse in cinv2.
cc     The robust distances are placed in ndist.
cc     In the main loop we count the total number of iteration steps
cc     with the variable tottimes.
cc
cc     The algorithm returns here twice when the dataset is divided
cc     at the beginning of the program. According to the situation,
cc     new initializations are made. The second stage, where the subdatasets
cc     are merged, is indicated by the logical value fine and
cc     the last stage, when the whole dataset is considered, by the logical
cc     variable final. In the last stage, the number of iterations nrep
cc     is determined according to the total number of observations
cc     and the dimension.
cc
	tottimes=0
 5555	object=10.D25

C	CALL INTPR('MAIN LOOP - NUMBER of TRIALS NREP: ',-1,NREP,1)

	if(.not. part .or. final) then
	  nn=n
	else
c	  (part .and. .not. final)
	   if (fine) then
	      nn=minigr
	   endif
	endif

	if(fine .or.(.not.part.and.final)) then
	  nrep=10
	  nsel=nhalf
	  kstep=k2
	  if (final) then
	    nhalf=nhalff
	    ngroup=1
	    if (n*nvar .le.100000) then
	      kstep=k3
	    else
	      if (n*nvar .gt.100000 .and. n*nvar .le.200000) then
		kstep=10
	      else
		if (n*nvar .gt.200000 .and. n*nvar
     *		  .le.300000) then
		  kstep=9
		else
		  if (n*nvar .gt.300000 .and. n*nvar
     *		    .le.400000) then
		    kstep=8
		  else
		    if (n*nvar .gt.400000 .and. n*nvar
     *		      .le.500000) then
		      kstep=7
		    else
		      if (n*nvar .gt.500000 .and. n*nvar
     *			.le.600000) then
			kstep=6
		      else
			if (n*nvar .gt.600000 .and. n*nvar
     *			  .le.700000) then
			  kstep=5
			else
			  if (n*nvar .gt.700000 .and. n*nvar
     *			    .le.800000) then
			    kstep=4
			  else
			    if (n*nvar .gt.800000 .and. n*nvar
     *			      .le.900000) then
			      kstep=3
			    else
			      if (n*nvar .gt.900000 .and. n*nvar
     *				.le.1000000) then
				kstep=2
			      else
				kstep =1
			      endif
			    endif
			  endif
			endif
		      endif
		    endif
		  endif
		endif
	      endif
	    endif
	    if (n.gt.5000) then
	      nrep=1
	    endif
	  else
	    nhalf=int(minigr*percen)
	  endif
	endif

	do 81 i=1,nsel-1
	  index1(i)=i
 81	continue
	index1(nsel)=nsel-1

cc
cc  Initialization of the matrices to store partial results. For the
cc  first stage of the algorithm, the currently best covariance matrices and
cc  means are stored in the matrices c1stock and m1stock initialized earlier.
cc  The corresponding objective values and the number of the trial subset
cc  are stored in the matrix mcdndex.
cc  For the second stage of the algorithm or for small datasets, only the
cc  currently best objective values are stored in the same matrix mcdndex
cc  and the corresponding covariance matrices and mean vectors are stored in
cc  the matrices cstock and mstock initialized earlier.
cc
	if(.not. final) then
	  do 83 i=1,10
	    do 85 j=1,ngroup
	      mcdndex(i,1,j)=10.D25
 85	      mcdndex(i,2,j)=10.D25
 83	   continue
	endif
	if(.not.fine.and..not.final) then
	  do 82 j=1,nvar
	    do 84 i=1,n
	      am(i)=dat(i,j)
 84	      am2(i)=dat(i,j)
	    if(2*n/2 .eq. n) then
	      med1=rffindq(am,n,n/2,index2)
	      med2=rffindq(am2,n,(n+2)/2,index2)
	      med(j)=(med1+med2)/2
	    else
	      med(j)=rffindq(am,n,(n+1)/2,index2)
	    endif
	    do 86 i=1,n
 86	      ndist(i)=dabs(dat(i,j)-med(j))
	    mad(j)=rffindq(ndist,n,nhalff,index2)
	    if(mad(j)-0.D0 .lt. eps) then
	      do 80,k=1,j-1
		do 79,i=1,n
 79		  dat(i,k)=dat(i,k)*mad(k)+med(k)
 80	      continue
	      call rfcovinit(sscp1,nvar+1,nvar+1)
	      do 88 k=1,nsel
		do 89 m=1,nvar
 89		  rec(m)=dat(index2(k),m)
		call rfadmit(rec,nvar,nvar+1,sscp1)
 88	      continue
	      call rfcovar(nsel,nvar,nvar+1,sscp1,cova1,means,sd)
	      call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
	      if(z(j).ne.1) then
		do 77, kk=1,nvar
		  if(z(kk*nvar+j).eq.1) then
		    do 75, l=1,nvar
 75		       z(l)=z(kk*nvar+l)
		    goto 76
		  endif
 77		continue
	      endif
 76	      continue
	      call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
	      call rfexact(kount,n,ndist, nvmax1,nvar,
     *		   sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
	      call rfcovcopy(cova1,initcov,nvar,nvar)
	      call rfcovcopy(means,initmean,nvar,1)
	      do 78 jjj=1,nvar
 78		coeff(1,jjj)=z(jjj)
	      fit=2
	      goto 9999
	    endif
	    do 87 i=1,n
 87	      dat(i,j)=(dat(i,j)-med(j))/mad(j)
 82	  continue
	endif
cc
cc  The matrix dath contains the observations to be used in the
cc  algorithm. In the first stage of the split-up procedure dath contains
cc  nmini objects, corresponding to the original observations, with the index
cc  of the processed group in the array subdat. For the second stage, the
cc  data points of all the subdatasets are merged in dath.
cc  The variable kount indicates the occurrence of a singular subsample leading
cc  to the corresponding plane. In some situations the variable kount counts
cc  the number of observations on that plane.
cc
	if (fine .and. .not. final) then
	  do 91, j=1,minigr
	    do 93, k=1,nvar
 93	      dath(j,k)=dat(subdat(1,j),k)
 91	  continue
	endif
	kount=0

c---- For-Loop over groups  - - - - - - - - - - - - - - - - - - - - -
	do 1111 ii= 1,ngroup
	  if(.not.fine) kount=0
	  if(part .and. .not. fine) nn=mini(ii)
	  do 101 i=1,nn
	    index2(i)=i
 101	  continue
	  if(part .and. .not. fine) then
	    jndex=0
	    do 103 j=1,minigr
	      if(subdat(2,j).eq.ii) then
		jndex=jndex+1
		subndex(jndex)=subdat(1,j)
	      endif
 103	    continue
	    do 105 j=1,mini(ii)
	      do 107 k=1,nvar
 107		dath(j,k)=dat(subndex(j),k)
 105	    continue
	  endif

cc  The number of trial subsamples is represented by nrep, which depends
cc  on the data situation.
cc  When all (p+1)-subsets out of n can be drawn, the subroutine rfgenpn
cc  is used. Otherwise, random subsamples are drawn by the routine
cc  rfrangen. The trial subsamples are put in the array index1. The
cc  same thing happens for large datasets, except that the number of
cc  observations is nmini instead of n.
cc
cc  When a trial subsample is singular, the algorithm counts the number of
cc  observations that lie on the hyperplane corresponding to this sample.
cc  If, for small datasets, this number is larger than nhalff, the program
cc  stops (exact fit) and gives the mean and the covariance matrix
cc  of the observations on the hyperplane, together with the equation
cc  of the hyperplane.
cc  For large datasets, the algorithm first checks whether there are more
cc  than nhalff observations on the hyperplane. If this is the case, the
cc  program stops for the same reason of exact fit and gives the covariance
cc  matrix and mean of the observations on the hyperplane. If not, the
cc  algorithm counts the number of observations that lie on the hyperplane.
cc  When this number is smaller than the current nhalf in the subdataset, these
cc  observations are extended to nhalf observations by adding those
cc  observations that have smallest orthogonal distances to the hyperplane
cc  and the algorithm continues.
cc  When larger, the coefficients of the hyperplane are stored in the matrix
cc  m1stock for use as starting value in the next stage, and the flag of this
cc  estimate gets the value zero.
cc
cc  In the second stage of the algorithm, when the subdatasets are merged,
cc  the array index2 contains the indices of the observations
cc  corresponding to the nhalf observations with minimal relative distances
cc  with respect to the best estimates of the first stage.
cc  When the estimate of the first stage is a hyperplane, the algorithm
cc  investigates whether there are more than the current nhalf observations of
cc  the merged subdataset on that hyperplane. If so, the coefficients of the
cc  hyperplane are again stored, now in the matrix mstock, for the final
cc  stage of the algorithm.
cc  If not, the observations on the hyperplane are extended to nhalf
cc  observations by adding the observations in the merged dataset with
cc  smallest orthogonal distances to that hyperplane.
cc  For small datasets or for larger datasets with n <= nmini*kmini,
cc  the algorithm already stops when one solution becomes singular,
cc  since we then have an exact fit.
cc
cc  In the third stage, the covariance matrices and means of the best
cc  solutions of the second stage are used as starting values.
cc  Again, when a solution becomes singular, the subroutine 'exact'
cc  determines the hyperplane through at least nhalff observations and stops
cc  because of the exact fit.
cc
cc  When the program stops because of an exact fit, the covariance matrix and
cc  mean of the observations on the hyperplane will always be given.
cc
	do 1000 i=1,nrep
	  pnsel=nsel
	  tottimes=tottimes+1
	  deti=0.D0
	  detimin1=0.D0
	  step=0
	  call rfcovinit(sscp1,nvar+1,nvar+1)
	  if((part.and..not.fine).or.(.not.part.and..not.final)) then
	    if(part) then
		call rfrangen(mini(ii),nsel,index1,seed)
	    else
	      if(all) then
		call rfgenpn(n,nsel,index1)
	      else
		call rfrangen(n,nsel,index1,seed)
	      endif
	    endif
	  endif
cc
cc  The covariance matrix and mean of the initial subsamples are
cc  calculated with the subroutine covar and represented by
cc  the variables cova1 and means.
cc
cc  In the following stages of the algorithm, the covariance matrices and means
cc  used as starting values are already stored in the matrices c1stock
cc  and m1stock (for the second stage), and in the matrices cstock and mstock
cc  (for the third stage).
cc
cc  The inverse cinv1 of the covariance matrix is calculated by the
cc  subroutine rfcovsweep, together with its determinant det.
cc
 9550	  call rfcovinit(sscp1,nvar+1,nvar+1)
	  if(.not.fine.and.part) then
	    do 121 j=1,pnsel
	      do 123 m=1,nvar
 123		rec(m)=dath(index1(j),m)
	      call rfadmit(rec,nvar,nvar+1,sscp1)
 121	    continue
	    call rfcovar(pnsel,nvar,nvar+1,sscp1,cova1,means,sd)
	  endif
	  if(.not.part.and..not.final) then
	    do 122 j=1,pnsel
	      do 124 m=1,nvar
 124		rec(m)=dat(index1(j),m)
	      call rfadmit(rec,nvar,nvar+1,sscp1)
 122	    continue
	    call rfcovar(pnsel,nvar,nvar+1,sscp1,cova1,means,sd)
	  endif
	  if (final) then
	    if(mstock(i,1).ne.1000000.D0) then
	      do 125 jj=1,nvar
		means(jj)=mstock(i,jj)
		do 127 kk=1,nvar
		  cova1((jj-1)*nvar+kk)=cstock(i,(jj-1)*nvar+kk)
 127		continue
 125	      continue
	    else
	      goto 1111
	    endif
	    if(flag(i).eq.0) then
	      qorder=1.D0
	      do 129,jjj=1,nvar
 129		z(jjj)=coeff(1,jjj)
	      call rfdis(dat,z,ndist,n,nvar,nn,nvar, means)
	      dist2=rffindq(ndist,nn,nhalf,index2)
	      goto 9555
	    endif
	  endif
	  if (fine .and. .not.final) then
	      if(m1stock((ii-1)*10+i,1).ne.1000000.D0) then
		do 131 jj=1,nvar
		  means(jj)=m1stock((ii-1)*10+i,jj)
		  do 133 kk=1,nvar
		    cova1((jj-1)*nvar+kk)=c1stock((ii-1)*10+i,
     *			  (jj-1)*nvar+kk)
 133		 continue
 131		continue
	      else
		goto 1111
	      endif
	    if(flag((ii-1)*10+i).eq.0) then
	      qorder=1.D0
	      do 135,jjj=1,nvar
 135		z(jjj)=coeff(ii,jjj)
	      call rfdis(dath,z,ndist,nmaxi,nvmax,nn,nvar, means)
	      call rfshsort(ndist,nn)
	      qorder=ndist(nhalf)
	      if(dabs(qorder-0.D0).lt.10.D-8.and.kount.eq.0
     *		   .and.n.gt.nmini*kmini) then
		kount=nhalf
		do 137,kkk=nhalf+1,nn
		  if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
		    kount=kount+1
		  endif
 137		continue
		flag(1)=0
		do 139,kkk=1,nvar
 139		  coeff(1,kkk)=z(kkk)
		call rfstore2(nvar,cstock,mstock,nvmax2,nvmax,
     *		       kmini,cova1,means,i,mcdndex,kount)
		kount=1
		goto 1000
	      else
		if(dabs(qorder-0.D0).lt.10.D-8.and.
     *		      kount.ne.0.and.n.gt.nmini*kmini) then
		  goto 1000
		else
		  flag(1)=1
		  dist2=rffindq(ndist,nn,nhalf,index2)
		  goto 9555
		endif
	      endif
	    endif
	  endif
	  call rfcovcopy(cova1,cinv1,nvar,nvar)
	  det=1.D0
	  do 200 j=1,nvar
	    pivot=cinv1((j-1)*nvar+j)
	    det=det*pivot
	    if(pivot.lt.eps) then
	      call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
	      qorder=1.D0
	      if(.not.part.or.final) then
		call rfdis(dat,z,ndist,n,nvar,nn,nvar,means)
	      else
		call rfdis(dath,z,ndist,nmaxi,nvmax,nn,nvar,means)
	      endif
	      call rfshsort(ndist,nn)
	      qorder=ndist(nhalf)
	      if(dabs(qorder-0.D0).lt. 10.D-8 .and. .not.part) then
		call transfo(cova1,means,dat,med,mad,nvar,n)
		call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		call rfdis(dat,z,ndist,n,nvar,nn,nvar,means)
		call rfexact(kount,n,ndist, nvmax1,nvar,
     *		     sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
		call rfcovcopy(cova1,initcov,nvar,nvar)
		call rfcovcopy(means,initmean,nvar,1)
		do 140,jjj=1,nvar
 140		   coeff(1,jjj)=z(jjj)
		fit=2
		goto 9999
	      else
		if(dabs(qorder-0.D0).lt. 10.D-8 .and. part .and.
     *		      kount.eq.0) then
		  call rfdis(dat,z,ndist,n,nvar,n,nvar, means)
		  call rfshsort(ndist,n)
		  if(dabs(ndist(nhalff)-0.D0).lt.10.D-8) then
		    call transfo(cova1,means,dat,med,mad,nvar,n)
		    call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		    call rfdis(dat,z,ndist,n,nvar,nn,nvar,means)
		    call rfexact(kount,n,ndist, nvmax1,nvar,
     *			sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
		    call rfcovcopy(cova1,initcov,nvar,nvar)
		    call rfcovcopy(means,initmean,nvar,1)
		    do 142,jjj=1,nvar
 142		      coeff(1,jjj)=z(jjj)
		    fit=2
		    goto 9999
		  endif
		  call rfdis(dath,z,ndist,nmaxi,nvmax,nn,nvar, means)
		  call rfshsort(ndist,nn)
		  kount=nhalf
		  do 141,kkk=nhalf+1,nn
		    if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
		      kount=kount+1
		    endif
 141		  continue
		  flag((ii-1)*10+1)=0
		  do 143,kkk=1,nvar
 143		    coeff(ii,kkk)=z(kkk)
		  call rfstore1(nvar,c1stock,m1stock,nvmax2,nvmax,
     *			 kmini,cova1,means,i,km10,ii,mcdndex, kount)
		  kount=1
		  goto 1000
		else
		  if(dabs(qorder-0.D0).lt. 10.D-8 .and. part .and.
     *			kount.ne.0) then
		    goto 1000
		  else
		    call rfishsort(index1,pnsel)
		    call prdraw(index1,pnsel, nn)
		    pnsel=pnsel+1
		    goto 9550
		  endif
		endif
	      endif
	    endif
	    call rfcovsweep(cinv1,nvar,j)
 200	  continue
cc
cc  Mahalanobis distances are computed with the subroutine rfmahad
cc  and stored in the array ndist.
cc  The k-th order statistic of the mahalanobis distances is stored
cc  in dist2. The array index2 containes the indices of the
cc  corresponding observations.
cc
	  do 151 j=1,nn
	    if(.not.part.or.final) then
	    do 152 mm=1,nvar
 152	      rec(mm)=dat(j,mm)
	    else
	    do 153 mm=1,nvar
 153	      rec(mm)=dath(j,mm)
	    endif
	    t=rfmahad(rec,nvar,means,cinv1)
	    ndist(j)=t
 151	  continue
	  dist2=rffindq(ndist,nn,nhalf,index2)
cc
cc  The variable kstep represents the number of iterations. They depend on
cc  the situation of the program (k1, k2, or k3). Within each
cc  iteration the mean and covariance matrix of nhalf observations are
cc  calculated. The nhalf smallest corresponding mahalanobis distances
cc  determine the subset for the next iteration.
cc  The best subset for the whole data is stored in the array inbest.
cc  The iteration stops when two subsequent determinants become equal.
cc
 9555	  do 400 step=1,kstep
	    tottimes=tottimes+1
	    call rfcovinit(sscp1,nvar+1,nvar+1)
	    do 155 j=1,nhalf
 155	      temp(j)=index2(j)
	    call rfishsort(temp,nhalf)
	    do 157 j=1,nhalf
	      if(.not.part.or.final) then
		do 158 mm=1,nvar
 158		   rec(mm)=dat(temp(j),mm)
	      else
		do 159 mm=1,nvar
 159		   rec(mm)=dath(temp(j),mm)
	      endif
	      call rfadmit(rec,nvar,nvar+1,sscp1)
 157	    continue
	    call rfcovar(nhalf,nvar,nvar+1,sscp1,cova1,means,sd)
	    call rfcovcopy(cova1,cinv1,nvar,nvar)
	    det=1.D0
	    do 600 j=1,nvar
	      pivot=cinv1((j-1)*nvar+j)
	      det=det*pivot
	      if(pivot.lt.eps) then
		if(final .or. .not.part .or.
     *		   (fine.and. .not.final .and. n .le. nmini*kmini)) then
		  call transfo(cova1,means,dat,med,mad,nvar,n)
		  call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		  if(final.or..not.part) then
		    call rfdis(dath,z,ndist,nmax,nvmax,nn,nvar, means)
		  else
		    call rfdis(dath,z,ndist,nmaxi,nvmax,nn,nvar, means)
		  endif
		  call rfexact(kount,n,ndist, nvmax1,nvar,
     *		       sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
		  call rfcovcopy(cova1,initcov,nvar,nvar)
		  call rfcovcopy(means,initmean,nvar,1)
		  do 160 jjj=1,nvar
 160		     coeff(1,jjj)=z(jjj)
		  fit=2
		  goto 9999
		endif
		if(part.and..not.fine.and.kount.eq.0) then
		  call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		  call rfdis(dat,z,ndist,n,nvar,n,nvar, means)
		  call rfshsort(ndist,n)
		  if(dabs(ndist(nhalff)-0.D0).lt.10.D-8) then
		    call transfo(cova1,means,dat,med,mad,nvar,n)
		    call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		    call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
		    call rfexact(kount,n,ndist, nvmax1,nvar,
     *			sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
		    call rfcovcopy(cova1,initcov,nvar,nvar)
		    call rfcovcopy(means,initmean,nvar,1)
		    do 161 jjj=1,nvar
 161		       coeff(1,jjj)=z(jjj)
		    fit=2
		    goto 9999
		  endif
		  call rfdis(dath,z,ndist,nmaxi,nvmax,nn,nvar, means)
		  call rfshsort(ndist,nn)
		  kount=nhalf
		  do 162,kkk=nhalf+1,nn
		    if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
		      kount=kount+1
		    endif
 162		  continue
		  flag((ii-1)*10+1)=0
		  do 164 kkk=1,nvar
 164		    coeff(ii,kkk)=z(kkk)
		  call rfstore1(nvar,c1stock,m1stock,nvmax2,nvmax,
     *			 kmini,cova1,means,i,km10,ii,mcdndex, kount)
		  kount=1
		  goto 1000
		else
		  if(part.and..not.fine.and.kount.ne.0) then
		    goto 1000
		  endif
		endif
		if(fine.and..not.final.and.kount.eq.0) then
		  call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		  call rfdis(dat,z,ndist,n,nvar,n,nvar, means)
		  call rfshsort(ndist,n)
		  if(dabs(ndist(nhalff)-0.D0).lt.10.D-8) then
		    call transfo(cova1,means,dat,med,mad,nvar,n)
		    call rs(nvar,nvar,cova1,w,matz,z,fv1,fv2,ierr)
		    call rfdis(dat,z,ndist,n,nvar,n,nvar,means)
		    call rfexact(kount,n,ndist, nvmax1,nvar,
     *			sscp1,rec,dat, cova1,means,sd,nvar+1,weight)
		    call rfcovcopy(cova1,initcov,nvar,nvar)
		    call rfcovcopy(means,initmean,nvar,1)
		    do 165 jjj=1,nvar
 165		       coeff(1,jjj)=z(jjj)
		    fit=2
		    goto 9999
		  endif
		  call rfdis(dath,z,ndist,nmaxi,nvmax,nn,nvar, means)
		  call rfshsort(ndist,nn)
		  kount=nhalf
		  do 166,kkk=nhalf+1,nn
		    if(dabs(ndist(kkk)-0.D0).lt.10.D-8) then
		      kount=kount+1
		    endif
 166		  continue
		  flag(1)=0
		  do 168,kkk=1,nvar
 168		    coeff(1,kkk)=z(kkk)
		  call rfstore2(nvar,cstock,mstock,nvmax2,nvmax,
     *			 kmini,cova1,means,i,mcdndex,kount)
		  kount=1
		  goto 1000
		else
		  if(fine.and..not.final.and.kount.ne.0) then
		    goto 1000
		  endif
		endif
	      endif
	      call rfcovsweep(cinv1,nvar,j)
 600	    continue
	    if(step.ge.2 .and. det.eq.detimin1) then
	      goto 5000
	    endif
	    detimin1=deti
	    deti=det
	    do 171 j=1,nn
	      if(.not.part.or.final) then
		do 172 mm=1,nvar
 172		  rec(mm)=dat(j,mm)
	      else
		do 173 mm=1,nvar
 173		  rec(mm)=dath(j,mm)
	      endif
	      t=rfmahad(rec,nvar,means,cinv1)
	      ndist(j)=t
 171	    continue
	    dist2=rffindq(ndist,nn,nhalf,index2)
	    dist=dsqrt(dist2)
	    if(((i.eq.1.and.step.eq.1.and..not.fine)
     *		 .or.det.lt.object).and.(final)) then
	      medi2=rffindq(ndist,nn,int(n/2),index1)
	      object=det
	      do 175 jjj=1,nhalf
		inbest(jjj)=index2(jjj)
 175	      continue
	      call rfcovcopy(cova1,cova2,nvar,nvar)
	      call rfcovcopy(cinv1,cinv2,nvar,nvar)
	      call rfcovcopy(means,bmeans,nvar,1)
	    endif
 400	  continue

cc  After each iteration, it has to be checked whether the new solution
cc  is better than some previous one and therefore needs to be stored. This
cc  isn't necessary in the third stage of the algorithm, where only the best
cc  solution is kept.

 5000	  if(.not. final) then
	    if(part .and. .not. fine) then
	      iii=ii
	    else
	      iii=1
cc	      At the end of the algorithm, only the ten
cc	      best solutions need to be stored.
	    endif

cc  For each data group :
cc    If the objective function is lower than the largest value in the
cc    matrix mcdndex :
cc    A distinction is made between different stages of the algorithm:
cc	* At the first stage of the split-up situation:
cc	  -If the new objective value did not yet occur in mcdndex
cc	   its value and corresponding covariance matrix and mean are
cc	   stored at the right place in the matrices mcdndex, c1stock and
cc	   m1stock, and other values are shifted to their new position
cc	   in these arrays.
cc	  -If the new objective value already occurs in mcdndex, a
cc	   comparison is made between the new mean vector and covariance matrix
cc	   and those estimates with the same determinant.
cc	   When for an equal determinant, the mean vector or covariance matrix
cc	   do not correspond, both of them are kept in the matrices mcdndex
cc	   and nbest.
cc	* In the second stage of the algorithm, the covariances and means
cc	  are stored :
cc	  - If the new objective value did not yet occur
cc	    in the matrix mcdndex, it is inserted by shifting the greater
cc	    determinants upwards and doing the same in the arrays mstock
cc	    and cstock.
cc	  - If the new objective value already occurs in the array mcdndex,
cc	    it is compared with all solutions with the same determinant.
cc	    In the case of an equality, the means and covariances
cc	    are compared to determine whether or not to insert the
cc	    new solution.
cc    Otherwise nothing happens. When a singularity occurs,
cc    the determinant in the matrix mcdndex is zero and the
cc    corresponding flag is zero too, so the search in the arrays mcdndex,
cc    m1stock, c1stock, mstock and cstock is done on the rows with flag one.
cc

	    if( flag((iii-1)*10+1).eq.1) then
	      lll=1
	    else
	      lll=2
	    endif
	    do 201, j=lll,10
	      if (det .le. mcdndex(j,2,iii)) then
		if(det.ne.mcdndex(j,2,iii)) then
		  if(.not.fine.and.part)			goto 203
								goto 205
		else
		  do 207 kkk=j,10
		    if(det.eq.mcdndex(kkk,2,iii)) then
		      do 209, jjj=1,nvar
			if(part.and..not.fine) then
			  if(means(jjj).ne.m1stock((iii-1)*10+ kkk,jjj))
     *								goto 203
			else
			  if(means(jjj).ne.mstock(kkk,jjj))	goto 205
			endif
 209		      continue
		      do 211, jjj=1,nvar*nvar
			if(part.and..not.fine) then
			  if(cova1(jjj).ne.c1stock((iii-1)*10+ kkk,jjj))
     *								goto 203
			else
			  if(cova1(jjj).ne.cstock(kkk,jjj))	goto 205
			endif
 211		      continue
		    endif
 207		  continue
		endif
		goto 1000

 203		do 221 k=10,j+1,-1
		  do 223 kk=1,nvar*nvar
 223		    c1stock((iii-1)*10+k,kk)=
     *		    c1stock((iii-1)*10+k-1,kk)
		  do 225 kk=1,nvar
 225		    m1stock((iii-1)*10+k,kk)=
     *		    m1stock((iii-1)*10+k-1,kk)
		  mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
		  mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
 221		continue
		do 227 kk=1,nvar
		  do 229 kkk=1,nvar
 229		     c1stock((iii-1)*10+j,(kk-1)*nvar+kkk)=
     *				    cova1((kk-1)*nvar+kkk)
		  m1stock((iii-1)*10+j,kk)=means(kk)
 227		continue
		mcdndex(j,1,iii)=i
		mcdndex(j,2,iii)=det
		goto 1000

 205		do 231 k=10,j+1,-1
		  do 233 kk=1,nvar*nvar
 233		    cstock(k,kk)=
     *		    cstock(k-1,kk)
		  do 235 kk=1,nvar
 235		    mstock(k,kk)=
     *		    mstock(k-1,kk)
		  mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
		  mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
 231		continue
		do 237 kk=1,nvar
		  do 239 kkk=1,nvar
 239		     cstock(j,(kk-1)*nvar+kkk)=
     *			cova1((kk-1)*nvar+kkk)
		  mstock(j,kk)=means(kk)
 237		continue
		mcdndex(j,1,iii)=i
		mcdndex(j,2,iii)=det
		goto 1000

	      endif
 201	    continue

	  endif
c		(not final)

 1000	continue
 1111	continue
c---- - - - - - end [ For-loop -- ii= 1, ngroup	 ]  - - - - - - - - -

cc  Determine whether the algorithm needs to be run again or not.
cc
	if(part .and. .not. fine) then
	  fine= .true.
	  goto 5555
	else if(.not. final .and. ((part.and.fine).or. .not.part)) then
	  final= .true.
	  goto 5555
	endif

cc******** end { Main Loop } ************** --------------------------------


	do 261, j=1,nhalf
	  temp(j)=inbest(j)
 261	continue
	call rfishsort(temp,nhalf)

C	CALL INTPR('BEST SUBSAMPLE (SORTED): ',-1,TEMP,NHALF)

	do 271,j=1,nvar
 271	  means(j)=bmeans(j)*mad(j)+med(j)
	call rfcovcopy(means,initmean,nvar,1)

C	CALL DBLEPR('CENTER: ',-1,initmean,nvar)
cc
	do 9145 i=1,nvar
	  do 9147 j=1,nvar
 9147	    cova1((i-1)*nvar+j)=cova2((i-1)*nvar+j)*mad(i)*mad(j)
 9145	continue
	call rfcovcopy(cova1,initcov,nvar,nvar)
	det=object
	do 9149 j=1,nvar
	  det=det*mad(j)*mad(j)
 9149	continue
cc

cc	VT::chimed is passed now as a parameter
cc	  call rfcovmult(cova1,nvar,nvar,medi2/chimed(nvar))
cc	  call rfcovmult(cova2,nvar,nvar,medi2/chimed(nvar))
cc	  call rfcovmult(cinv2,nvar,nvar,1.D0/(medi2/chimed(nvar)))

	call rfcovmult(cova1,nvar,nvar,medi2/chimed)
	call rfcovmult(cova2,nvar,nvar,medi2/chimed)
	call rfcovmult(cinv2,nvar,nvar,1.D0/(medi2/chimed))
	call rfcovcopy(cova1,adcov,nvar,nvar)
cc
cc	The MCD location is in bmeans.
cc	The MCD scatter matrix is in cova2,
cc	and its inverse in cinv2.
cc
cc	For every observation we compute its MCD distance
cc	and compare it to a cutoff value.
cc
	call rfcovinit(sscp1,nvar+1,nvar+1)
	nin=0

cc VT:: no need - the cutoff now is passed as a parameter
cc	cutoff=chi2(nvar)

	do 280 i=1,n
	  do 282 mm=1,nvar
 282	  rec(mm)=dat(i,mm)
	  dist2=rfmahad(rec,nvar,bmeans,cinv2)
	  if(dist2.le.cutoff) then
	    nin=nin+1
	    weight(i)=1
	  else
	    weight(i)=0
	  endif
 280	continue
cc
	call transfo(cova2,bmeans,dat,med,mad,nvar,n)
	goto 9999
cc ******************************************************************

 9999   continue
        call rndend
C            ------ == PutRNGstate() in C
	return
	end
ccccc   end {rffastmcd}
ccccc
ccccc
ccccc
ccccc
      subroutine rfexact(kount,nn,ndist, nvmax1,nvar,sscp1,
     *	rec,dat, cova1,means,sd,nvar1,weight)
cc
cc Determines how many objects lie on the hyperplane with equation
cc z(1,1)*(x_i1 - means_1)+ ... + z(p,1)* (x_ip - means_p) = 0
cc and computes their mean and their covariance matrix.
cc
      double precision ndist(nn)
      double precision sscp1(nvar1,nvar1)
      double precision rec(nvmax1)
      double precision dat(nn,nvar)
      double precision cova1(nvar,nvar)
      double precision means(nvar)
      double precision sd(nvar)
      integer weight(nn)
cc
      call rfcovinit(sscp1,nvar+1,nvar+1)
      kount=0
      do 10,kk=1,nn
	if(dabs(ndist(kk)-0.D0).lt.10.D-8) then
	  kount=kount+1
	  weight(kk)=1
	  do 20,j=1,nvar
 20	    rec(j)=dat(kk,j)
	  call rfadmit(rec,nvar,nvar+1,sscp1)
	else
	  weight(kk)=0
	endif
 10   continue
      call rfcovar(kount,nvar,nvar+1,sscp1,cova1,means,sd)
      return
      end
ccccc
ccccc
      subroutine transfo(cova,means,dat,med,mad,nvar,n)
cc
      double precision cova(nvar,nvar)
      double precision means(nvar)
      double precision dat(n,nvar)
      double precision med(nvar),mad(nvar)
cc
      do 5,j=1,nvar
	means(j)=means(j)*mad(j)+med(j)
	do 10, k=1,nvar
 10	  cova(j,k)=cova(j,k)*mad(j)*mad(k)
	do 20,i=1,n
 20	  dat(i,j)=dat(i,j)*mad(j)+med(j)
 5    continue
      return
      end
ccccc
ccccc
      subroutine rfcovmult(a,n1,n2,fac)
cc
cc  Multiplies the matrix a by the real factor fac.
cc
      double precision a(n1,n2)
      double precision fac
cc
      do 100 i=1,n1
	do 90 j=1,n2
	  a(i,j)=a(i,j)*fac
 90	continue
 100  continue
      return
      end
ccccc
ccccc
      subroutine rfadmit(rec,nvar,nvar1,sscp)
cc
cc  Updates the sscp matrix with the additional case rec.
cc
      double precision rec(nvar)
      double precision sscp(nvar1,nvar1)
cc
      sscp(1,1)=sscp(1,1)+1.D0
      do 10 j=1,nvar
	sscp(1,j+1)=sscp(1,j+1)+rec(j)
	sscp(j+1,1)=sscp(1,j+1)
 10   continue
      do 100 i=1,nvar
	do 90 j=1,nvar
	  sscp(i+1,j+1)=sscp(i+1,j+1)+rec(i)*rec(j)
 90	continue
 100  continue
      return
      end
ccccc
ccccc
      subroutine rfcovar(n,nvar,nvar1,sscp,cova,means,sd)
cc
cc  Computes the classical mean and covariance matrix.
cc
      double precision sscp(nvar1,nvar1)
      double precision cova(nvar,nvar)
      double precision means(nvar)
      double precision sd(nvar)
      double precision f
cc
      do 100 i=1,nvar
	means(i)=sscp(1,i+1)
	sd(i)=sscp(i+1,i+1)
	f=(sd(i)-means(i)*means(i)/n)/(n-1)
	if(f.gt.0.D0) then
	  sd(i)=dsqrt(f)
	else
	  sd(i)=0.D0
	endif
	means(i)=means(i)/n
 100  continue
      do 200 i=1,nvar
	do 190 j=1,nvar
	  cova(i,j)=sscp(i+1,j+1)
 190	continue
 200  continue
      do 300 i=1,nvar
	do 290 j=1,nvar
	  cova(i,j)=cova(i,j)-n*means(i)*means(j)
	  cova(i,j)=cova(i,j)/(n-1)
 290	continue
 300  continue
      return
      end
ccccc
ccccc
	subroutine rfcorrel(nvar,a,b,sd)
cc
cc  Transforms the scatter matrix a to the correlation matrix b.
cc
	double precision a(nvar,nvar)
	double precision b(nvar,nvar)
	double precision sd(nvar)
cc
	do 10,j=1,nvar
 10	  sd(j)=1/sqrt(a(j,j))
	do 100 i=1,nvar
	  do 90 j=1,nvar
	    if(i.eq.j) then
	      b(i,j)=1.0
	    else
	      b(i,j)=a(i,j)*sd(i)*sd(j)
	    endif
 90	  continue
 100	continue
	return
	end

	subroutine prdraw(a,pnsel, nn)
cc
        implicit none
	integer nn, a(nn), pnsel
c
        double precision unifrnd
	integer jndex, nrand, i,j
cc
	jndex=pnsel
cOLD 	nrand=int(uniran(seed)*(nn-jndex))+1
	nrand=int(unifrnd() * (nn-jndex))+1
C         if(nrand .gt. nn-jndex) then
C            call intpr(
C      1          '** prdraw(): correcting nrand > nn-jndex; nrand=',
C      2          -1, nrand, 1)
C            nrand=nn-jndex
C         endif

	jndex=jndex+1
	a(jndex)=nrand+jndex-1
	do 5, i=1,jndex-1
           if(a(i).gt.nrand+i-1) then
              do 6,j=jndex,i+1,-1
                 a(j)=a(j-1)
 6            continue
              a(i)=nrand+i-1
              goto 10
c             ------- break
           endif
 5	continue
 10	continue
	return
	end
ccccc
ccccc
	function rfmahad(rec,nvar,means,sigma)
cc
cc  Computes a Mahalanobis-type distance.
cc
	double precision rec(nvar), means(nvar), sigma(nvar,nvar)
	double precision rfmahad, t
cc
	t=0
	do 100 j=1,nvar
	  do 90 k=1,nvar
	    t=t+(rec(j)-means(j))*(rec(k)-means(k))*sigma(j,k)
 90	  continue
 100	continue
	rfmahad=t
	return
	end
ccccc
ccccc
       subroutine rfdis(da,z,ndist,nm,nv,nn,nvar, means)
cc
cc Computes the distance between the objects of da and a hyperplane with
cc equation z(1,1)*(x_i1 - means_1) + ... + z(p,1)*(x_ip - means_p) = 0
cc
       double precision da(nm,nv)
       double precision z(nvar,nvar)
       double precision ndist(nn)
       double precision means(nvar)
cc
       do 10, i=1,nn
	 ndist(i)=0
	 do 20, j=1,nvar
 20	   ndist(i)=z(j,1)*(da(i,j)-means(j))+ndist(i)
 10	   ndist(i)=dabs(ndist(i))
       return
       end
ccccc
ccccc
      subroutine rfstore2(nvar,cstock,mstock,nvmax2,nvmax,
     * kmini,cova1,means,i,mcdndex,kount)
cc
cc  Stores the coefficients of a hyperplane
cc  z(1,1)*(x_i1 - means_1) + ... +  z(p,1)*(x_ip - means_p) = 0
cc  into the first row of the matrix mstock, and shifts the other
cc  elements of the arrays mstock and cstock.
cc
      double precision cstock(10,nvmax2)
      double precision mstock(10,nvmax)
      double precision mcdndex(10,2,kmini)
      double precision cova1(nvar,nvar)
      double precision means(nvar)
cc
      do 10,k=10,2,-1
	do 20 kk=1,nvar*nvar
 20	  cstock(k,kk)=
     *	cstock(k-1,kk)
	do 30 kk=1,nvar
 30	  mstock(k,kk)=
     *	mstock(k-1,kk)
	mcdndex(k,1,1)=mcdndex(k-1,1,1)
	mcdndex(k,2,1)=mcdndex(k-1,2,1)
 10   continue
      do 40 kk=1,nvar
	mstock(1,kk)=means(kk)
	do 50 jj=1,nvar
 50	  cstock(1,(kk-1)*nvar+jj)=cova1(kk,jj)
 40   continue
      mcdndex(1,1,1)=i
      mcdndex(1,2,1)=kount
      return
      end
ccccc
ccccc
      subroutine rfstore1(nvar,c1stock,m1stock,nvmax2,nvmax,
     * kmini,cova1,means,i,km10,ii,mcdndex,kount)
cc
      double precision c1stock(km10,nvmax2)
      double precision m1stock(km10,nvmax)
      double precision mcdndex(10,2,kmini)
      double precision cova1(nvar,nvar)
      double precision means(nvar)
cc
      do 10,k=10,2,-1
	do 20 kk=1,nvar*nvar
 20	  c1stock((ii-1)*10+k,kk)=
     *	c1stock((ii-1)*10+k-1,kk)
	do 30 kk=1,nvar
 30	  m1stock((ii-1)*10+k,kk)=
     *	m1stock((ii-1)*10+k-1,kk)
	mcdndex(k,1,ii)=mcdndex(k-1,1,ii)
	mcdndex(k,2,ii)=mcdndex(k-1,2,ii)
 10   continue
      do 40 kk=1,nvar
	m1stock((ii-1)*10+1,kk)=means(kk)
	do 50 jj=1,nvar
 50	  c1stock((ii-1)*10+1,(kk-1)*nvar+jj)=
     *	  cova1(kk,jj)
 40   continue
      mcdndex(1,1,ii)=i
      mcdndex(1,2,ii)=kount
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccccc
ccccc
      subroutine rfcovinit(a,n1,n2)
cc
cc  Initializes the matrix a by filling it with zeroes.
cc
      double precision a(n1,n2)
cc
      do 100 i=1,n1
	do 90 j=1,n2
	  a(i,j)=0.D0
 90	continue
 100  continue
      return
      end
ccccc
ccccc
	subroutine rfcovsweep(a,nvar,k)
cc
	double precision a(nvar,nvar)
	double precision b
	double precision d
cc
	d=a(k,k)
	do 100 j=1,nvar
	  a(k,j)=a(k,j)/d
 100	continue
	do 1000 i=1,nvar
	  if(i.ne.k) then
	    b=a(i,k)
	    do 200 j=1,nvar
	      a(i,j)=a(i,j)-b*a(k,j)
 200	    continue
	    a(i,k)=-b/d
	  endif
 1000	continue
	a(k,k)=1/d
	return
	end
ccccc
ccccc
