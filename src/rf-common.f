c
c--   Routines common to
c--   fastLTS ( ./rfltsreg.f )  and
c--   fastMCD ( ./rffastmcd.f )
c
c
      subroutine rfrangen(n, nsel, index)
c
c     Randomly draws nsel cases out of n cases.
c     Here, index is the index set.
c
      implicit none
      integer n, nsel, index(nsel)
c     real uniran
      double precision unifrnd
      integer i,j, num
c
      do 100 i=1,nsel
c     OLD 10      num=int(uniran(seed)*n)+1
 10      num=int(unifrnd()*n)+1
C        if(num .gt. n) then
C           call intpr('** rfrangen(): num > n; num=', -1, num, 1)
C           num=n
C        endif
         if(i.gt.1) then
            do 50 j=1,i-1
               if(index(j).eq.num) goto 10
 50         continue
         endif
         index(i)=num
 100  continue
      return
      end
ccccc
ccccc
cOLD    function uniran(seed)
cOLD cc
cOLD cc  Draws a random number from the uniform distribution on [0,1].
cOLD cc
cOLD    real uniran
cOLD    integer seed
cOLD    integer quot
cOLD cc
cOLD    seed=seed*5761+999
cOLD    quot=seed/65536
cOLD    seed=seed-quot*65536
cOLD    uniran=float(seed)/65536.D0
cOLD    return
cOLD    end
ccccc
ccccc
      subroutine rfgenpn(n,nsel,index)
cc
cc    Constructs all subsets of nsel cases out of n cases.
cc
      implicit none
      integer n,nsel,index(nsel)
cc
      integer k,i

      k=nsel
      index(k)=index(k)+1
c while
 10   if(k.eq.1 .or. index(k).le.(n-(nsel-k))) goto 100
      k=k-1
      index(k)=index(k)+1
      do 50 i=k+1,nsel
         index(i)=index(i-1)+1
 50   continue
      goto 10
c end{while}
 100  return
      end
ccccc
ccccc
      subroutine rfshsort(a,n)
cc
cc  Sorts the array a of length n.
cc
      implicit none
      integer n
      double precision a(n)
c
      double precision t
      integer gap, i,j, nextj

      gap=n
 100  gap=gap/2
      if(gap.eq.0) goto 200
      do 180 i=1,n-gap
         j=i
 120     if(j.lt.1) goto 180
         nextj=j+gap
         if(a(j).gt.a(nextj)) then
            t=a(j)
            a(j)=a(nextj)
            a(nextj)=t
         else
            j=0
         endif
         j=j-gap
         goto 120
 180  continue
      goto 100
 200  return
      end
ccccc
ccccc
      subroutine rfishsort(a,kk)
cc
cc  Sorts the integer array a of length kk.
cc
      implicit none
      integer kk, a(kk)
c
      integer t, gap, i,j, nextj

      gap=kk
 100  gap=gap/2
      if(gap.eq.0) goto 200
      do 180 i=1,kk-gap
         j=i
 120     if(j.lt.1) goto 180
         nextj=j+gap
         if(a(j).gt.a(nextj)) then
            t=a(j)
            a(j)=a(nextj)
            a(nextj)=t
         else
            j=0
         endif
         j=j-gap
         goto 120
 180  continue
      goto 100
 200  return
      end
ccccc
ccccc
      function replow(k)
cc
cc    Find out which combinations of n and p are
cc    small enough in order to perform exaustive search
cc    Returns the maximal n for a given p, for which
cc    exhaustive search is to be done
cc
cc    k is the number of variables (p)
cc
      implicit none
      integer replow, k
c
      integer irep(6)
      data irep/500,50,22,17,15,14/
c
      if(k .le. 6) then
         replow = irep(k)
      else
         replow = 0
      endif
      return
      end
ccccc
ccccc
      integer function rfncomb(k,n)
cc
cc  Computes the number of combinations of k out of n.
cc  (To avoid integer overflow during the computation,
cc  ratios of reals are multiplied sequentially.)
cc  For comb > 1E+009 the resulting 'comb' may be too large
cc  to be put in the integer 'rfncomb', but the main program
cc  only calls this function for small enough n and k.
cc
      implicit none
      integer k,n
c
      double precision comb,fact
      integer j
c
      comb=dble(1.0)
      do 10 j=1,k
         fact=(dble(n-j+1.0))/(dble(k-j+1.0))
         comb=comb*fact
 10   continue
c     Should give error now instead of integer overflow!
c     Don't know how to get .Machine$integer.max in Fortran, portably
      if(comb .gt. 2147483647) then
         comb=2147483647.
         call
     +    dblepr('** too many combinations; using max.integer instead:',
     +           -1,comb,1)
      endif
      rfncomb=int(comb+0.5D0)
      return
      end
ccccc
ccccc
      subroutine rfcovcopy(a,b,n1,n2)
cc
cc  Copies matrix a to matrix b.
cc
      double precision a(n1,n2)
      double precision b(n1,n2)
c
      do 100 i=1,n1
         do 90 j=1,n2
            b(i,j)=a(i,j)
 90      continue
 100  continue
      return
      end
ccccc
ccccc
      function rffindq(aw,ncas,k,index)
cc
cc  Finds the k-th order statistic of the array aw of length ncas.
cc
cc MM{FIXME}: "rather" use R's C API   rPsort (double* X, int N, int K)

      implicit none
      integer ncas,k,index(ncas)
      double precision rffindq, aw(ncas)
c
      double precision ax,wa
      integer i,j,l,lr,jnc
cc
      do 10 j=1,ncas
         index(j)=j
 10   continue
      l=1
      lr=ncas
c    while(l < lr)
 20   if(l.ge.lr) goto 90
      ax=aw(k)
      jnc=l
      j=lr

 30   if(jnc.gt.j) goto 80
 40   if(aw(jnc).ge.ax) goto 50
      jnc=jnc+1
      goto 40
 50   if(aw(j).le.ax) goto 60
      j=j-1
      goto 50
 60   if(jnc .le. j) then
         i=index(jnc)
         index(jnc)=index(j)
         index(j)=i
         wa=aw(jnc)
         aw(jnc)=aw(j)
         aw(j)=wa
         jnc=jnc+1
         j=j-1
      endif
      goto 30
 80   if(j.lt.k) l=jnc
      if(k.lt.jnc) lr=j
      goto 20
 90   rffindq=aw(k)
      return
      end
ccccc
ccccc
      subroutine rfrdraw(a,n,ntot,mini,ngroup,kmini)
cc
cc  Draws ngroup nonoverlapping subdatasets out of a dataset of size n,
cc  such that the selected case numbers are uniformly distributed from 1 to n.
cc
      implicit none
      integer n, ntot, kmini, a(2,ntot), mini(kmini), ngroup
c
      double precision unifrnd
c
      integer jndex, nrand, k,m,i,j
cc
      jndex=0
      do 10 k=1,ngroup
         do 20 m=1,mini(k)
cOLD        nrand=int(uniran(seed)*(n-jndex))+1
            nrand=int(unifrnd()*(n-jndex))+1
C           if(nrand .gt. n-jndex) then
C              call intpr(
C      1         '** rfrdraw(): need to correct nrand > n-jndex; nrand=',
C      2                          -1, nrand, 1)
C              nrand=n-jndex
C           endif

            jndex=jndex+1
            if(jndex.eq.1) then
               a(1,jndex)=nrand
               a(2,jndex)=k
            else
               a(1,jndex)=nrand+jndex-1
               a(2,jndex)=k
               do 5,i=1,jndex-1
                  if(a(1,i).gt.nrand+i-1) then
                     do 6, j=jndex,i+1,-1
                        a(1,j)=a(1,j-1)
                        a(2,j)=a(2,j-1)
 6                   continue
                     a(1,i)=nrand+i-1
                     a(2,i)=k
                     goto 20
c                    ------- break
                  endif
 5             continue
            endif
 20      continue
 10   continue
      return
      end
ccccc
ccccc
      function rfodd(n)
cc
      logical rfodd
cc
      rfodd=.true.
      if(2*(n/2).eq.n) rfodd=.false.
      return
      end

ccccc
c unused        function rfnbreak(nhalf,n,nvar)
c unused cc
c unused cc  Computes the breakdown value - in percent! - of the MCD estimator
c unused cc
c unused         implicit none
c unused        integer rfnbreak, nhalf, n, nvar
c unused
c unused        if (nhalf.le.(n+nvar+1)/2) then
c unused          rfnbreak=(nhalf-nvar)*100/n
c unused        else
c unused          rfnbreak=(n-nhalf+1)*100/n
c unused        endif
c unused        return
c unused        end
ccccc

      subroutine rfmcduni(w,ncas,jqu,slutn,bstd,aw,aw2,factor,len)
cc
cc  rfmcduni : calculates the MCD in the univariate case.
cc           w contains the ordered observations
cc
c This version returns the index (jint) in 'len'
c which is used in rfltreg.f

      implicit double precision (a-h,o-z), integer(i-n)
      integer ncas, jqu, len
      double precision w(ncas), aw(ncas), aw2(ncas)
      double precision slutn(len)
cc
      sq=0.D0
      sqmin=0.D0
      ndup=1
      do 5 j=1,ncas-jqu+1
         slutn(j)=0.D0
 5    continue

      do 20 jint=1,ncas-jqu+1
         aw(jint)=0.D0
         do 10 j=1,jqu
            aw(jint)=aw(jint)+w(j+jint-1)
            if (jint.eq.1) sq=sq+w(j)*w(j)
 10      continue
         aw2(jint)=aw(jint)*aw(jint)/jqu
         if (jint.eq.1) then
            sq=sq-aw2(jint)
            sqmin=sq
            slutn(ndup)=aw(jint)
            len=jint
         else
            sq=sq - w(jint-1)*w(jint-1) + w(jint+jqu-1)*w(jint+jqu-1)
     *           - aw2(jint) + aw2(jint-1)
            if(sq.lt.sqmin) then
               ndup=1
               sqmin=sq
               slutn(ndup)=aw(jint)
               len=jint
            else
               if(sq.eq.sqmin) then
                  ndup=ndup+1
                  slutn(ndup)=aw(jint)
               endif
            endif
         endif
 20   continue
      slutn(1)=slutn(int((ndup+1)/2))/jqu
      bstd=factor*sqrt(sqmin/jqu)
      return
      end
ccccc
ccccc
