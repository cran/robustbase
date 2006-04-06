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

	subroutine rfltsreg(dat,n,nvar,nhalff,krep,inbest,objfct,
     *   intercept,intadjust,nvad,datt,iseed,
     *   weights,temp,index1,index2,aw2,aw,residu,y,nmahad,ndist,
     *   am,am2,slutn,
     *   jmiss,xmed,xmad,a,da,h,hvec,c,cstock,mstock,c1stock,
     *   m1stock,dath,sd,means,bmeans)
     
cc
	implicit integer(i-n), double precision(a-h,o-z)
cc
ccc	parameter (nvmax=115)
ccc	parameter (nmax=57000)
cc
	parameter (kmini=5)
	parameter (nmini=300)
	parameter (k1=2)
	parameter (k2=2)
	parameter (k3=100)
cc
ccc	parameter (nvmax1=nvmax+1)
ccc	parameter (nvmax2=nvmax*nvmax)
ccc	parameter (nvm11=nvmax*(nvmax+1))

	parameter (km10=10*kmini)
	parameter (nmaxi=nmini*kmini)
C--VT   parameter (maxmini=int((3*nmini-1)/2)+1)
	parameter (maxmini=450)
cc
	integer inbest(nhalff)
	double precision dat(n,nvad)
	double precision datt(n,nvad)

        double precision weights(n)
	integer temp(n)
	integer index1(n)
	integer index2(n)
        double precision aw2(n),aw(n)
	double precision residu(n)
	double precision y(n)
	double precision nmahad(n)
	double precision ndist(n)
	double precision am(n),am2(n),slutn(n)

	integer matz,seed,tottimes,step
	integer pnsel
        integer xreplow
	integer krep,n,nvar,nhalff,nvad
	double precision objfct
	logical all,part,fine,final,xrfodd,more1,more2
	integer intercept,intadjust
	integer xrfnbreak,xrfncomb
	
        integer flag(km10)
	integer mini(kmini)
	integer subdat(2,nmaxi)
	integer subndex(maxmini)
	double precision faclts(11)
	double precision mcdndex(10,2,kmini)


ccc     integer jmiss(nvmax1)
ccc	double precision xmed(nvmax1)
ccc	double precision xmad(nvmax1)
ccc	double precision a(nvmax1), da(nvmax1)
ccc	double precision h(nvmax,nvmax1),hvec(nvm11)
ccc	double precision c(nvmax,nvmax1)
ccc	double precision cstock(10,nvmax2)
ccc	double precision mstock(10,nvmax)
ccc	double precision c1stock(km10,nvmax2)
ccc	double precision m1stock(km10,nvmax)
ccc	double precision dath(nmaxi,nvmax1)
ccc	double precision sd(nvmax)
ccc	double precision means(nvmax)
ccc	double precision bmeans(nvmax)

        integer jmiss(nvad)
	double precision xmed(nvad)
	double precision xmad(nvad)
	double precision a(nvad), da(nvad)
	double precision h(nvar,nvad),hvec(nvar*nvad)
	double precision c(nvar,nvad)
	double precision cstock(10,nvar*nvar)
	double precision mstock(10,nvar)
	double precision c1stock(km10,nvar*nvar)
	double precision m1stock(km10,nvar)
	double precision dath(nmaxi,nvad)
	double precision sd(nvar)
	double precision means(nvar)
	double precision bmeans(nvar)

        data faclts/2.6477,2.5092,2.3826,2.2662,2.1587,
     *  2.0589,1.9660,1.879,1.7973,1.7203,1.6473/ 
cc
cc
CDDD	CALL INTPR('>>> Enter RFLTSREG ... nvar=',-1,nvar,1)

CCCC    10.10.2005 - substitute the parameters nmax and nvmax
        nmax = n
        nvmax = nvar        
        nvmax1 = nvmax+1
        nvm11 = nvmax*(nvmax+1)
        
        nrep = krep

	if(nvar.lt.5) then
		eps=1.0D-12
	else
		if(nvar.ge.5.and.nvar.le.8) then
			eps=1.0D-14
		else
			eps=1.0D-16
		endif
	endif
	prec=1.0D-6

cc	nhalff=int((n+nvar+1)/2)

	jmin=(n/2)+1
	jmax=max((3*n/4)+(nvar+1)/4,nhalff)
	nquant=min(nint(real(((nhalff*1.0/n)-0.5)*40))+1,11)
	factor=faclts(nquant)
	jbreak=xrfnbreak(nhalff,n,nvar)
	jdefaul=(n+nvar+1)/2
	percen = (1.D0*nhalff)/(1.D0*n)
	if(nvad.eq.1) goto 9000
cc
CDDD	CALL INTPR('>>> Enter RFLTSREG ... iseed=',-1,iseed,1)
        seed=iseed
        matz=1
        nsel=nvar
        nstop=0
        ngroup=1
        part=.false.
        fine=.false.
        final=.false.
        all=.true.
        do 21,i=1,nmaxi
          subdat(1,i)=1000000
          subdat(2,i)=1000000
 21     continue
cc
        mini(1)=0
        mini(2)=0
        mini(3)=0
        mini(4)=0
        mini(5)=0
        if(n.gt.(2*nmini-1)) then
          kstep=k1
          part=.true.
          ngroup=int(n/(nmini*1.D0))
          if(n.ge.(2*nmini) .and. n.le.(3*nmini-1)) then
            if(xrfodd(n)) then
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
cccc	CALL INTPR('>>> RFLTSREG ... minigr=',-1,iseed,1)
          call xrfrdraw(subdat,n,seed,minigr,mini,ngroup,kmini)
        else
          minigr=n
          nhalf=nhalff
          kstep=k1
          if(n.le.xreplow(nsel)) then
            nrep=xrfncomb(nsel,n)
          else
            nrep = krep
            all=.false.
          endif
        endif

        seed=iseed
cc

CDDD	CALL INTPR('>>> Start initialization ... nrep=',-1,nrep,1)

        do 31, j=1,nvmax
          do 33, k=1,10
            mstock(k,j)=1000000.D0
            do 35, kk=1,kmini
 35           m1stock((kk-1)*10+k,j)=1000000.D0
            do 37 i=1,nvmax
              do 39,kk=1,kmini
 39             c1stock((kk-1)*10+k,(j-1)*nvmax+i)=1000000.D0
 37         cstock(k,(j-1)*nvmax+i)=1000000.D0
 33       continue
          means(j)=0.D0
          bmeans(j)=0.D0
          sd(j)=0.D0
          do 46, k=1,nvmax1
              c(j,k)=0.D0
 46       h(j,k)=0.D0
 31     continue

        do 41, j=1,nmax
          nmahad(j)=0.D0
          ndist(j)=0.D0
          index1(j)=1000000
          index2(j)=1000000
          temp(j)=1000000
          weights(j)=0.D0
          aw(j)=0.D0
          aw2(j)=0.D0
          residu(j)=0.D0
          y(j)=0.D0
          am(j)=0.D0
          am2(j)=0.D0
          slutn(j)=0.D0
 41     continue

        do 43,j=1,km10
 43       flag(j)=1
        do 45, j=1,nvmax1
          jmiss(j)=0
          xmed(j)=0.D0
          xmad(j)=0.D0
          a(j)=0.D0
          da(j)=0.D0
          do 48,k=1,nmaxi
 48         dath(k,j)=0.D0
 45     continue

        do 44, j=1,maxmini
          subndex(j)=0.D0
 44     continue
        do 42,j=1,nhalff
          inbest(j)=0
 42     continue

        do 47,j=1,nvm11
 47       hvec(j)=0.D0

CDDD	CALL INTPR('>>> Initialization ready',-1,0,0)

 9000   continue
        if(nvad.eq.1) then
          do 23, jj=1,n
 23         ndist(jj)=dat(jj,1)
          call xrfshsort(ndist,n)
          call xrfmcduni(ndist,n,nhalff,slutn,bstd,am,am2,factor,
     *      n-nhalff+1)
          goto 9999
        endif
cc
        if(.not.fine.and..not.final) then
          call rfstatis(dat,xmed,xmad,aw,aw2,intercept,nvad,
     *      nvmax1,nmax,n,
     *      nstop,prec,weights,y,nvar,index2)
          if(nstop.eq.1) goto 9999
        endif

cc
        jreg=1
        call rflsreg(nvmax1,nmax,nvmax,nvar,n,a,dat,y,weights,da,
     *    h,fckw,hvec,nvm11,jmiss,nvad,n)
CC        nfac=nvad-1
        nfac=nvar-1
        call rfrtran(nvar,intercept,nfac,nvad,nvmax1,xmed,
     *    xmad,a,nvad,fckw)
        call rftrc(h,da,nvmax,nvmax1,nvar,intercept,nfac,nvad,
     *    xmed,xmad)
        jerd=0

        tottimes=0
 5555   object=10.D25
        if(.not. part .or. final) then
          nn=n
        endif
        if(part .and. fine .and. .not. final) nn=minigr
 
        if(fine.or.(.not.part.and.final)) then
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
     *            .le.300000) then
                  kstep=9
                else
                  if (n*nvar .gt.300000 .and. n*nvar
     *              .le.400000) then
                    kstep=8
                  else
                    if (n*nvar .gt.400000 .and. n*nvar
     *                .le.500000) then
                      kstep=7
                    else
                      if (n*nvar .gt.500000 .and. n*nvar
     *                  .le.600000) then
                        kstep=6
                      else
                        if (n*nvar .gt.600000 .and. n*nvar
     *                    .le.700000) then
                          kstep=5
                        else
                          if (n*nvar .gt.700000 .and. n*nvar
     *                      .le.800000) then
                            kstep=4
                          else
                            if (n*nvar .gt.800000 .and. n*nvar
     *                        .le.900000) then
                              kstep=3
                            else
                              if (n*nvar .gt.900000 .and. n*nvar
     *                          .le.1000000) then
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
 81     continue
        index1(nsel)=nsel-1
cc
        if(.not. final) then
          do 83 i=1,10
            do 85 j=1,ngroup
              mcdndex(i,1,j)=10.D25
 85         mcdndex(i,2,j)=10.D25
 83       continue
        endif
        if (fine .and. .not. final) then
          do 91, j=1,minigr
            do 93, k=1,nvad
 93           dath(j,k)=dat(subdat(1,j),k)
 91       continue
        endif
        kount=0

CDDD	CALL INTPR('>>> MAIN LOOP BY GROUPS: NGROUP= ',-1,ngroup,1)

        do 1111 ii=1,ngroup
CDDD	     CALL INTPR('>>> LOOPING BY GROUPS...II: ',-1,ii,1)
          
          if(.not.fine) kount=0
          if(part .and. .not. fine) nn=mini(ii)
          do 101 i=1,nn
              index2(i)=i
 101      continue
          if(part .and. .not. fine) then
              jndex=0
              do 103 j=1,minigr
                  if(subdat(2,j).eq.ii) then
                      jndex=jndex+1
                      subndex(jndex)=subdat(1,j)
                  endif
 103          continue
              do 105 j=1,mini(ii)
                  do 107 k=1,nvad
 107                  dath(j,k)=dat(subndex(j),k)
 105          continue
          endif
cc
CDDD      CALL INTPR('>>> MAIN LOOP: NREP=',-1,nrep,1)
          do 1000 i=1,nrep
CDDD               CALL INTPR('>>> LOOPING...I: ',-1,i,1)

          pnsel=nsel
          tottimes=tottimes+1
          fckwi=0.D0
          fckw1=0.D0
          step=0
 132      if((part.and..not.fine).or.(.not.part.and..not.final)) then
              if(part) then
                  call xrfrangen(mini(ii),nsel,index1,seed)
              else
                  if(all) then
                      call xrfgenpn(n,nsel,index1)
                  else
                      call xrfrangen(n,nsel,index1,seed)
                  endif
              endif
          endif

 9550     continue
          if(.not.fine.and.part) then
              do 121 j=1,pnsel
                  do 123 m=1,nvad
 123                  c(j,m)=dath(index1(j),m)
 121          continue
          endif
          if(.not.part.and..not.final) then
              do 122 j=1,pnsel
                  do 124 m=1,nvad
 124                  c(j,m)=dat(index1(j),m)
 122          continue
          endif
          if((.not.part.and..not.final).or.(.not.fine.and.part)) then
              if(nvar.gt.1) then
                  call rfequat(c,nvmax,nvmax1,hvec,nvm11,nvar,1,nerr)
                  if(nerr.ge.0) goto 126
                  jerd=jerd+1
                  if(.not.all.and.i.gt.2) goto 132
                  goto 1000
              else
                  if(c(1,1).ne.0.D0) c(1,1)=c(1,2)/c(1,1)
              endif

 126          continue
 
              do 136 jnc=1,nvar
 136              a(jnc)=c(jnc,1)
          endif
          if (final) then
              if(mstock(i,1).ne.1000000.D0) then
                  do 125 jj=1,nvar
                      a(jj)=mstock(i,jj)
 125              continue
              else
                  goto 1111
              endif
          endif
          if (fine.and..not.final) then
              if(m1stock((ii-1)*10+i,1).ne.1000000.D0) then
                  do 131 jj=1,nvar
                      a(jj)=m1stock((ii-1)*10+i,jj)
 131              continue
              else
                  goto 1111
              endif
          endif
            
 151      do 152 jnc=1,nn
              residu(jnc)=0.D0
              do 153 j=1,nvar
                  if(part.and..not.final) then
                      residu(jnc)=residu(jnc)+a(j)*dath(jnc,j)
                  else
                      residu(jnc)=residu(jnc)+a(j)*dat(jnc,j)
                  endif
 153          continue
              if(part.and..not.final) then
                  residu(jnc)=dath(jnc,nvad)-residu(jnc)
              else
                  residu(jnc)=dat(jnc,nvad)-residu(jnc)
              endif
              aw(jnc)=residu(jnc)
 152      continue
          
          more1=.false.
          more2=.false.
          nmore=200 
          nmore2=nmore/2

          if(intadjust.eq.1) then

CDDD		CALL INTPR('>>> INTERCEPT ADJUSTMENT 1',-1,i,1)
          if(intercept.eq.1.and.((.not.fine.and.part).or.
     *      .not.part.or.((nn-nhalf).le.nmore))) then
              call xrfshsort(aw,nn)
              call xrfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *            factor,nn-nhalf+1)
              a(nvar)=a(nvar)+slutn(1)
              do 154 jnc=1,nn
 154              residu(jnc)=residu(jnc)-slutn(1)
          else
              if(intercept.eq.1) then
                  call xrfshsort(aw,nn)
                  do 184 jj=1,nn
 184                  am2(jj)=abs(aw(jj))
                  dist2=xrffindq(am2,nn,nhalf,index1)
                  do 174, jj=1,nhalf
 174                  aw2(jj)=aw(index1(jj))
                  dist2=xrffindq(aw2,nhalf,1,index2) 
                  jnc=index1(index2(1))
                  if(jnc+nmore-nmore2+nhalf-1.gt.nn.or.jnc-nmore2.lt.1) 
     *                then
                      call xrfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *                    factor,nn-nhalf+1)
                      a(nvar)=a(nvar)+slutn(1)
                      do 169 jnc=1,nn
 169                      residu(jnc)=residu(jnc)-slutn(1)
                  else
 555                  do 178 jj=0,nhalf-1+nmore
 178                      aw2(jj+1)=aw(jnc-nmore2+jj) 
                      nlen=nmore+1
                      call xrfmcduni(aw2,nhalf+nmore,nhalf,slutn,
     *          bstd,am,am2,factor,nlen)
                      if(nlen.eq.1.and..not.more1) then
                          if(.not.more2) then
                              nmore=nmore2
                              nmore2=nmore2+nmore2
                              more1=.true.
                              if(jnc-nmore2.ge.1) goto 555
                          endif
                      else
                          if(nlen.eq.(nmore+1).and..not.more2) then
                              if(.not.more1) then
                                  nmore=nmore2
                                  nmore2=-nmore2
                                  more2=.true.
                                  if(jnc+nmore-nmore2+nhalf-1.le.nn) 
     *                                goto 555
                              endif
                          else
                              if(nlen.eq.1.and.more1) then
                                  if(.not.more2) then
                                      nmore2=nmore2+100
                                      if(jnc-nmore2.ge.1) goto 555
                                  endif
                              else
                                  if(nlen.eq.(nmore+1).and.more2) then
                                      if(.not.more1) then
                                          nmore2=nmore2+100
                            if(jnc+nmore-nmore2+nhalf-1.le.nn) goto 555
                                      endif
                                  endif
                              endif
                          endif
                      endif
                      a(nvar)=a(nvar)+slutn(1)
                      do 170 jnc=1,nn
 170                      residu(jnc)=residu(jnc)-slutn(1)
                  endif
              endif
          endif
          endif
          
          do 156 jnc=1,nn
 156        residu(jnc)=abs(residu(jnc))
          dist2=xrffindq(residu,nn,nhalf,index2)
 9555     do 400 step=1,kstep
            tottimes=tottimes+1
            do 155, j=1,nhalf
 155          temp(j)=index2(j)
            call xrfishsort(temp,nhalf)
            do 157 j=1,nhalf
              if(.not.part.or.final) then
                do 158 mm=1,nvad
  158             datt(j,mm)=dat(temp(j),mm)
              else
                do 159 mm=1,nvad
  159             datt(j,mm)=dath(temp(j),mm)
              endif
 157        continue
            call rflsreg(nvmax1,nmax,nvmax,nvar,n,a,datt,y,weights,
     *        da,h,fckw,hvec,nvm11,jmiss,nvad,nn)
 171        do 172 jnc=1,nn
              residu(jnc)=0.D0
            do 173 j=1,nvar
              if(part.and..not.final) then
                residu(jnc)=residu(jnc)+a(j)*dath(jnc,j)
              else
                residu(jnc)=residu(jnc)+a(j)*dat(jnc,j)
              endif
 173        continue
            if(part.and..not.final) then
              residu(jnc)=dath(jnc,nvad)-residu(jnc)
            else
              residu(jnc)=dat(jnc,nvad)-residu(jnc)
            endif
            aw(jnc)=residu(jnc)
 172      continue

          more1=.false.
          more2=.false.
          nmore=200 
          nmore2=nmore/2
          
          if(intadjust.eq.1) then

CDDD		CALL INTPR('>>> INTERCEPT ADJUSTMENT 2',-1,step,1)
          if(intercept.eq.1.and.((.not.fine.and.part).or.
     *      .not.part.or.((nn-nhalf).le.nmore))) then
              call xrfshsort(aw,nn)
              call xrfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *            factor,nn-nhalf+1)
              a(nvar)=a(nvar)+slutn(1)
              do 179 jnc=1,nn
 179              residu(jnc)=residu(jnc)-slutn(1)
          else
              if(intercept.eq.1) then
                  call xrfshsort(aw,nn)
                  do 185 jj=1,nn
 185                  am2(jj)=abs(aw(jj))
                  dist2=xrffindq(am2,nn,nhalf,index1)
                  do 180, jj=1,nhalf
 180                  aw2(jj)=aw(index1(jj))
                  dist2=xrffindq(aw2,nhalf,1,index2) 
                  jnc=index1(index2(1))
                  if(jnc+nmore-nmore2+nhalf-1.gt.nn.or.jnc-nmore2.lt.1)
     *                then
                      call xrfmcduni(aw,nn,nhalf,slutn,bstd,am,am2,
     *                    factor,nn-nhalf+1)
                      a(nvar)=a(nvar)+slutn(1)
                      do 168 jnc=1,nn
 168                      residu(jnc)=residu(jnc)-slutn(1)
                  else
 666                  do 181 jj=0,nhalf-1+nmore
 181                      aw2(jj+1)=aw(jnc-nmore2+jj) 
                      nlen=nmore+1
                      call xrfmcduni(aw2,nhalf+nmore,nhalf,slutn,bstd,
     *                     am,am2,factor,nlen)
                      if(nlen.eq.1.and..not.more1) then
                          if(.not.more2) then
                              nmore=nmore2
                              nmore2=nmore2+nmore2
                              more1=.true.
                              if(jnc-nmore2.ge.1) goto 666 
                          endif
                      else
                          if(nlen.eq.(nmore+1).and..not.more2) then
                              if(.not.more1) then
                                  nmore=nmore2
                                  nmore2=-nmore2
                                  more2=.true.
                                  if(jnc+nmore-nmore2+nhalf-1.le.nn) 
     *                                goto 666 
                              endif
                          else
                              if(nlen.eq.1.and.more1) then
                                  if(.not.more2) then
                                      nmore2=nmore2+100
                                      if(jnc-nmore2.ge.1) goto 666 
                                  endif
                              else
                                  if(nlen.eq.(nmore+1).and.more2) then
                                      if(.not.more1) then
                                          nmore2=nmore2+100
                            if(jnc+nmore-nmore2+nhalf-1.le.nn) goto 666 
                                      endif
                                  endif
                              endif
                          endif
                      endif
                      a(nvar)=a(nvar)+slutn(1)
                      do 182 jnc=1,nn
 182                      residu(jnc)=residu(jnc)-slutn(1)
                  endif
              endif
          endif
          endif
          
          do 177 jnc=1,nn
 177        residu(jnc)=abs(residu(jnc))
            dist2=xrffindq(residu,nn,nhalf,index2)
            fckw=0.D0
            do 176 jnc=1,nhalf
 176          fckw=fckw+residu(jnc)**2

            if(step.ge.2 .and. fckw.eq.fckw1) then
              goto 5000
            endif
            fckw1=fckwi
            fckwi=fckw
            if(((i.eq.1.and.step.eq.1.and..not.fine)
     *        .or.fckw.lt.object).and.(final)) then
              object=fckw
              objfct=fckw
              do 175 jjj=1,nhalf
                inbest(jjj)=index2(jjj)
 175          continue
              call xrfcovcopy(a,bmeans,nvar,1)
            endif
 400      continue
cc
 5000     if(.not. final) then
            if(part .and. .not. fine) then
              iii=ii
            else
              iii=1
cc            At the end of the algorithm, only the ten
cc            best solutions need to be stored.
            endif
 
          if( flag((iii-1)*10+1).eq.1) then
            lll=1
          else
            lll=2
          endif
          do 201, j=lll,10
            if (fckw .le. mcdndex(j,2,iii)) then
              if(fckw.ne.mcdndex(j,2,iii)) then
                if(.not.fine.and.part) goto 203
                goto 205
              else
                do 207 kkk=j,10
                  if(fckw.eq.mcdndex(kkk,2,iii)) then
                    do 209, jjj=1,nvar
                      if(part.and..not.fine) then
                        if(a(jjj).ne.m1stock((iii-1)*10+
     *                    kkk,jjj)) then
                            goto 203
                        endif
                      else
                        if(a(jjj).ne.mstock(kkk,jjj)) then
                          goto 205
                        endif
                      endif
  209               continue
                  endif
  207           continue
              endif
              goto 1000
  203         do 221,k=10,j+1,-1
                do 225 kk=1,nvar
  225             m1stock((iii-1)*10+k,kk)=
     *            m1stock((iii-1)*10+k-1,kk)
                mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
  221         continue
              do 227 kk=1,nvar
                m1stock((iii-1)*10+j,kk)=a(kk)
  227         continue
              mcdndex(j,1,iii)=i
              mcdndex(j,2,iii)=fckw
              goto 1000
  205         do 231,k=10,j+1,-1
                do 235 kk=1,nvar
  235             mstock(k,kk)=
     *            mstock(k-1,kk)
                mcdndex(k,1,iii)=mcdndex(k-1,1,iii)
                mcdndex(k,2,iii)=mcdndex(k-1,2,iii)
  231         continue
              do 237 kk=1,nvar
                mstock(j,kk)=a(kk)
  237         continue
              mcdndex(j,1,iii)=i
              mcdndex(j,2,iii)=fckw
              goto 1000
            endif
  201     continue
        endif
 1000 continue
 1111 continue
cc
      if(part .and. .not. fine) then
        fine=.true.
        goto 5555
      endif
 5550 if((part.and.fine.and..not.final).or.
     *  (.not.part.and..not.final)) then
        final=.true.
        goto 5555
      endif
cc
      goto 9999
 9999 return
      end
ccccc
ccccc

      subroutine rfstatis(x,xmed,xmad,aw,aw2,intercept,nvad,nvmax1,
     *  nmax,n,nstop,prec,weights,y,nvar,index2)
cc
cc
      double precision xmed(nvmax1), x(n,nvad), xmad(nvmax1)
      double precision aw(nmax), aw2(nmax)
      double precision weights(nmax)
      double precision y(nmax)
      double precision prec
      double precision rfamdan
      integer index2(nmax)
cc
      if (intercept.eq.1) goto 60
cc regression without intercept
      do 50 j=1,nvad
        xmed(j)=0.0
        do 10 jnc=1,n
 10       aw2(jnc)=abs(x(jnc,j))
        xmad(j)=rfamdan(aw,nmax,aw2,n,index2)*1.4826
        if(abs(xmad(j)).gt.prec) goto 30
        xmad(j)=0.0
        do 20 jnc=1,n
 20       xmad(j)=xmad(j)+aw2(jnc)
        xmad(j)=(xmad(j)/n)*1.2533
        if(abs(xmad(j)).gt.prec) goto 30
        nstop=1
        return
 30     do 40 jnc=1,n
          x(jnc,j)=x(jnc,j)/xmad(j)
 40     continue
 50   continue
      goto 150
cc regression with intercept
 60   xmed(nvar)=0.D0
      xmad(nvar)=1.D0
      do 120 j=1,nvad
      if(j.eq.nvar) goto 120
      do 70 jnc=1,n
 70     aw2(jnc)=x(jnc,j)
      xmed(j)=rfamdan(aw,nmax,aw2,n,index2)
      do 80 jnc=1,n
 80     aw2(jnc)=abs(aw2(jnc)-xmed(j))
      xmad(j)=rfamdan(aw,nmax,aw2,n,index2)*1.4826
      if(abs(xmad(j)).gt.prec) goto 100
      xmad(j)=0.0
      do 90 jnc=1,n
 90     xmad(j)=xmad(j)+aw2(jnc)
      xmad(j)=(xmad(j)/n)*1.2533
      if(dabs(xmad(j)).gt.prec) goto 100
      nstop=1
      return
 100  do 110 jnc=1,n
        x(jnc,j)=(x(jnc,j)-xmed(j))/xmad(j)
 110  continue
 120  continue
 150  continue
      do 270, jnc=1,n
        weights(jnc)=1.0
        y(jnc)=x(jnc,nvad)
 270  continue
      return
      end
cc
      function rfamdan(aw,nmax,aa,n,index2)
cc
      double precision aa(n),aw(nmax)
      integer index2(nmax)
      double precision xrffindq
      double precision rfamdan
cc
      jndl=int(n/2.0)
      if(mod(n,2).eq.0) then
        rfamdan=(xrffindq(aa,n,jndl,index2)+
     *    xrffindq(aa,n,jndl+1,index2))/2.0
      else
        rfamdan=xrffindq(aa,n,jndl+1,index2)
      endif
      return
      end
cc
      subroutine rflsreg(nvmax1,nmax,nvmax,k,n,f,x,y,w,da,h,fckw,
     *  hvec,nvm11,jmiss,nvad,nnn)
cc
      double precision x(n,nvad),f(k),y(n),w(n),da(k)
      double precision hvec(nvm11),h(nvmax,nvmax1)
      double precision fckw,dfckw,dfact
      double precision dwjnc,dyj,dfka
      integer jmiss(nvmax1)
cc
      kplus=k+1
      do 10 jnc=1,k
        do 20 j=1,kplus
          h(jnc,j)=0.D0
 20     continue
 10   continue
      anul=0.0
      do 30 jnc=1,nnn
        call rffcn(k,f,x,jnc,nmax,nvmax1,n,nvad)
        dwjnc=dble(w(jnc))
        anul=anul+w(jnc)
        dyj=dble(x(jnc,kplus))
        do 40 ka=1,k
          dfka=dble(f(ka))
          h(ka,k+1)=h(ka,k+1)+dwjnc*dfka*dyj
          do 50 l=1,ka
            h(ka,l)=h(ka,l)+dwjnc*dfka*dble(f(l))
 50       continue
 40     continue
 30   continue
      do 60 j=1,k
        do 70 jnc=1,j
          h(jnc,j)=h(j,jnc)
 70     continue
 60   continue
      call rfmatnv(h,nvmax,nvmax1,hvec,nvm11,k,1,jmiss)
      mm=k+1
      fckw=rfqlsrg(k,n,nvmax1,nmax,nvmax,f,x,y,w,h,mm,nvad,nnn)
      do 80 jnc=1,k
        f(jnc)=h(jnc,k+1)
 80   continue
      dfckw=dble(fckw)
      ank=anul-k
      dfact=dble(ank)
      dfact=dfckw/dfact
      do 90 jnc=1,k
        do 100 j=1,k
          h(jnc,j)=h(jnc,j)*dfact
 100    continue
 90   continue
      do 110 jnc=1,k
        hda=h(jnc,jnc)
        da(jnc)=sqrt(hda)
 110  continue
      return
      end
ccccc
ccccc
      subroutine rffcn(k,f,x,jnc,nmax,nvmax1,n,nvad)
cc
      double precision f(k),x(n,nvad)
cc
      do 10,j=1,k
        f(j)=x(jnc,j)
 10   continue
      return
      end
ccccc
ccccc
      subroutine rfmatnv(an,nvmax,nvmax1,hvec,nvm11,na,nb,
     *  jmiss)
cc
      double precision deter,turn,swap
      double precision hvec(nvm11),an(nvmax,nvmax1)
      integer jmiss(nvmax1)
      integer ldel
cc
      deter=1.0D0
      n=na
      npnb=n+nb
      jnk=0
      do 10 j=1,npnb
        jnk=(j-1)*nvmax
        do 10 nc=1,nvmax
          jnk=jnk+1
          hvec(jnk)=an(nc,j)
 10   continue
      ldel=0
      jdm=nvmax
      nma=n-1
      jdelc=1-jdm
      do 130 jhfd=1,n
        turn=0.0D0
        jdelc=jdelc+jdm
        jdla=jdelc+jhfd-1
        jdlb=jdelc+nma
        do 40 jncb=jdla,jdlb
          if(dabs(hvec(jncb))-dabs(turn)) 40,40,30
 30           turn=hvec(jncb)
              ldel=jncb
 40     continue
        if (turn) 50,180,50
 50     jpaal=ldel-jdelc+1
        jmiss(jhfd)=jpaal
        if(jpaal-jhfd) 80,80,60
 60     deter=-deter
        jpaal=jpaal-jdm
        jncd=jhfd-jdm
        do 70 jnc=1,npnb
          jpaal=jpaal+jdm
          jncd=jncd+jdm
          swap=hvec(jncd)
          hvec(jncd)=hvec(jpaal)
 70     hvec(jpaal)=swap
 80     deter=deter*turn
        turn=1.0D0/turn
        jncd=jdelc+nma
        do 90 jnc=jdelc,jncd
 90       hvec(jnc)=-hvec(jnc)*turn
        hvec(jdla)=turn
        jncb=jhfd-jdm
        jpaal=1-jdm
        do 120 jnc=1,npnb
          jpaal=jpaal+jdm
          jncb=jncb+jdm
          if(jnc-jhfd) 100,120,100
 100      jcl=jpaal+nma
          swap=hvec(jncb)
          jncd=jdelc-1
          do 110 jncc=jpaal,jcl
            jncd=jncd+1
 110        hvec(jncc)=hvec(jncc)+swap*hvec(jncd)
          hvec(jncb)=swap*turn
 120    continue
 130  continue
      do 160 jncb=1,n
        jhfd=n+1-jncb
        ldel=jmiss(jhfd)
        if(ldel-jhfd) 140,160,140
 140    jpaal=(ldel-1)*jdm+1
        jcl=jpaal+nma
        jdelc=(jhfd-1)*jdm+1-jpaal
        do 150 jncc=jpaal,jcl
          jncd=jncc+jdelc
          swap=hvec(jncc)
          hvec(jncc)=hvec(jncd)
 150      hvec(jncd)=swap
 160  continue
      goto 180
 180  jnk=0
      do 190 j=1,npnb
        do 190 nc=1,nvmax
          jnk=jnk+1
          an(nc,j)=hvec(jnk)
 190  continue
      return
      end
ccccc
ccccc
      function rfqlsrg(k,n,nvmax1,nmax,nvmax,f,x,y,w,h,mm,nvad,nnn)
cc
      double precision f(k), x(n,nvad), y(n), w(n)
      double precision q,hsum,h(nvmax,nvmax1)
cc
      q=0.D0
      do 30 jnc=1,nnn
        call rffcn(k,f,x,jnc,nmax,nvmax1,n,nvad)
        hsum=0.D0
        do 20 jncb=1,k
 20       hsum=h(jncb,mm)*f(jncb)+hsum
 30   q=(hsum-x(jnc,mm))*(hsum-x(jnc,mm))*w(jnc)+q
      rfqlsrg=q
      return
      end
ccccc
ccccc
      subroutine rfrtran(nvar,jcst,nfac,nvad,nvmax1,xmed,xmad,
     *  aa,jal,fckw)
cc
      double precision aa(jal), xmed(nvmax1), xmad(nvmax1)
      double precision fckw
cc
      if(nvar.le.1) then
        aa(1)=aa(1)*xmad(nvad)/xmad(1)
        goto 30
      endif
      do 10 j=1,nfac
 10     aa(j)=aa(j)*xmad(nvad)/xmad(j)
      if(jcst.eq.0) then
        aa(nvar)=aa(nvar)*xmad(nvad)/xmad(nvar)
      else
        aa(nvar)=aa(nvar)*xmad(nvad)
        do 20 j=1,nfac
 20       aa(nvar)=aa(nvar)-aa(j)*xmed(j)
        aa(nvar)=aa(nvar)+xmed(nvad)
      endif
 30   fckw=fckw*(xmad(nvad)*xmad(nvad))
      return
      end
ccccc
ccccc
      subroutine rftrc(h,da,nvmax,nvmax1,nvar,jcst,nfac,nvad,
     *  xmed,xmad)
cc
      double precision da(nvmax)
      double precision xmed(nvmax1),xmad(nvmax1)
      double precision h(nvmax,nvmax1),xmp2,hnn
cc
      xmp2=dble(xmad(nvad))*dble(xmad(nvad))
      if(jcst.eq.0) then
        do 10 j=1,nvar
          do 20 k=1,j
            h(j,k)=h(j,k)*(xmp2/(dble(xmad(j))*dble(xmad(k))))
 20       continue
          da(j)=dsqrt(h(j,j))
 10     continue
      else
        do 25 j=1,nvar
          h(j,nvad)=h(j,j)
 25     continue
        do 30, j=1,nvar
          do 40 k=1,j
            h(j,k)=h(j,k)*xmp2/(dble(xmad(j))*dble(xmad(k)))
 40       continue
          da(j)=dsqrt(h(j,j))
 30     continue
        do 50 k=1,nfac
          h(nvar,k)=h(k,nvar)*xmp2/dble(xmad(k))
          do 60 k2=1,nvar
            if(k.eq.k2) then
              h(nvar,k)=h(nvar,k)-dble(xmed(k))*xmp2/
     *          (dble(xmad(k2))*dble(xmad(k)))*h(k2,nvad)
              goto 60
            endif
            if(k.lt.k2) then
              h(nvar,k)=h(nvar,k)-(dble(xmed(k2))*xmp2)/
     *          (dble(xmad(k2))*dble(xmad(k)))*h(k,k2)
            else
              h(nvar,k)=h(nvar,k)-dble(xmed(k2))*xmp2/
     *          (dble(xmad(k2))*dble(xmad(k)))*h(k2,k)
            endif
 60       continue
 50     continue
        h(nvar,nvar)=h(nvar,nvad)*xmp2
        do 70 k=1,nvar
          h(nvar,nvar)=h(nvar,nvar)+
     *      (dble(xmed(k))*dble(xmed(k)))*xmp2/
     *      (dble(xmad(k))*dble(xmad(k)))*h(k,nvad)
 70      continue
         do 80 k=1,nvar
           if(k.ne.nvar) then
             h(nvar,nvar)=h(nvar,nvar)-2.0D0*xmp2*dble(xmed(k))/
     *         (dble(xmad(k)))*h(k,nvar)
           else
             h(nvar,nvar)=h(nvar,nvar)-2.0D0*xmp2*dble(xmed(k))/
     *         (dble(xmad(k)))*h(nvar,nvad)
           endif
 80      continue
         do 90 j=1,nfac
           ju=j+1
           do 90 k=ju,nvar
             hnn=2.0D0*dble(xmed(j))*dble(xmed(k))*xmp2
             h(nvar,nvar)=h(nvar,nvar)+hnn/
     *         (dble(xmad(j))*dble(xmad(k)))*h(j,k)
 90      continue
         da(nvar)=dsqrt(h(nvar,nvar))
       endif
       return
       end
ccccc
ccccc
      subroutine rfequat(am,nvmax,nvmax1,hvec,nvm11,na,nb,nerr)
      double precision am(nvmax,nvmax1)
      double precision hvec(nvm11),turn,swap,deter
      integer ldel

      ldel=0
      jdm=nvmax
      deter=1.0D0
      n=na
      jmat=n+nb
      jnk=0
      do 10 j=1,jmat
         jnk=(j-1)*nvmax
         do 10 nc=1,nvmax
            jnk=jnk+1
            hvec(jnk)=am(nc,j)
 10   continue

      nznde=n-1
      lclpl=-jdm
      do 120 jhfd=1,n
         turn=0.D0
         lclpl=lclpl+jdm+1
         jdel=lclpl+n-jhfd
         do 40 jncb=lclpl,jdel
         if(dabs(hvec(jncb))-dabs(turn)) 40,40,30
 30         turn=hvec(jncb)
            ldel=jncb
 40      continue
         if(dabs(turn).le.1D-8) goto 170
         if(ldel-lclpl) 60,80,60
 60         deter=-deter
            ldel=ldel-jdm
            jncb=lclpl-jdm
            do 70 jncc=jhfd,jmat
               ldel=ldel+jdm
               jncb=jncb+jdm
               swap=hvec(jncb)
               hvec(jncb)=hvec(ldel)
 70         hvec(ldel)=swap
 80         deter=deter*turn
         if(jhfd.eq.n) goto 120
         turn=1./turn
         jncb=lclpl+1
         do 90 jncc=jncb,jdel
 90         hvec(jncc)=hvec(jncc)*turn
         jncd=lclpl
         jrow=jhfd+1
         do 110 jncb=jrow,n
            jncd=jncd+1
            jnce=lclpl
            jncf=jncd
            do 100 jncc=jrow,jmat
               jnce=jnce+jdm
               jncf=jncf+jdm
 100        hvec(jncf)=hvec(jncf)-hvec(jnce)*hvec(jncd)
 110     continue
 120  continue
      
      nerr=0
      neqa=n+1
      jbegx=nznde*jdm+1
      do 150 jnc=neqa,jmat
      jbegx=jbegx+jdm
      jendx=jbegx+n
      jbegc=n*jdm+1
      jendc=jbegc+nznde
      do 140 jncb=1,nznde
      jendx=jendx-1
      jbegc=jbegc-jdm
      jendc=jendc-jdm-1
      hvec(jendx)=hvec(jendx)/hvec(jendc+1)
      swap=hvec(jendx)
      jncd=jbegx-1
      do 130 jncc=jbegc,jendc
      jncd=jncd+1
 130  hvec(jncd)=hvec(jncd)-hvec(jncc)*swap
 140  continue
      hvec(jbegx)=hvec(jbegx)/hvec(1)
 150  continue
      jnc=-jdm
      jbegx=nznde*jdm+1
      jendx=jbegx+nznde
      do 160 jncb=neqa,jmat
      jbegx=jbegx+jdm
      jendx=jendx+jdm
      jnc=jnc+jdm
      jncd=jnc
      do 160 jncc=jbegx,jendx
      jncd=jncd+1
 160  hvec(jncd)=hvec(jncc)
      goto 180
 170  nerr=-1
 180  jnk=0
      do 190 j=1,jmat
      do 190 nc=1,nvmax
      jnk=jnk+1
      am(nc,j)=hvec(jnk)
 190  continue
      return
      end
ccccc
C--VT-- The following functions were added
ccccc
ccccc
      function xreplow(k)
cc
cc    Find out which combinations of n and p are
cc    small enough in order to perform exaustive search
cc    Returns the maximal n for a given p, for which
cc    exhaustive search is to be done
cc    
cc    k is the number of variables (p)
cc
      integer xreplow, k
      integer irep(6)
      data irep/500,50,22,17,15,14/

      iret=0
      if(k.le.6) iret = irep(k)

      xreplow = iret
      return
      end 
ccccc
ccccc
      function xrfncomb(k,n)
cc
cc  Computes the number of combinations of k out of n.
cc  (To avoid integer overflow during the computation,
cc  ratios of reals are multiplied sequentially.)
cc  For comb > 1E+009 the resulting 'comb' may be too large
cc  to be put in the integer 'ncomb', but the main program
cc  only calls this function for small enough n and k.
cc
      integer xrfncomb,k,n
      double precision comb,fact
cc
      comb=dble(1.0)
      do 10 j=1,k
      fact=(dble(n-j+1.0))/(dble(k-j+1.0))
      comb=comb*fact
 10   continue
      xrfncomb=int(comb+0.5D0)

      return
      end
ccccc
ccccc
	subroutine xrfrangen(n,nsel,index,seed)
cc
cc    Randomly draws nsel cases out of n cases.
cc    Here, index is the index set.
cc
	integer seed
	integer index(nsel)
	real uniran
cc
	do 100 i=1,nsel
 10       num=int(uniran(seed)*n)+1
	  if(i.gt.1) then
	    do 50 j=1,i-1
	      if(index(j).eq.num) goto 10
 50         continue
	  endif
	  index(i)=num
 100    continue
	return
	end
ccccc
ccccc
	function uniran(seed)
cc
cc  Draws a random number from the uniform distribution on [0,1].
cc
	real uniran
	integer seed
	integer quot
cc
	seed=seed*5761+999
	quot=seed/65536
	seed=seed-quot*65536
	uniran=float(seed)/65536.D0
	return
	end
ccccc
ccccc
	subroutine xrfshsort(a,n)
cc
cc  Sorts the array a of length n.
cc
	double precision a(n)
	double precision t
	integer gap
cc
	gap=n
 100    gap=gap/2
	if(gap.eq.0) goto 200
	do 180 i=1,n-gap
	  j=i
 120      if(j.lt.1) goto 180
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
 180    continue
	goto 100
 200    return
	end
ccccc
ccccc
	subroutine xrfishsort(a,kk)
cc
cc  Sorts the integer array a of length kk.
cc
	integer a(kk)
	integer t
	integer gap
cc
	gap=kk
 100    gap=gap/2
	if(gap.eq.0) goto 200
	do 180 i=1,kk-gap
	  j=i
 120      if(j.lt.1) goto 180
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
 180    continue
	goto 100
 200    return
	end
ccccc
ccccc
	subroutine xrfcovcopy(a,b,n1,n2)
cc
cc  Copies matrix a to matrix b.
cc
	double precision a(n1,n2)
	double precision b(n1,n2)
cc
	do 100 i=1,n1
	  do 90 j=1,n2
	    b(i,j)=a(i,j)
 90       continue
 100    continue
	return
	end
ccccc
ccccc
	subroutine xrfgenpn(n,nsel,index)
cc
cc    Constructs all subsets of nsel cases out of n cases.
cc
	integer index(nsel)
cc
	k=nsel
	index(k)=index(k)+1
 10     if(k.eq.1.or.index(k).le.(n-(nsel-k))) goto 100
	k=k-1
	index(k)=index(k)+1
	do 50 i=k+1,nsel
	  index(i)=index(i-1)+1
 50     continue
	goto 10
 100    return
	end
ccccc
ccccc
	function xrffindq(aw,ncas,k,index)
cc
cc  Finds the k-th order statistic of the array aw of length ncas.
cc
	double precision xrffindq
	double precision aw(ncas)
	double precision ax,wa
	integer index(ncas)
cc
	do 10 j=1,ncas
	  index(j)=j
 10     continue
	l=1
	lr=ncas
 20     if(l.ge.lr) goto 90
	ax=aw(k)
	jnc=l
	j=lr
 30     if(jnc.gt.j) goto 80
 40     if(aw(jnc).ge.ax) goto 50
	jnc=jnc+1
	goto 40
 50     if(aw(j).le.ax) goto 60
	j=j-1
	goto 50
 60     if(jnc.gt.j) goto 70
	i=index(jnc)
	index(jnc)=index(j)
	index(j)=i
	wa=aw(jnc)
	aw(jnc)=aw(j)
	aw(j)=wa
	jnc=jnc+1
	j=j-1
 70     goto 30
 80     if(j.lt.k) l=jnc
	if(k.lt.jnc) lr=j
	goto 20
 90     xrffindq=aw(k)
	return
	end
ccccc
ccccc
	subroutine xrfrdraw(a,n,seed,ntot,mini,ngroup,kmini)
cc
cc  Draws ngroup nonoverlapping subdatasets out of a dataset of size n,
cc  such that the selected case numbers are uniformly distributed from 1 to n.
cc
	integer a(2,ntot)
	integer mini(kmini)
	integer seed
cc
	jndex=0
	do 10 k=1,ngroup
	  do 20 m=1,mini(k)
	    nrand=int(uniran(seed)*(n-jndex))+1 
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
 6                  continue
		    a(1,i)=nrand+i-1
		    a(2,i)=k
		    goto 20
		  endif
 5              continue
	    endif
 20       continue
 10     continue
	return
	end
ccccc
ccccc
	function xrfodd(n)
cc
	logical xrfodd
cc
	xrfodd=.true.
	if(2*(n/2).eq.n) xrfodd=.false.
	return
	end
ccccc
ccccc
	function xrfnbreak(nhalf,n,nvar)
	integer xrfnbreak
cc
cc  Computes the breakdown value of the MCD estimator 
cc
	if (nhalf.lt.(n+nvar+1)/2) then
	  xrfnbreak=(nhalf-nvar+1)*100/n
	else
	  xrfnbreak=(n-nhalf+1)*100/n
	endif
	return
	end
ccccc
ccccc
        subroutine xrfmcduni(w,ncas,jqu,slutn,bstd,aw,aw2,factor,len)
cc
cc  mcduni : calculates the MCD in the univariate case.
cc           w contains the ordered observations 
cc
       implicit double precision (a-h,o-z), integer(i-n)
       double precision w(ncas),aw(ncas),aw2(ncas)
       double precision slutn(len)
cc
       sq=0.D0
       sqmin=0.D0
       ndup=1
       do 5 j=1,ncas-jqu+1
 5       slutn(j)=0.D0
       do 20 jint=1,ncas-jqu+1
       aw(jint)=0.D0
       do 10 j=1,jqu
       aw(jint)=aw(jint)+w(j+jint-1)
       if (jint.eq.1) sq=sq+w(j)*w(j)
 10    continue 
       aw2(jint)=aw(jint)*aw(jint)/jqu
       if (jint.eq.1) then 
	 sq=sq-aw2(jint)
         sqmin=sq
	 slutn(ndup)=aw(jint)
         len=jint
		      else 
         sq=sq-w(jint-1)*w(jint-1)+
     *        w(jint+jqu-1)*w(jint+jqu-1)
     *        -aw2(jint) + aw2(jint-1)
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
 20    continue 
       slutn(1)=slutn(int((ndup+1)/2))/jqu
       bstd=factor*sqrt(sqmin/jqu)
       return 
       end 
ccccc
