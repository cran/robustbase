C=======================================================================
      SUBROUTINE rlSTORm2(Y,N,J,YJ)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(N)
C-----------------------------------------------------------------------
C     rlSTORm2 SEARCHES THE J-TH VALUE IN ORDER OF MAGNITUDE IN
C     A VECTOR OF LENGTH N.
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobmm.f -------------------------
      L=1
      LR=N
 20   IF (L.GE.LR) GOTO 90
      AX=Y(J)
      JNC=L
      JJ=LR
 30   IF(JNC.GT.JJ) GOTO 80
 40   IF (Y(JNC).GE.AX) GOTO 50
      JNC=JNC+1
      GOTO 40
 50   IF(Y(JJ).LE.AX) GOTO 60
      JJ=JJ-1
      GOTO 50
 60   IF(JNC.GT.JJ) GOTO 70
      WA=Y(JNC)
      Y(JNC)=Y(JJ)
      Y(JJ)=WA
      JNC=JNC+1
      JJ=JJ-1
 70   GOTO 30
 80   IF(JJ.LT.J) L=JNC
      IF(J.LT.JNC) LR=JJ
      GOTO 20
 90   YJ=Y(J)
      RETURN
      END
C=======================================================================
      SUBROUTINE rlCOLbi(V1,V2,MLT,M,IOUT)
C.......................................................................
      DOUBLE PRECISION V1(M),V2(M),MLT
C-----------------------------------------------------------------------
C     AUXILIARY ROUTINE FOR rlLARSbi
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobbi.f -------------------------
      DO 220 I=1,M
         IF (I .EQ. IOUT) GOTO 220
         V1(I)=V1(I)-V2(I)*MLT
 220  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE rlICHGbi(A,B)
C.......................................................................
C     AUXILIARY ROUTINE FOR rlLARSbi
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobbi.f -------------------------
      DOUBLE PRECISION A,B,C
      C=A
      A=B
      B=C
      RETURN
      END
C=======================================================================
      SUBROUTINE rlLARSbi(X,Y,N,NP,MDX,MDT,TOL,NIT,K,
     +     KODE,SIGMA,THETA,RS,SC1,SC2,SC3,SC4,BET0)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(MDX,NP),Y(N),THETA(MDT),RS(N),SC1(N),SC2(NP),
     +     SC3(NP),SC4(NP)
      INTEGER OUT
      LOGICAL STAGE,TEST
      DATA ZERO,TWO,EPS,BIG/0.D0,2.D0,1.0D-10,3.401D38/
C      DATA ZERO,TWO,EPS,BIG/0.D0,2.D0,2.22D-16,1.796D308/
C-----------------------------------------------------------------------
C     LEAST ABSOLUTE RESIDUALS -- aka  L_1 - Regression
C      --> Result in THETA[1:NP]
C-----------------------------------------------------------------------
C--- copied from robust package: src/lmrobbi.f -------------------------
      SUM=ZERO
      DO 10 J=1,NP
         SC4(J)=DBLE(J)
         SC2(J)=ZERO
 10   CONTINUE
      DO 40 I=1,N
         SC1(I)=DBLE(NP+I)
         THETA(I)=Y(I)
         IF (Y(I) .GE. ZERO) GOTO 30
         DO 20 J=1,NP
            X(I,J)=-X(I,J)
 20      CONTINUE
         THETA(I)=-THETA(I)
         SC1(I)=-SC1(I)
 30      SUM=SUM+THETA(I)
 40   CONTINUE
C-----------------------------------------------------------------------
C     COMPUTE THE MARGINAL COSTS.
C-----------------------------------------------------------------------
      SUMIN=SUM
      DO 60 J=1,NP
         SUM=ZERO
         DO 50 I=1,N
            SUM=SUM+X(I,J)
 50      CONTINUE
         SC3(J)=SUM
 60   CONTINUE
C-----------------------------------------------------------------------
C     STAGE I. DETERMINE THE VECTOR TO ENTER THE BASIS.
C-----------------------------------------------------------------------
      STAGE=.TRUE.
      KOUNT=0
      KR=1
      KL=1
 70   VMAX=-1.D0
      DNP=DBLE(NP)
      DO 80 J=KR,NP
         IF (DABS(SC4(J)) .GT. DNP) GOTO 80
         D=DABS(SC3(J))
         IF (D-VMAX .LE. ZERO) GOTO 80
         IF (D-VMAX .LE. EPS)  GOTO 80
         VMAX=D
         IN=J
 80   CONTINUE
      IF (SC3(IN) .GE. ZERO) GOTO 100
      DO 90 I=1,N
         X(I,IN)=-X(I,IN)
 90   CONTINUE
      SC3(IN)=-SC3(IN)
      SC4(IN)=-SC4(IN)
C-----------------------------------------------------------------------
C     DETERMINE THE VECTOR TO LEAVE THE BASIS.
C-----------------------------------------------------------------------
 100  K=0
      DO 110 I=KL,N
         D=X(I,IN)
         IF (D .LE. TOL) GOTO 110
         K=K+1
         Y(K)=THETA(I)/D
         RS(K)=DBLE(I)
         TEST=.TRUE.
 110  CONTINUE
 120  IF (K .GT. 0) GOTO 130
      TEST=.FALSE.
      GOTO 150
 130  VMIN=BIG
      DO 140 I=1,K
         IF (Y(I)-VMIN .GE. ZERO) GOTO 140
         IF (VMIN-Y(I) .LE. EPS)  GOTO 140
         J=I
         VMIN=Y(I)
         OUT=INT(RS(I))
 140  CONTINUE
      Y(J)=Y(K)
      RS(J)=RS(K)
      K=K-1
C-----------------------------------------------------------------------
C     CHECK FOR LINEAR DEPENDENCE IN STAGE I.
C-----------------------------------------------------------------------
 150  IF (TEST .OR. .NOT.STAGE) GOTO 170
      DO 160 I=1,N
         CALL rlICHGbi(X(I,KR),X(I,IN))
 160  CONTINUE
      CALL rlICHGbi(SC3(KR),SC3(IN))
      CALL rlICHGbi(SC4(KR),SC4(IN))
      KR=KR+1
      GOTO 260
 170  IF (TEST) GOTO 180
      KODE=2
      GOTO 350
 180  PIVOT=X(OUT,IN)
      IF (SC3(IN)-PIVOT-PIVOT .LE. TOL) GOTO 200
      DO 190 J=KR,NP
         D=X(OUT,J)
         SC3(J)=SC3(J)-D-D
         X(OUT,J)=-D
 190  CONTINUE
      D=THETA(OUT)
      SUMIN=SUMIN-D-D
      THETA(OUT)=-D
      SC1(OUT)=-SC1(OUT)
      GOTO 120
C-----------------------------------------------------------------------
C     PIVOT ON X(OUT,IN).
C-----------------------------------------------------------------------
 200  DO 210 J=KR,NP
         IF (J.EQ.IN) GOTO 210
         X(OUT,J)=X(OUT,J)/PIVOT
 210  CONTINUE
      THETA(OUT)=THETA(OUT)/PIVOT
      DO 230 J=KR,NP
         IF (J .EQ. IN) GOTO 230
         D=X(OUT,J)
         SC3(J)=SC3(J)-D*SC3(IN)
         CALL rlCOLbi(X(1,J),X(1,IN),D,N,OUT)
 230  CONTINUE
      SUMIN=SUMIN-SC3(IN)*THETA(OUT)
      DO 240 I=1,N
         IF (I .EQ. OUT) GOTO 240
         D=X(I,IN)
         THETA(I)=THETA(I)-D*THETA(OUT)
         X(I,IN)=-D/PIVOT
 240  CONTINUE
      SC3(IN)=-SC3(IN)/PIVOT
      X(OUT,IN)=1.D0/PIVOT
      CALL rlICHGbi(SC1(OUT),SC4(IN))
      KOUNT=KOUNT+1
      IF (.NOT.STAGE) GOTO 270
C-----------------------------------------------------------------------
C     INTERCHANGE ROWS IN STAGE I.
C-----------------------------------------------------------------------
      KL=KL+1
      DO 250 J=KR,NP
         CALL rlICHGbi(X(OUT,J),X(KOUNT,J))
 250  CONTINUE
      CALL rlICHGbi(THETA(OUT),THETA(KOUNT))
      CALL rlICHGbi(SC1(OUT),SC1(KOUNT))
 260  IF (KOUNT+KR .NE. NP+1) GOTO 70
C-----------------------------------------------------------------------
C     STAGE II. DETERMINE THE VECTOR TO ENTER THE BASIS.
C-----------------------------------------------------------------------
      STAGE=.FALSE.
 270  VMAX=-BIG
      DO 290 J=KR,NP
         D=SC3(J)
         IF (D .GE. ZERO) GOTO 280
         IF (D+TWO .GT. ZERO) GOTO 290
         D=-D-TWO
 280     IF (D-VMAX .LE. ZERO) GOTO 290
         IF (D-VMAX .LE. EPS)  GOTO 290
         VMAX=D
         IN=J
 290  CONTINUE
      IF (VMAX .LE. TOL) GOTO 310
      IF (SC3(IN) .GT. ZERO) GOTO 100
      DO 300 I=1,N
         X(I,IN)=-X(I,IN)
 300  CONTINUE
      SC3(IN)=-SC3(IN)-2.D0
      SC4(IN)=-SC4(IN)
      GOTO 100
C-----------------------------------------------------------------------
C     PREPARE OUTPUT
C-----------------------------------------------------------------------
 310  L=KL-1
      DO 330 I=1,N
         RS(I)=ZERO
         IF (I .GT. L .OR. THETA(I) .GE. ZERO) GOTO 330
         DO 320 J=KR,NP
            X(I,J)=-X(I,J)
 320     CONTINUE
         THETA(I)=-THETA(I)
         SC1(I)=-SC1(I)
 330  CONTINUE
      KODE=0
      IF (KR .NE. 1) GOTO 350
      DO 340 J=1,NP
         D=DABS(SC3(J))
         IF (D .LE. TOL .OR. TWO-D .LE. TOL) GOTO 350
 340  CONTINUE
      KODE=1
 350  DO 380 I=1,N
         K=INT(SC1(I))
         D=THETA(I)
         IF (K .GT. 0) GOTO 360
         K=-K
         D=-D
 360     IF (I .GE. KL) GOTO 370
         SC2(K)=D
         GOTO 380
 370     K=K-NP
         RS(K)=D
 380  CONTINUE
      K=NP+1-KR
      SUM=ZERO
      DO 390 I=KL,N
         SUM=SUM+THETA(I)
 390  CONTINUE
      SUMIN=SUM
      NIT=KOUNT
      DO 400 J=1,NP
         THETA(J)=SC2(J)
 400  CONTINUE
      DO 500 I=1,N
         Y(I)=DABS(RS(I))
 500  CONTINUE
      N2=N/2+1
      CALL RLSTORM2(Y,N,N2,SIGMA)
      SIGMA=SIGMA/BET0
      RETURN
      END
