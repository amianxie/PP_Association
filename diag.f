c=======================================================================
c* Desk :Find the eignvalues and eignvactors: i.e. Modes:  
c !Modes start from NDD
c
      subroutine CharmDiagq(N,hessian,NDD,mode_need,ev,vec)
      implicit none
      INTEGER N,mode_need,NDD,NADD,i !N=3*npart
c      parameter (n_dim=225000000)
      real*8 hessian((N*N+N)/2),hes((N*N+N)/2)
      REAL*8 VEC(N,mode_need),EV(N)
      REAL*8 A(N+1),B(N+1),P(N+1),W(N+1),TA(N+1),TB(N+1),Y(N+1)
c
c      order=1
c      do i=1,N
c         do j=i,N
c            hessian(order)=hes(i,j)
c            order=order+1
c         enddo
c      enddo
      do i=1,(N*N+N)/2
         hes(i)=hessian(i)
      enddo
      NADD=NDD-1                !NDD is the first mode to find.
      call DIAGQ(N,hes,mode_need,VEC,A,B,P,W,EV,TA,TB,Y,NADD)
c
      return
      end
c
c-----------------------------------------------------------------------
c
c**********************************************************************
cc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      SUBROUTINE DIAGQ(NX,DD,NFRQX,VEC,A,B,P,W,EV,TA,TB,Y,NADD)
C
C     THIS ROUTINE IS A CONGLOMERATION OF GIVEN, HOUSEC, AND EIGEN
C     WHERE THE BEST FEATURES OF EACH WERE KEPT AND SEVERAL OTHER
C     MODIFICATIONS HAVE BEEN MADE TO INCREASE EFFICIENCY AND ACCURACY.
C
C   By Bernard R. Brooks   1981
C
C   NX      - ORDER OF MATRIX
C   NFRQX   - NUMBER OF ROOTS DESIRED
C   DD      - SECOND DERIVATIVE MATRIX IN UPPER TRIANGULAR FORM
C   VEC     - EIGENVECTORS RETURNED (NX,NFRQX)
C   EV      - EIGENVALUES RETURNED (NX)
C   A,B,P,W,TA,TB,Y - ALL SCRATCH VECTORS (NX+1)
C   NADD    - NUMBER OF LOWEST ROOTS TO SKIP
C
C
      implicit none
! This is from  include 'number.fcm'
CHARMM Element source/fcm/number.fcm 1.1
C
C This file contains floating point numbers.
C
C positive numbers
      REAL*8     ZERO, ONE, TWO, THREE, FOUR, FIVE, SIX,
     &           SEVEN, EIGHT, NINE, TEN, ELEVEN, TWELVE, THIRTN,
     &           FIFTN, NINETN, TWENTY, THIRTY
      PARAMETER (ZERO   =  0.D0, ONE    =  1.D0, TWO    =  2.D0,
     &           THREE  =  3.D0, FOUR   =  4.D0, FIVE   =  5.D0,
     &           SIX    =  6.D0, SEVEN  =  7.D0, EIGHT  =  8.D0,
     &           NINE   =  9.D0, TEN    = 10.D0, ELEVEN = 11.D0,
     &           TWELVE = 12.D0, THIRTN = 13.D0, FIFTN  = 15.D0,
     &           NINETN = 19.D0, TWENTY = 20.D0, THIRTY = 30.D0)
C
      REAL*8     FIFTY, SIXTY, SVNTY2, EIGHTY, NINETY, HUNDRD,
     &           ONE2TY, ONE8TY, THRHUN, THR6TY, NINE99, FIFHUN, THOSND,
     &           FTHSND,MEGA
      PARAMETER (FIFTY  = 50.D0,  SIXTY  =  60.D0,  SVNTY2 =   72.D0,
     &           EIGHTY = 80.D0,  NINETY =  90.D0,  HUNDRD =  100.D0,
     &           ONE2TY = 120.D0, ONE8TY = 180.D0,  THRHUN =  300.D0,
     &           THR6TY=360.D0,   NINE99 = 999.D0,  FIFHUN = 1500.D0,
     &           THOSND = 1000.D0,FTHSND = 5000.D0, MEGA   =   1.0D6)
C
C negative numbers
      REAL*8     MINONE, MINTWO, MINSIX
      PARAMETER (MINONE = -1.D0,  MINTWO = -2.D0,  MINSIX = -6.D0)
C
C common fractions
      REAL*8     TENM20,TENM14,TENM8,TENM5,PT0001,PT0005,PT001,PT005,
     &           PT01, PT02, PT05, PTONE, PT125, PT25, SIXTH, THIRD,
     &           PTFOUR, PTSIX, HALF, PT75, PT9999, ONEPT5, TWOPT4
      PARAMETER (TENM20 = 1.0D-20,  TENM14 = 1.0D-14,  TENM8  = 1.0D-8,
     &           TENM5  = 1.0D-5,   PT0001 = 1.0D-4, PT0005 = 5.0D-4,
     &           PT001  = 1.0D-3,   PT005  = 5.0D-3, PT01   = 0.01D0,
     &           PT02   = 0.02D0,   PT05   = 0.05D0, PTONE  = 0.1D0,
     &           PT125  = 0.125D0,  SIXTH  = ONE/SIX,PT25   = 0.25D0,
     &           THIRD  = ONE/THREE,PTFOUR = 0.4D0,  HALF   = 0.5D0,
     &           PTSIX  = 0.6D0,    PT75   = 0.75D0, PT9999 = 0.9999D0,
     &           ONEPT5 = 1.5D0,    TWOPT4 = 2.4D0)
C
C others
      REAL*8 ANUM,FMARK
      REAL*8 RSMALL,RBIG
      PARAMETER (ANUM=9999.0D0, FMARK=-999.0D0)
      PARAMETER (RSMALL=1.0D-10,RBIG=1.0D20)
C
C Machine constants (these are very machine dependent).
C
C RPRECI should be the smallest number you can add to 1.0 and get a number
C that is different than 1.0.  Actually: the following code must pass for
C for all real numbers A (where no overflow or underflow conditions exist).
C 
C         B = A * RPRECI
C         C = A + B
C         IF(C.EQ.A) STOP 'precision variable is too small'
C 
C The RBIGST value should be the smaller of:
C 
C         - The largest real value               
C         - The reciprocal of the smallest real value 
C
C If there is doubt, be conservative.
C 
C!!!!! NOTE:  Some of the values have not been checked....
C!!!!! Please fix these values - BRB .......
C
      REAL*8 RPRECI,RBIGST
C
      PARAMETER (RPRECI = 2.22045D-16, RBIGST = 4.49423D+307)
C
c End of 'number.fcm'
c
c This is from:  include 'stream.fcm'
CHARMM Element source/fcm/stream.fcm 1.1
C
C     This is the STREAM data block.
C     It contains information abount the current runstream.
C
C     MXSTRM - Maximum number of active stream files.
C     POUTU  - Default output unit number.
C     NSTRM  - Number of active input streams.
C     ISTRM = JSTRM(NSTRM) - Current input unit number
C     JSTRM(*) - stack of input streams numbers.
C     OUTU   - Unit number for all standard CHARMM output.
C     PRNLEV - Print level control for all writing to OUTU
C     IOLEV  - -1 to 1  -1=write no files.  1= write all files.
C     WRNLEV - -5 TO 10  0=SEVERE ONLY, 10=LIST ALL WARNINGS
C     LOWER  - if .true. all files with names not in double quotes
C              will be opened in lower case for write. For read
C              UPPER case will be tried first and if not succesful
C              file name will be converted to lower case.
C     QLONGL - Use long lines in the output where appropriate.
C              (Otherwise, restrict output lines to 80 characters)
C
      LOGICAL LOWER,QLONGL
      INTEGER MXSTRM,POUTU
      PARAMETER (MXSTRM=20,POUTU=6)
      INTEGER   NSTRM,ISTRM,JSTRM,OUTU,PRNLEV,WRNLEV,IOLEV
C
      COMMON /CASE/   LOWER, QLONGL
      COMMON /STREAM/ NSTRM,ISTRM,JSTRM(MXSTRM),OUTU,PRNLEV,WRNLEV,IOLEV
C
C End of 'stream.fcm'
c
      INTEGER NX,NFRQX,NADD
      REAL*8 VEC(*)
      REAL*8 A(NX),B(NX),P(NX),W(NX),EV(NX),TA(NX),TB(NX),Y(NX)
      Real*8 DD((NX*NX+NX)/2)      !For system with no more than 10000 points
C 
      REAL*8 ETA,THETA,DEL1,DELTA,SMALL,DELBIG,THETA1,TOLER,ETAR
      REAL*8 RPOWER,RPOW1,RAND1,FACTOR,ANORM,U,ANORMR
      REAL*8 SUM1,BX,S,SGN,TEMP,XKAP,EXPR,ALIMIT,ROOTL,ROOTX,TRIAL,F0
      REAL*8 AROOT,ELIM1,ELIM2,T,EPR,XNORM,XNORM1,EVDIFF
      INTEGER N,NEV,NEVADD,NTOT,I,IPT,J,IJ,NN,MI,MI1,JI,JI2,II
      INTEGER ML,ML1,L,M,K,MJ,MJ1,NOMTCH,NOM,IA,ITER
      INTEGER J1,MK,MK1,KK
C
CCC
      real*8 anumx
c      common /Eigen/ DD
CCC
      anumx=zero
      do i=1,nx
         A(I)=anumx
         B(I)=anumx
         P(I)=anumx 
         W(I)=anumx
         EV(I)=anumx
         TA(I)=anumx
         TB(I)=anumx
         Y(I)=anumx
      enddo
CCC
C
      ETA=RPRECI
      THETA=RBIGST
C
      N=NX
      NEV=NFRQX
      NEVADD=NEV+NADD
C
      DEL1=ETA/100.0
      DELTA=ETA**2*100.0
      SMALL=ETA**2/100.0
      DELBIG=THETA*DELTA/1000.0
      THETA1=1000.0/THETA
      TOLER=100.0*ETA
      ETAR=1.0/ETA
      RPOWER=8388608.0
      RPOW1=RPOWER*0.50
      RAND1=RPOWER-3.0
C
C Find largest element.
      FACTOR=ZERO
      NTOT=(N*(N+1))/2
      DO I=1,NTOT
         FACTOR=MAX(FACTOR,ABS(DD(I)))
      ENDDO
C
C Check for zero matrix.
      IF(FACTOR.LE.THETA1) THEN
         IF(WRNLEV.GE.2) WRITE(OUTU,811)
 811     FORMAT(' WARNING FROM <DIAGQ>. Zero matrix passed.',
     1     ' Identity matrix returned.')
         DO I=1,NEV
            EV(I)=ZERO
            IPT=(I-1)*N
            DO J=1,N
               IPT=IPT+1
               VEC(IPT)=ZERO
               IF(I+NADD.EQ.J) VEC(IPT)=ONE
            ENDDO
         ENDDO
         RETURN
      ENDIF
C
C Compute norm of matrix
      FACTOR=ONE/FACTOR
      IJ=0
      ANORM=ZERO
      DO I=1,N
         DO J=I,N
            IJ=IJ+1
            U=(DD(IJ)*FACTOR)**2
            IF(I.EQ.J) U=U*HALF
            ANORM=ANORM+U
         ENDDO
      ENDDO
C
C Scale the matrix
      ANORM=SQRT(ANORM+ANORM)/FACTOR
      ANORMR=ONE/ANORM
      DO I=1,NTOT
         DD(I)=DD(I)*ANORMR
      ENDDO
C
      NN=N-1
      MI=0
      MI1=N-1
C
C Perform trigiagonalization
      DO I=1,NN
         SUM1=ZERO
         B(I)=ZERO
         JI=I+1
         IPT=MI+I
         A(I)=DD(IPT)
         IPT=IPT+1
         BX=DD(IPT)
         JI2=JI+1
         DO J=JI2,N
            IPT=IPT+1
            SUM1=SUM1+DD(IPT)*DD(IPT)
         ENDDO
         IF(SUM1.LT.SMALL) THEN
            B(I)=BX
            DD(MI+JI)=ZERO
         ELSE
            S=SQRT(SUM1+BX**2)
            SGN=SIGN(ONE,BX)
            TEMP=ABS(BX)
            W(JI)=SQRT(HALF*(ONE+(TEMP/S)))
            IPT=MI+JI
            DD(IPT)=W(JI)
            II=I+2
            IF(II.LE.N) THEN
               TEMP=SGN/(TWO*W(JI)*S)
               DO J=II,N
                  IPT=IPT+1
                  W(J)=TEMP * DD(IPT)
                  DD(IPT)=W(J)
               ENDDO
            ENDIF
            B(I)=-SGN*S
C
            DO J=JI,N
               P(J)=ZERO
            ENDDO
            ML=MI + MI1
            ML1=MI1-1
            DO L=JI,N
               IPT=ML+L
               DO M=L,N
                  BX=DD(IPT)
                  P(L)=P(L)+BX*W(M)
                  IF(L.NE.M) P(M)=P(M)+BX*W(L)
                  IPT=IPT+1
               ENDDO
               ML=ML +ML1
               ML1=ML1-1
            ENDDO
C
C
            XKAP=ZERO
            DO K=JI,N
               XKAP=XKAP+W(K)*P(K)
            ENDDO
            DO L=JI,N
               P(L)=P(L)-XKAP*W(L)
            ENDDO
            MJ=MI+MI1
            MJ1=MI1-1
            DO J=JI,N
               DO K=J,N
                  EXPR=(P(J)*W(K))+(P(K)*W(J))
                  DD(MJ+K)=DD(MJ+K)-EXPR-EXPR
               ENDDO
               MJ=MJ+MJ1
               MJ1=MJ1-1
            ENDDO
         ENDIF
         MI=MI+MI1
         MI1=MI1-1
      ENDDO
C
C Begin sturm bisection method.
C
      A(N)=DD(MI+N)
      B(N)=ZERO
C
      ALIMIT=ONE
      DO I=1,N
         W(I)=B(I)
         B(I)=B(I)*B(I)
      ENDDO
      DO I=1,NEVADD
         EV(I)=ALIMIT
      ENDDO
      ROOTL=-ALIMIT
C
      DO I=1,NEVADD
         ROOTX=ALIMIT
         DO J=I,NEVADD
            ROOTX=MIN(ROOTX,EV(J))
         ENDDO
         EV(I)=ROOTX
C
 130     CONTINUE
            TRIAL=(ROOTL+EV(I))*HALF
c##IF CRAY (oldcode_for_cray)
c            IF(TRIAL.EQ.ROOTL.OR.TRIAL.EQ.EV(I)) GOTO 200
c##ELSE (oldcode_for_cray)
            EVDIFF=ABS(ROOTL-EV(I))
            IF(EVDIFF.LT.THETA1) GOTO 200
            IF(EVDIFF*ETAR.LT.ABS(TRIAL)) GOTO 200
c##ENDIF (oldcode_for_cray)
            NOMTCH=N
            J=1
 150     CONTINUE
            F0=A(J)-TRIAL
 160     CONTINUE
            IF(ABS(F0).LT.THETA1) GOTO 170
            IF(F0.GE.ZERO) NOMTCH=NOMTCH-1
            J=J+1
            IF(J.GT.N) GOTO 180
            F0=A(J)-TRIAL-B(J-1)/F0
            GOTO160
 170     CONTINUE
            J=J+2
            NOMTCH=NOMTCH-1
            IF(J.LE.N) GOTO 150
 180     CONTINUE
            IF(NOMTCH.GE.I) GOTO 190
            ROOTL=TRIAL
            GOTO 130
 190     CONTINUE
            EV(I)=TRIAL
            NOM=MIN(NEVADD,NOMTCH)
            EV(NOM)=TRIAL
            GOTO 130
 200     CONTINUE
      ENDDO
C
C Finished computing requested eigenvalues
      DO I=1,NEV
         EV(I)=EV(I+NADD)
      ENDDO
C
C Compute eigenvectors (backtransformation)
      DO I=1,NEV
         AROOT=EV(I)
         DO J=1,N
            Y(J)=ONE
         ENDDO
         IA=IA+1
         IF(I.EQ.1) THEN
            IA=0
         ELSE
            IF(ABS(EV(I-1)-AROOT).GE.TOLER) IA=0
         ENDIF
         ELIM1=A(1)-AROOT
         ELIM2=W(1)
         DO J=1,NN
            IF(ABS(ELIM1).LE.ABS(W(J))) THEN
               TA(J)=W(J)
               TB(J)=A(J+1)-AROOT
               P(J)=W(J+1)
               TEMP=ONE
               IF(ABS(W(J)).GT.THETA1) TEMP=ELIM1/W(J)
               ELIM1=ELIM2-TEMP*TB(J)
               ELIM2=-TEMP*W(J+1)
            ELSE
               TA(J)=ELIM1
               TB(J)=ELIM2
               P(J)=ZERO
               TEMP=W(J)/ELIM1
               ELIM1=A(J+1)-AROOT-TEMP*ELIM2
               ELIM2=W(J+1)
            ENDIF
            B(J)=TEMP
         ENDDO
C
         TA(N)=ELIM1
         TB(N)=ZERO
         P(N)=ZERO
         P(NN)=ZERO
         ITER=1
         IF(IA.NE.0) GOTO 460
C
 320     L=N+1
         DO J=1,N
            L=L-1
 330        CONTINUE
            IF(L.EQ.N) THEN
               ELIM1=Y(L)
            ELSE IF(L.EQ.N-1) THEN
               ELIM1=Y(L)-Y(L+1)*TB(L)
            ELSE
               ELIM1=Y(L)-Y(L+1)*TB(L)-Y(L+2)*P(L)
            ENDIF
C
C Overflow check
            IF(ABS(ELIM1).GT.DELBIG) THEN
               DO K=1,N
                  Y(K)=Y(K)/DELBIG
               ENDDO
               GOTO 330
            ENDIF
            TEMP=TA(L)
            IF(ABS(TEMP).LT.DELTA) TEMP=DELTA
            Y(L)=ELIM1/TEMP
         ENDDO
C
         IF(ITER.EQ.2) GOTO 500
         ITER=ITER+1
C
 420     CONTINUE
         ELIM1=Y(1)
         DO J=1,NN
            IF(TA(J).EQ.W(J)) THEN
               Y(J)=Y(J+1)
               ELIM1=ELIM1-Y(J+1)*B(J)
            ELSE
               Y(J)=ELIM1
               ELIM1=Y(J+1)-ELIM1*B(J)
            ENDIF
         ENDDO
         Y(N)=ELIM1
         GOTO 320
C
 460     CONTINUE
         DO J=1,N
            RAND1=MOD(4099.0*RAND1,RPOWER)
            Y(J)=RAND1/RPOW1-ONE
         ENDDO
         GOTO 320
C
C Orthog to previous
 500     IF(IA.EQ.0) GOTO 550
         DO J1=1,IA
            K=I-J1
            TEMP=ZERO
            IPT=(K-1)*N
            DO J=1,N
               IPT=IPT+1
               TEMP=TEMP+Y(J)*VEC(IPT)
            ENDDO
            IPT=(K-1)*N
            DO J=1,N
               IPT=IPT+1
               Y(J)=Y(J)-TEMP*VEC(IPT)
            ENDDO
          ENDDO
 550      CONTINUE
          IF(ITER.EQ.1) GOTO 420
C
C Normalize
 560     CONTINUE
         ELIM1=ZERO
         DO J=1,N
            ELIM1=MAX(ELIM1,ABS(Y(J)))
         ENDDO
         TEMP=ZERO
         DO J=1,N
            ELIM2=Y(J)/ELIM1
            TEMP=TEMP+ELIM2*ELIM2
         ENDDO
         TEMP=ONE/(SQRT(TEMP)*ELIM1)
         DO J=1,N
            Y(J)=Y(J)*TEMP
            IF(ABS(Y(J)).LT.DEL1) Y(J)=ZERO
         ENDDO
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=Y(J)
         ENDDO
      ENDDO
C
      DO I=1,NEV
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            Y(J)=VEC(IPT)
         ENDDO
C
         L=N-2
         MK=(N*(N-1))/2-3
         MK1=3
C
         DO J=1,L
            T=ZERO
            K=N-J-1
            M=K+1
            DO KK=M,N
               T=T+DD(MK+KK)*Y(KK)
            ENDDO
            DO KK=M,N
               EPR=T*DD(MK+KK)
               Y(KK)=Y(KK)-EPR-EPR
            ENDDO
            MK=MK-MK1
            MK1=MK1+1
         ENDDO 
C
         T=ZERO
         DO J=1,N
            T=T+Y(J)*Y(J)
         ENDDO
         XNORM=SQRT(T)
         XNORM1=ONE/XNORM
         DO J=1,N
            Y(J)=Y(J)*XNORM1
         ENDDO
C
         IPT=(I-1)*N
         DO J=1,N
            IPT=IPT+1
            VEC(IPT)=Y(J)
         ENDDO
      ENDDO
C
      DO I=1,N
         EV(I)=EV(I)*ANORM
      ENDDO
C
      RETURN
      END
C*Desk :CharmDiagq: End!
c-------------------------------------------------------------------------------------
