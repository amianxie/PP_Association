CHARMM Element source/manip/rotlsq.src 1.1
      SUBROUTINE rotlsq(X1,Y1,Z1,NATOM1,X2,Y2,Z2,NATOM2,ATOMPR,NPAIR)
ccc     $                  LMASS,AMASS1,AMASS2,QWGHT,KWMULT,LNOROT,LPRINT)
C-----------------------------------------------------------------------
C     THE PROGRAM ROTATES COORDINATE SET 2 RESULTING IN A 2 SUCH THAT
C     THE SUM OF THE SQUARE OF THE DISTANCE BETWEEN EACH COORDINATE IN 1
C     AND 2 IS A MINIMUM.
C     THE LEAST SQUARE MINIMIZATION IS DONE ONLY WITH RESPECT TO THE
C     ATOMS REFERRED TO IN THE PAIR ARRAY. THE ROTATION MATRIX THAT IS
C     CALCULATED IS APPLIED TO THE ENTIRE SECOND COORDINATE SET.
C     BERNARD R. BROOKS
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/exfunc.fcm'
c##INCLUDE '~/charmm_fcm/number.fcm'
c##INCLUDE '~/charmm_fcm/stack.fcm'
C
      INTEGER NATOM1,NATOM2,NPAIR
      REAL*8 X1(*),Y1(*),Z1(*),X2(*),Y2(*),Z2(*)
c      REAL*8 AMASS1(*),AMASS2(*)
c      REAL*8 KWMULT(*)
c      LOGICAL LMASS,QWGHT,LNOROT,LPRINT
      INTEGER ATOMPR(2,*)
c      INTEGER OLDLST,MASS
C
c      OLDLST=LSTUSD
c      MASS=ALLSTK(IREAL8(NPAIR))
c      CALL PKMASS(AMASS1,AMASS2,STACK(MASS),ATOMPR,NPAIR,LMASS,
c     $     QWGHT,KWMULT)
c      print*,"@"
c      INTEGER i,j
c      DO j=1,NATOM1
c         write(*,2000) j,NATOM1,NATOM2,NPAIR,X1(j),Y1(j),Z1(j)
c         write(*,2000) j,NATOM1,NATOM2,NPAIR,X2(j),Y2(j),Z2(j)
c      ENDDO
c      DO i=1,2
c         DO j=1,NPAIR
c            write(*,1000) i,j,ATOMPR(i,j)
c         ENDDO
c      ENDDO
c 1000 format(3I3)     
c 2000 format(4I3,3F8.3)
 
      CALL ROTLS1(X1,Y1,Z1,X2,Y2,Z2,NATOM2,ATOMPR,NPAIR)
c     $     STACK(MASS),LPRINT,LNOROT)
c      CALL FRESTK(LSTUSD-OLDLST)
c
c      DO j=1,NATOM1
c         write(*,2000) j,NATOM1,NATOM2,NPAIR,X1(j),Y1(j),Z1(j)
c         write(*,2000) j,NATOM1,NATOM2,NPAIR,X2(j),Y2(j),Z2(j)
c      ENDDO


      RETURN
      END

      SUBROUTINE ROTLS1(XA,YA,ZA,XB,YB,ZB,NATOMB,ATOMPR,NPAIR)
c     $                  BMASS,LPRINT,LNOROT)
C-----------------------------------------------------------------------
C     THIS ROUTINE DOES THE ACTUAL ROTATION OF B TO MATCH A. THE NEW
C     ARRAY B IS RETURNED IN X.
C     THIS ROUTINE WRITTEN BY B. BROOKS , ADAPTED FROM ACTA CRYST
C     (1976) A32,922 W. KABSCH
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/number.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
C
      INTEGER NATOMB,NPAIR
      REAL*8 BMASS(NPAIR)
      REAL*8 XA(*),YA(*),ZA(*),XB(*),YB(*),ZB(*)
      INTEGER ATOMPR(2,NPAIR)
c      LOGICAL LPRINT,LNOROT
C
c##INCLUDE '~/charmm_fcm/corman.fcm'
C
      REAL*8 R(9),U(9),EVA(3),DEVA(3,3)
      REAL*8 RN(3), PHI
      REAL*8 CMXA,CMYA,CMZA,CMXB,CMYB,CMZB,CMXC,CMYC,CMZC
      REAL*8 XI,YI,ZI,XJ,YJ,ZJ
      REAL*8 TMASS,RMST,RMSV,A0
      INTEGER K,KA,KB,I,j
c      LOGICAL LPRINT2
C
c      LPRINT2=(LPRINT.AND.(PRNLEV.GT.6))
C
c      print*,"@@"

c      DO j=1,NPAIR
c         write(*,2000) j,NPAIR,XA(j),YA(j),ZA(j)
c         write(*,2000) j,NPAIR,XB(j),YB(j),ZB(j)
c      ENDDO
c 2000 format(2I3,3F8.3)
      PHI=0
      do i=1,9
         U(i)=0
      enddo
      do i=1,3
         EVA(i)=0
         RN(i)=0
         do j=1,3
            DEVA(i,j)=0
         enddo
      enddo

c      print*,xa(1),xb(1)

      do j=1,npair
        bmass(j)=1
      enddo


      CMXA=0.0
      CMYA=0.0
      CMZA=0.0
      CMXB=0.0
      CMYB=0.0
      CMZB=0.0
      TMASS=0.0
      DO K=1,NPAIR
        KA=ATOMPR(1,K)
        KB=ATOMPR(2,K)
C        IF (XA(KA).NE.ANUM .AND. XB(KB).NE.ANUM) THEN
          CMXA=CMXA+XA(KA)*BMASS(K)
          CMYA=CMYA+YA(KA)*BMASS(K)
          CMZA=CMZA+ZA(KA)*BMASS(K)
          CMXB=CMXB+XB(KB)*BMASS(K)
          CMYB=CMYB+YB(KB)*BMASS(K)
          CMZB=CMZB+ZB(KB)*BMASS(K)
          TMASS=TMASS+BMASS(K)
C        ENDIF
      ENDDO
      CMXA=CMXA/TMASS
      CMYA=CMYA/TMASS
      CMZA=CMZA/TMASS
      CMXB=CMXB/TMASS
      CMYB=CMYB/TMASS
      CMZB=CMZB/TMASS
      DO K=1,NATOMB
c        IF(XB(K).NE.ANUM) THEN
          XB(K)=XB(K)-CMXB
          YB(K)=YB(K)-CMYB
          ZB(K)=ZB(K)-CMZB
c        ENDIF
      ENDDO
C
      CMXC=CMXA-CMXB
      CMYC=CMYA-CMYB
      CMZC=CMZA-CMZB
c      CALL SETMSR('XCEN',CMXA)
c      CALL SETMSR('YCEN',CMYA)
c      CALL SETMSR('ZCEN',CMZA)
c      CALL SETMSR('XMOV',CMXC)
c      CALL SETMSR('YMOV',CMYC)
c      CALL SETMSR('ZMOV',CMZC)
c      IF (LPRINT) THEN
c        WRITE(OUTU,44) CMXB,CMYB,CMZB
c        WRITE(OUTU,45) CMXA,CMYA,CMZA
c        WRITE(OUTU,46) CMXC,CMYC,CMZC
c      ENDIF
c 44   FORMAT(' CENTER OF ATOMS BEFORE TRANSLATION',3F12.5)
c 45   FORMAT(' CENTER OF REFERENCE COORDINATE SET',3F12.5)
c 46   FORMAT(' NET TRANSLATION OF ROTATED ATOMS  ',3F12.5)
C
c      IF (LNOROT) THEN
C
C       USE A UNIT ROTATION MATRIX. NO ROTATION IS SPECIFIED
C
c        DO K=1,NATOMB
c          IF (XB(K).NE.ANUM) THEN
c            XB(K)=XB(K)+CMXA
c            YB(K)=YB(K)+CMYA
c            ZB(K)=ZB(K)+CMZA
c          ENDIF
c        ENDDO
C
c      ELSE
C
C       COMPUTE ROTATION MATRIX FROM LAGRANGIAN
C

        DO I=1,9
          R(I)=0.0
        ENDDO
        DO K=1,NPAIR
          KA=ATOMPR(1,K)
          KB=ATOMPR(2,K)
c          IF (XA(KA).NE.ANUM .AND. XB(KB).NE.ANUM) THEN
            XI=XB(KB)*BMASS(K)
            YI=YB(KB)*BMASS(K)
            ZI=ZB(KB)*BMASS(K)
            XJ=XA(KA)-CMXA
            YJ=YA(KA)-CMYA
            ZJ=ZA(KA)-CMZA
            R(1)=R(1)+XI*XJ
            R(2)=R(2)+XI*YJ
            R(3)=R(3)+XI*ZJ
            R(4)=R(4)+YI*XJ
            R(5)=R(5)+YI*YJ
            R(6)=R(6)+YI*ZJ
            R(7)=R(7)+ZI*XJ
            R(8)=R(8)+ZI*YJ
            R(9)=R(9)+ZI*ZJ
c          ENDIF
        ENDDO
C

        A0=0
        CALL FROTU(R,EVA,DEVA,U,A0)
c        ,LPRINT2)


C
c        do i=1,natomb
c        print*,XB(i),YB(i),ZB(i)
c       enddo
c        stop         
       
        DO K=1,NATOMB
c          IF (XB(K).NE.ANUM) THEN
            CMXC=U(1)*XB(K)+U(4)*YB(K)+U(7)*ZB(K)+CMXA
            CMYC=U(2)*XB(K)+U(5)*YB(K)+U(8)*ZB(K)+CMYA
            ZB(K)=U(3)*XB(K)+U(6)*YB(K)+U(9)*ZB(K)+CMZA
            XB(K)=CMXC
            YB(K)=CMYC
c          ENDIF
        ENDDO
C
c        IF (LPRINT) WRITE(OUTU,55) U
c 55     FORMAT(' ROTATION MATRIX',3(/1X,3F12.6))
        CALL FNDROT(U,RN,PHI)

c      print*,U(1)
c      stop
c        ,LPRINT)
c        CALL SETMSR('THET',PHI)
c        QAXISC=.TRUE.
c        AXISCX= ZERO
c        AXISCY= ZERO
c        AXISCZ= ZERO
c        AXISR = ONE
c        AXISX = RN(1)
c        AXISY = RN(2)
c        AXISZ = RN(3)
c        CALL SETMSR('XAXI',AXISX)
c        CALL SETMSR('YAXI',AXISY)
c        CALL SETMSR('ZAXI',AXISZ)
c        CALL SETMSR('RAXI',AXISR)
c        CALL SETMSR('XCEN',AXISCX)
c        CALL SETMSR('YCEN',AXISCY)
c        CALL SETMSR('ZCEN',AXISCZ)
c      ENDIF
C
      RMST=0.0
      DO K=1,NPAIR
        KA=ATOMPR(1,K)
        KB=ATOMPR(2,K)
        IF (XA(KA).NE.9999.AND. XB(KB).NE.9999) RMST=RMST+
     $       BMASS(K)*((XB(KB)-XA(KA))**2+(YB(KB)-YA(KA))**2+
     $       (ZB(KB)-ZA(KA))**2)
      ENDDO
      RMSV=SQRT(RMST/TMASS)
c      CALL SETMSR('RMS ',RMSV)
c      IF(LPRINT) WRITE(OUTU,14) RMST,TMASS,RMSV
c 14   FORMAT(' TOTAL SQUARE DIFF IS',F12.4,'  DENOMINATOR IS',F12.4,/
c     $    '       THUS RMS DIFF IS',F12.6)
C



      RETURN
      END

      SUBROUTINE FROTU(R2,EVA,DEVA,U,EVWID)
c      ,QPRINT)
C-----------------------------------------------------------------------
C     THIS ROUTINE SOLVES THE CONSTRAINED MINIMIZATION EQUATION
C     USING LAGRANGE MULTIPLIERS.
C     BERNARD R. BROOKS
C
C Input:
C     R(3,3)  - R matrix    ( R(1,3)=sum(Xai*Zbi*MASSi) )
C     QPRINT  - Flag indicating a verbose output
C     EVWID
C Output:
C     U(3,3)  - Resultant rotation matrix.
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/number.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
C
      REAL*8 R2(9),EVA(3),DEVA(3,3),U(9),EVWID
c      LOGICAL QPRINT
C
      real*8 R(3,3),U2(3,3),u3(3,3),A2(3,3)
      REAL*8 W(9),A(3,3),B(3,3),SCR(24),RTA(3,3)
      REAL*8 EV(3),RT(3,3),BT(3,3),S(3,3),T(3,3),DSW(3,3)
      real*8 w2(3,3),tem1(3),tem2(3),B1(3),B2(3),B3(3)
      REAL*8 NORM1,NORM2,NORM3
C
      INTEGER I,IPT,J,K,JP,JQ,KP,KQ
      REAL*8 TRACE,EVS,DET,DOT,DRAT,SWVAL,DSWVAL
      
      real*8 wyh1(3),wyh2(3),wyh3(3),wyh4(3),wyh5(3)
      real*8 wyh6(3),wyh7(3)
      real*8 a0
C
c      IF(QPRINT) write(6,56) ' The R matrix:',R
c  56  format(A/,(10X,3F12.5))
C
c      print*,"@@@"

       do i=1,3
          tem1(i)=0
          tem2(i)=0
          B1(i)=0
          B2(i)=0
          B3(i)=0
          EV(i)=0
          do j=1,3
             R(i,j)=0
             U2(i,j)=0
             U3(i,j)=0
             A2(i,j)=0
             A(i,j)=0
             B(i,j)=0
             RTA(i,j)=0
             RT(i,j)=0
             BT(i,j)=0
             S(i,j)=0
             T(i,j)=0
             DSW(i,j)=0
             W2(i,j)=0
          enddo
       enddo
       
       do i=1,9
          w(i)=0
       enddo

       do i=1,24
          scr(i)=0
       enddo


       R(1,1)=R2(1)
       R(1,2)=R2(2)
       R(1,3)=R2(3)
       R(2,1)=R2(4)
       R(2,2)=R2(5)
       R(2,3)=R2(6)
       R(3,1)=R2(7)
       R(3,2)=R2(8)
       R(3,3)=R2(9)


      CALL DETM33(R, DET)
c      IF(QPRINT) write(6,56) ' The Rdet value:',DET
C
      IPT=0
      DO I=1,3
        DO J=1,3
          IPT=IPT+1
          W(IPT)=0
          W2(I,J)=0
          DO K=1,3
            W(IPT)=W(IPT)+R(J,K)*R(I,K)
            W2(I,J)=W2(I,J)+R(K,I)*R(K,J)
          ENDDO
        ENDDO
      ENDDO

c       do i=1,3
c       print*,W2(i,1),W2(i,2),W2(i,3)
c       enddo

C
c      IF(QPRINT) write(6,56) ' The W matrix:',W
C
C Handle the case where just one atom is selected.
      TRACE=W2(1,1)+W2(2,2)+W2(3,3)
      IF(TRACE.LT.0.000001) THEN
        CALL VZERO(EVA,3)
        CALL VZERO(DEVA,9)
        CALL VZERO(U,9)
c        DO I=1,3
c          U(I,I)=ONE
c        ENDDO
         U(1)=1
         U(5)=1
         U(9)=1
        RETURN
      ENDIF
C
c      W2(1,1)=W(1)
c      W2(1,2)=W(2)
c      W2(1,3)=W(3)
c      W2(2,1)=W(4)
c      W2(2,2)=W(5)
c      W2(2,3)=W(6)
c      W2(3,1)=W(7)
c      W2(3,2)=W(8)
c      W2(3,3)=W(9)

c      CALL DIAGQ(3,3,W,A2,wyh1,wyh2,wyh3,wyh4,EV,
c     $           wyh5,wyh6,wyh7,0)

      a0=0

c      call rs(3,3,w2,EV,3,A,tem1,tem2,a0)
      W(1)=W2(1,1)
      W(2)=W2(1,2)
      W(3)=W2(1,3)
      W(4)=W2(2,2)
      W(5)=W2(2,3)
      W(6)=W2(3,3)
      call CharmDiagq(3,W,1,3,EV,A)
c      call CharmDiagq(3,3,W,EV,A)
c       A(1,1)=A2(1)
c       A(1,2)=A2(2)
c       A(1,3)=A2(3)
c       A(2,1)=A2(4)
c       A(2,2)=A2(5)
c       A(2,3)=A2(6)
c       A(3,1)=A2(7)
c       A(3,2)=A2(8)
c       A(3,3)=A2(9)
C
c      IF(QPRINT) write(6,56) ' The A matrix:',A
c      IF(QPRINT) write(6,56) ' The EV matrix:',EV
C
      DO I=1,3
        EVA(I)=SQRT(ABS(EV(I)))
      ENDDO
      IF(DET.LT.0)then
         EVA(1)=-EVA(1)
      endif   
C
c      CALL TRANSP(RT,R,3,3)
c      CALL MULNXN(RTA,RT,A,3)

       do i=1,3
          do j=1,3
             do k=1,3
                RTA(i,j)=RTA(i,j)+R(i,k)*A(k,j)
             enddo
          enddo
       enddo

C
      DO J=1,3
        EVS=0
        IF(ABS(EVA(J)).GT.0.00001) EVS=1/EVA(J)
        DO I=1,3
          B(I,J)=RTA(I,J)*EVS
        ENDDO
      ENDDO

c      print*,DET,eva(1),eva(2),eva(3),A(1,1),A(2,1),A(3,1)
c      stop
C
c      IF(QPRINT) write(6,56) ' The B matrix:',B
C
      CALL DETM33(A, DET)
c      IF(QPRINT) write(6,56) ' The Adet value:',DET
C
C handle the colinear case (two zero eigenvalues)
C
      B1(1)=B(1,1)
      B1(2)=B(2,1)
      B1(3)=B(3,1)
c
      B2(1)=B(1,2)
      B2(2)=B(2,2)
      B2(3)=B(3,2)

      B3(1)=B(1,3)
      B3(2)=B(2,3)
      B3(3)=B(3,3)

 

      
      CALL NORMALL(B3,3)
      IF(ABS(EVA(2)).LE.0.00001) THEN
          B2(1)=0
          B2(2)=0
          B2(3)=1
          CALL ORTHOG(B2,B3,3)
          CALL DOTPR(B2,B2,3,DOT)
          IF(DOT.LT.0.00001) THEN
            B2(1)=0
            B2(2)=1
            B2(3)=0
            CALL ORTHOG(B2,B3,3)
          ENDIF
      ENDIF
      CALL NORMALL(B2,3)

c handle the coplanar case (one zero eigenvalue)

      IF(ABS(EVA(1)).LE.0.00001) THEN

       B(1,2)=B2(1)
       B(2,2)=B2(2)
       B(3,2)=B2(3)
       B(1,3)=B3(1)
       B(2,3)=B3(2)
       B(3,3)=B3(3)

        DO K=1,3
          KP=K+1
          KQ=K+2
          IF(KP.GT.3) KP=KP-3
          IF(KQ.GT.3) KQ=KQ-3
          B(K,1)=B(KP,2)*B(KQ,3)-B(KP,3)*B(KQ,2)
c          A(K,1)=A(KP,2)*A(KQ,3)-A(KP,3)*A(KQ,2)
          IF(DET.LT.0) B(K,1)=-B(K,1)
        ENDDO

       B1(1)=B(1,1)
       B1(2)=B(2,1)
       B1(3)=B(3,1)

      ENDIF

      CALL NORMALL(B1,3)

c      IF(QPRINT) write(6,56) ' The B matrix:',B

       B(1,1)=B1(1)
       B(2,1)=B1(2)
       B(3,1)=B1(3)
       B(1,2)=B2(1)
       B(2,2)=B2(2)
       B(3,2)=B2(3)
       B(1,3)=B3(1)
       B(2,3)=B3(2)
       B(3,3)=B3(3)


c       do i=1,3
c      print*,B1(i),B2(i),B3(i)
c      enddo


c      CALL TRANSP(BT,B,3,3)
c      CALL MULNXN(U2,A,BT,3)
 
      do i=1,3
         do j=1,3
            do k=1,3
               U2(i,j)=U2(i,j)+B(i,k)*A(j,k)
            enddo
         enddo
      enddo

  
      U(1)=U2(1,1)
      U(2)=U2(1,2)
      U(3)=U2(1,3)
      U(4)=U2(2,1)
      U(5)=U2(2,2)
      U(6)=U2(2,3)
      U(7)=U2(3,1)
      U(8)=U2(3,2)
      U(9)=U2(3,3)

c      do i=1,9
c         print*,u(i)
c         enddo
c      stop
C
c      IF(QPRINT) write(6,56) ' The U matrix:',U
C
C Check to insure unity (as opposed to anti-unitary)
      CALL DETM33(U, DET)

c      print*,det
c      IF(ABS(DET-1).GT.0.00001 .AND. WRNLEV.GT.0) WRITE(OUTU,55) DET
c 55   FORMAT(/' ***** WARNING ***** FROM FROTU. ROTATION MATRIX IS ',
c     $       'NOT UNITARY.'/,'  THE DETERMINANT IS',F14.8/)
C
c      CALL COPYR8(U,DEVA,9)
C
C--------------------------------------------------------------------
C Adjust the results if near a "fault"
C compute both solutions and use a switching function to get
C a linear combination.
C
c      IF(EVA(1)+EVA(2).GT.EVWID) RETURN
c      IF(EVA(1)*EVA(2).GT.0.0) RETURN
c      IF(EVA(2)-EVA(1).LT.0.00005) RETURN

c     DRAT=(EVWID-EVA(1)-EVA(2))/EVWID
C quadratic formula
c      SWVAL=HALF*DRAT**2
c      DSWVAL= TWO*(EVA(2)+EVA(1)) * DRAT/EVWID
CL      SWVAL=HALF*DRAT
CL      DSWVAL= (EVA(2)+EVA(1))/EVWID
C
c     IF(QPRINT) write(6,56) ' The SWVAL value:',SWVAL
C
c      DO J=1,2
c        EVS=1.0/EVA(J)
c        DO I=1,3
c          B(I,J)=RTA(I,J)*EVS
c        ENDDO
c      ENDDO
C
c      CALL VZERO(B(1,3),3)
c      CALL NORMALL(B(1,2),3)
c      CALL NORMALL(B(1,1),3)
C
c      IF(QPRINT) write(6,56) ' The BSW matrix:',B
C
c      CALL TRANSP(BT,B,3,3)
c      CALL MULNXN(DSW,A,BT,3)
C
c      IF(QPRINT) write(6,56) ' The DSW matrix:',DSW
C
C Now average the solutions
C
c      DO I=1,3
c        DO J=1,3
c          DEVA(I,J) = DEVA(I,J) + DSW(I,J)*(DSWVAL-TWO*SWVAL)
c        ENDDO
c      ENDDO
c      DO I=1,2
c        EVA(I)=EVA(I)*(ONE-TWO*SWVAL)
c      ENDDO
C--------------------------------------------------------------------
C
      RETURN
      END

c      SUBROUTINE PKMASS(AMASS1,AMASS2,BMASS,ATOMPR,NPAIR,LMASS,
c     $                  QWGHT,KWMULT)
C-----------------------------------------------------------------------
C     THIS ROUTINE PACKS THE BMASS ARRAY WITH MASSES FROM SELECTED
C     ATOMS. 1.0 IS USED IF LMASS IS FALSE
C     IF QWGHT IS TRUE THE PRODUCT OF BMASS AND KWMULT IS FORMED AND
C     STORED IN BMASS
C
C     Author: Bernie Brooks
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/number.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
C
c      INTEGER NPAIR
c      REAL*8 AMASS1(*),AMASS2(*)
c      REAL*8 KWMULT(*)
c      INTEGER ATOMPR(2,*)
c      REAL*8 BMASS(*)
c      LOGICAL LMASS,QWGHT
c      INTEGER NWARN,K,KX,KY
cC
c      IF (LMASS) THEN
c        NWARN=0
c        DO K=1,NPAIR
c          KX=ATOMPR(1,K)
c          KY=ATOMPR(2,K)
c          BMASS(K)=AMASS1(KX)
c          IF(AMASS1(KX).NE.AMASS2(KY)) THEN
c            NWARN=NWARN+1
c            BMASS(K)=SQRT(AMASS1(KX)*AMASS2(KY))
c          ENDIF
c        ENDDO
c        IF(NWARN.GT.0 .AND. WRNLEV.GE.2) WRITE(OUTU,45) NWARN
c      ELSE
c        DO K=1,NPAIR
c          BMASS(K)=ONE
c        ENDDO
c      ENDIF
c      IF (QWGHT) THEN
c        DO K=1,NPAIR
c          BMASS(K)=BMASS(K)*KWMULT(ATOMPR(1,K))
c        ENDDO
c      ENDIF
cC
c 45   FORMAT(/' *** WARNING *** MASSES DONT MATCH FOR THIS',
c     $ ' HOMOLOGY FOR',I5,' ATOMS. RESULTS WILL USE GEOMETRIC MEAN.')
cC
c      RETURN
c      END

      SUBROUTINE FNDROT(UX2,RN,PHI)
C      ,LPRINT)
C-----------------------------------------------------------------------
C     This routine finds the direction and magnitude of a 3x3
C     rotation matrix
C     BERNARD R. BROOKS (overhauled Feb. 1998 - BRB)
C
C   UX(3,3)  - Input 3x3 rotation matrix.  Checked for Det(UX)=1
C   RN(3)    - Returned normalized rotation vector (or zero if UX=1)
C   PHI      - Returned angle of rotation (in degrees)
C   LPRINT   - Logical flag requesting a print of RN and PHI.
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/exfunc.fcm'
c##INCLUDE '~/charmm_fcm/number.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
C
      real*8 ux2(9)
      REAL*8 UX(3,3),RN(3),PHI
c      LOGICAL LPRINT
C
c##INCLUDE '~/charmm_fcm/consta.fcm'
C
      REAL*8 RX(6),RY(6),U(3,3)
      LOGICAL LOK 
c      QPRINT
      INTEGER I,J
      REAL*8 DET,FACT
C
c      QPRINT = LPRINT .AND. (PRNLEV.GT.2)
C
C Test to make sure it is unitary
C
      ux(1,1)=ux2(1) 
      ux(1,2)=ux2(2) 
      ux(1,3)=ux2(3) 
      ux(2,1)=ux2(4) 
      ux(2,2)=ux2(5) 
      ux(2,3)=ux2(6) 
      ux(3,1)=ux2(7) 
      ux(3,2)=ux2(8) 
      ux(3,3)=ux2(9) 


      CALL DETM33(UX, DET)
      IF(ABS(DET-1).GT.0.00001) THEN
c        IF(WRNLEV.GE.2) WRITE(OUTU,15) DET
c 15     FORMAT(' **** WARNING FROM FNDROT **** INPUT MATRIX IS NOT',
c     $       ' UNITARY. DET=',F12.6)
        RETURN
      ENDIF
C
      RN(1)=(UX(3,2)-UX(2,3))*0.5
      RN(2)=(UX(1,3)-UX(3,1))*0.5
      RN(3)=(UX(2,1)-UX(1,2))*0.5
      CALL DOTPR(RN,RN,3,FACT)
      IF(FACT.LT.0.00000001) THEN
C     angle of rotation is either 0.0 or 180.0, use cosine method
        RN(1)=SQRT(ABS((1+UX(1,1)-UX(2,2)-UX(3,3))*0.5))
        RN(2)=SQRT(ABS((1-UX(1,1)+UX(2,2)-UX(3,3))*0.5))
        RN(3)=SQRT(ABS((1-UX(1,1)-UX(2,2)+UX(3,3))*0.5))
        IF(UX(2,1)+UX(1,2).LT.0) RN(2)=-RN(2)
        IF(UX(1,3)+UX(3,1).LT.0) RN(3)=-RN(3)
        IF(UX(3,2)+UX(2,3).LT.0 .AND.
     &           RN(2)*RN(3).GT.ZERO) RN(3)=-RN(3)
        CALL DOTPR(RN,RN,3,FACT)
        IF(FACT.LT.TENM5) THEN
c          IF(QPRINT) WRITE(6,'(A)')
c     $             ' VERY LITTLE ROTATION: NO AXIS FOUND'
          CALL VZERO(RN,3)
          PHI=ZERO
          RETURN
        ELSE
          CALL NORMALL(RN,3)
          PHI=180.0
        ENDIF
      ELSE
C       use sine method
        CALL NORMALL(RN,3)
        FACT=SQRT(FACT)
        IF(FACT.GT.1) FACT=1
        FACT=ASIN(FACT)
        PHI=RADDEG*FACT
C
C Determine if angle>90 degrees
        RY(1)=RN(2)*RN(3)
        RY(2)=RN(1)*RN(3)
        RY(3)=RN(1)*RN(2)
        RY(4)=RN(2)**2+RN(3)**2
        RY(5)=RN(1)**2+RN(3)**2
        RY(6)=RN(1)**2+RN(2)**2
        RX(1)=(UX(3,2)+UX(2,3))*0.5
        RX(2)=(UX(1,3)+UX(3,1))*0.5
        RX(3)=(UX(2,1)+UX(1,2))*0.5
        RX(4)=1-UX(1,1)
        RX(5)=1-UX(2,2)
        RX(6)=1-UX(3,3)
        DO I=1,6
          RX(I)=RX(I)-RY(I)
        ENDDO
        CALL DOTPR(RX,RY,6,FACT)
        IF(FACT.GT.0) PHI=180.0-PHI
      ENDIF
C
c      IF(QPRINT) WRITE(6,35) RN,PHI
c 35   FORMAT(' AXIS OF ROTATION IS',3F10.6,'  ANGLE IS',F8.2/)
C
      CALL FNDU(U,RN,PHI,LOK)
      DO I=1,3
        DO J=1,3
          LOK=LOK.AND.(ABS(U(I,J)-UX(I,J)).LT.0.00005)
        ENDDO
      ENDDO
c      IF (.NOT.(LOK)) THEN
c        IF(WRNLEV.GE.2) WRITE(OUTU,45)
c        IF(WRNLEV.GE.2) WRITE(OUTU,55) U
c      ENDIF
c 45   FORMAT(/' ** WARNING ** IN FNDROT. ROTATION AXIS POORLY',
c     $       ' RESOLVED')
c 55   FORMAT(' RESULTANT ROTATION MATRIX',3(/1X,3F12.6))
C
      RETURN
      END

      SUBROUTINE FNDU(U,RN,PHI,LOK)
C-----------------------------------------------------------------------
C     This routine finds the unitary transformation matrix
C     given an axis and an angle about that axis.  It also
C     normalizes the rotation axis vector.
C
C     BERNARD R. BROOKS   (overhauled Feb. 1998 - BRB)
C
C   U(3,3)   - Returnted 3x3 rotation matrix.
C   RN(3)    - Input rotation vector (!!It will be normalized!!)
C   PHI      - Input angle of rotation (in degrees)
C   LOK      - Returned logical flag indicating proper function
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/exfunc.fcm'
c##INCLUDE '~/charmm_fcm/number.fcm'
C
      REAL*8 U(3,3),RN(3),PHI
      LOGICAL LOK
C
      REAL*8 PHIR,FACT,T(3,3),T2(3,3),VAL1,VAL2
      INTEGER I,J
c##INCLUDE '~/charmm_fcm/consta.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
C
      PHIR=PHI*3.141592653589793/180
      do J=1,3
        FACT=RN(J)*RN(J)
      enddo
      IF(ABS(FACT).LT.TENM5) THEN
c        IF(WRNLEV.GE.2) WRITE(OUTU,35) FACT
c 35     FORMAT(/' **** WARNING ***** THE ROTATION AXIS MUST BE',
c     $       ' SPECIFIED.'/,' NO ROTATION CAN TAKE PLACE',D12.6)
        LOK=.FALSE.
        RETURN
      ENDIF
C
      CALL NORMALL(RN,3)
C
      DO I=1,3
        DO J=1,3
          U(I,J)=ZERO
        ENDDO
        U(I,I)=ONE
        T(I,I)=ZERO
      ENDDO
C
      T(3,2)= RN(1)
      T(2,3)=-RN(1)
      T(1,3)= RN(2)
      T(3,1)=-RN(2)
      T(2,1)= RN(3)
      T(1,2)=-RN(3)
C
      VAL1=SIN(PHIR)
      VAL2=1-COS(PHIR)
C
      CALL ADDCTV(U,T,9,VAL1)
      CALL MULNXN(T2,T,T,3)
      CALL ADDCTV(U,T2,9,VAL2)
C
      LOK=.TRUE.
C
      RETURN
      END

      SUBROUTINE DETM33(U, DET)
C
C     DET is the determinant of the double precision 3x3 matrix U.
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
C
      REAL*8  DET, U(3,3)
C
      DET = U(1,1) * (U(2,2)*U(3,3) - U(2,3)*U(3,2)) -
     *      U(1,2) * (U(2,1)*U(3,3) - U(3,1)*U(2,3)) +
     *      U(1,3) * (U(2,1)*U(3,2) - U(3,1)*U(2,2))
C
      RETURN
c##IF UNUSED (diagrs_unused)
      END
 
      SUBROUTINE MULNXN(RESULT, A, B, N)
C
C     Multiply the NxN matrices A and B to give RESULT. The matrices
C     are assumed to be stored column first and are double precision.
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
C
      INTEGER I, J, K, N
      REAL*8  A(N,N), B(N,N), RESULT(N,N)
C
      REAL*8  ZERO
      PARAMETER (ZERO = 0.D0)
C
      DO 30 I = 1,N
         DO 20 J = 1,N
            RESULT(J,I) = ZERO
            DO 10 K = 1,N
               RESULT(J,I) = RESULT(J,I) + A(J,K) * B(K,I)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RETURN
      END     
      
      SUBROUTINE TRANSP(RESULT, A, N1, N2)
C
C     Transpose a double precision N1xN2 matrix, A, into RESULT.
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
C
      INTEGER I, J, N1, N2
      REAL*8  A(N1,N2), RESULT(N2,N1)
C
      DO 20 I = 1,N1
         DO 10 J = 1,N2
            RESULT(J,I) = A(I,J)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END

C***********************************************************************
      SUBROUTINE VZERO(V,N)
C-----------------------------------------------------------------------
C     ZEROES OT A VECTOR OF LENGTH N
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
      INTEGER N,I
      REAL*8 V(N)
C
      DO 10 I=1,N
        V(I)=0.0
 10   CONTINUE
      RETURN
      END

      SUBROUTINE DOTPR(V,W,N,CX)
C-----------------------------------------------------------------------
C     COMPUTES THE INNER PRODUCT OF TWO VECTORS
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
      INTEGER N,I
      REAL*8 V(N),W(N),C,CX
C
      C=0.0
      DO 10 I=1,N
        C=C+V(I)*W(I)
 10   CONTINUE
      CX=C
      RETURN
      END

      SUBROUTINE ADDCTV(V,W,N,CX)
C-----------------------------------------------------------------------
C     ADDS A SCALAR TIMES VECTOR W TO VECTOR V.
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
      INTEGER N,I
      REAL*8 V(N),W(N),C,CX
C
      C=CX
      DO 10 I=1,N
        V(I)=V(I)+W(I)*C
 10   CONTINUE
      RETURN
      END  
     
      SUBROUTINE ORTHOG(V,W,N)
C-----------------------------------------------------------------------
C     ORTHOGONALIZES VECTOR V TO W WHEN W IS NORMALIZED
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
      INTEGER N
      REAL*8 V(N),W(N),C
C
      CALL DOTPR(V,W,N,C)
      C=-C
      CALL ADDCTV(V,W,N,C)
      RETURN
      END

      SUBROUTINE NORMALL(V,N)
C-----------------------------------------------------------------------
C     NORMALIZES VECTOR V OF LENGTH N
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
      INTEGER N,I
      REAL*8 V(N),C
C
      C=0.0
      DO 10 I=1,N
        C=C+V(I)*V(I)
 10   CONTINUE
      IF(C.LT.0.00000000001) THEN
c        IF(WRNLEV.GE.2) WRITE(OUTU,25) C
cCRCZ    CALL DIE
        DO 20 I=1,N
          V(I)=0.0
 20     CONTINUE
        RETURN
      ENDIF
c 25   FORMAT(' **** WARNING **** TRYING TO NORMALIZE A ZERO VECTOR',
c     $       ' NORM=',E12.4/' IT WILL BE ZEROED.')
C
      C=1.0/SQRT(C)
      DO 30 I=1,N
        V(I)=V(I)*C
 30   CONTINUE
C
      RETURN
      END
   
      SUBROUTINE COPYR8(A,COPY,N)
C-----------------------------------------------------------------------
C     Copies array A into COPY, applies to REAL*8 arrays
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
      REAL*8 A(*), COPY(*)
      INTEGER N, I
      IF(N.GT.0) THEN
        DO I=1,N
          COPY(I)=A(I)
        ENDDO
      ENDIF
      RETURN
      END

c      SUBROUTINE FRESTK(NUMWRD)
C-----------------------------------------------------------------------
C     Frees NUMWRD words off the stack
C
C      Author: Robert Bruccoleri
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c##INCLUDE '~/charmm_fcm/stack.fcm'
c##INCLUDE '~/charmm_fcm/stream.fcm'
c##INCLUDE '~/charmm_fcm/timer.fcm'
c      INTEGER NUMWRD,LEN
C
c      LEN=NUMWRD
c      IF (LEN.LT.0) THEN
c        CALL WRNDIE(-4,'<FRESTK>','Number of words less than zero')
c      ENDIF
c      IF(LEN.NE.(LEN/2)*2) LEN=LEN+1
c      LSTUSD=LSTUSD-LEN
c      IF(LSTUSD.LT.0) CALL WRNDIE(-4,'<FRESTK>','Stack underflow')
cC
c      IF(PRSTCK .AND. PRNLEV.GE.2) WRITE(OUTU,45) LSTUSD,LEN
c  45  FORMAT(' STACK> FRESTK CALLED:: LSTUSD=',I7,'  NUMWRD=',I7)
C
c      RETURN
c      END

c      SUBROUTINE SETMSR(WRD,VAL)
C-----------------------------------------------------------------------
C     THIS ROUTINE ADDS A NEW COMMAND PARAMETER TO THE
C     SET OF MISCELLANEOUS VALUES
C
C      Bernard R. Brooks   2/88, L.Nilsson MAY-90
C      now adds name and value if name not already present
C
c##INCLUDE '~/charmm_fcm/impnon.fcm'
c      CHARACTER*(*) WRD
c      REAL*8 VAL
c      INTEGER I
C
c##INCLUDE '~/charmm_fcm/mscpar.fcm'
C
c      DO I=1,MSCNUM
c        IF(WRD.EQ.MSCNAM(I)) THEN
c          MSCVAL(I)=VAL
c          RETURN
c        ENDIF
c      ENDDO
C
C     Couldn't find it. add a new one
C
c      IF(MSCNUM.EQ.MAXMSC) THEN
c        CALL WRNDIE(0,'<SETMSR>','TOO MANY SUBSTITUTION PARAMETERS')
c        RETURN
c      ENDIF
C
c      MSCNUM=MSCNUM+1
c      MSCNAM(MSCNUM)=WRD
c      MSCVAL(MSCNUM)=VAL
c      RETURN
c      END
c
c
c**********************************************************************
cc%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
