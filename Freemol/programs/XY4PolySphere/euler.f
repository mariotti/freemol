C-----------------------------------------------------------------------------
      SUBROUTINE EULER(NATOM,N,U,Y,XE,THETA,PHI,PSI,NTRIAL,IRC)
Cthis subroutine calculates the Euler angles for the rotation of any
Cnuclear frame given by the nuclear positions Y of the NATOM nuclei
Cto the Eckart reference frame (NATOM = N/3) and performs the rotation.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        PARAMETER(TOL=1.0D-06,MTRIAL=800)
        DIMENSION U(NATOM),Y(N),XE(N),XX(3),YY(3)
        DIMENSION CP(3,3),ANGLES(3),F(3),DF(3,3)
        ZABS(X)=ABS(X)
        IRC=0
Ctest Schwerpunktsystem und Haupttraegheitsachsen
        NTEST=N/3
        IF (NATOM.NE.NTEST) THEN
         WRITE(*,*) ' SUBROUTINE EULER; MISMATCH WITH NATOM'
         IRC=-1
         RETURN
        ENDIF
        SUM1=0.0
        SUM2=0.0
        SUM3=0.0
        DO 1 J=1,NATOM
         SUM1=SUM1+U(J)*Y(3*J-2)
         SUM2=SUM2+U(J)*Y(3*J-1)
    1    SUM3=SUM3+U(J)*Y(3*J)
        SUM=ZABS(SUM1)+ZABS(SUM2)+ZABS(SUM3)
        IF (SUM.GT.TOL) THEN
         WRITE(*,*) ' SUBROUTINE EULER; SUM(MJ*YJ) NOT ZERO',SUM
         IRC=-1
         RETURN
        ENDIF
        SUM1=0.0
        SUM2=0.0
        SUM3=0.0
        DO 2 J=1,NATOM
         SUM1=SUM1+U(J)*XE(3*J-2)
         SUM2=SUM2+U(J)*XE(3*J-1)
    2    SUM3=SUM3+U(J)*XE(3*J)
        SUM=ZABS(SUM1)+ZABS(SUM2)+ZABS(SUM3)
        IF (SUM.GT.TOL) THEN
         WRITE(*,*) ' SUBROUTINE EULER; SUM(MJ*XEJ) NOT ZERO'
         IRC=-1
         RETURN
        ENDIF
        SUM1=0.0
        SUM2=0.0
        SUM3=0.0
        DO 3 J=1,NATOM
         SUM1=SUM1+U(J)*XE(3*J-2)*XE(3*J-1)
         SUM2=SUM2+U(J)*XE(3*J-1)*XE(3*J)
    3    SUM3=SUM3+U(J)*XE(3*J)*XE(3*J-2)
        SUM=ZABS(SUM1)+ZABS(SUM2)+ZABS(SUM3)
        IF (SUM.GT.TOL) THEN
         WRITE(*,*) ' SUBROUTINE EULER; PRODUCTS OF INERTIA NOT ZERO;'
         IRC=-1
         RETURN
        ENDIF
Ccalculate cross products of XE with Y
        DO 4 K=1,3
        DO 4 I=1,3
         SUM=0.0
         DO 5 J=0,(NATOM-1)
    5     SUM=SUM+U(J+1)*XE(3*J+I)*Y(3*J+K)
    4    CP(I,K)=SUM
c        print '(3(x,g16.9))', (cp(1,k),k=1,3)
c        print '(3(x,g16.9))', (cp(2,k),k=1,3)
c        print '(3(x,g16.9))', (cp(3,k),k=1,3)
Cinital guesses for THETA, PHI and PSI
C     READ(*,*) (ANGLES(I),I=1,3)
C      READ(*,*) NTRIAL
        ANGLES(1) = THETA
        ANGLES(2) = PHI
        ANGLES(3) = PSI
      WRITE(*,*) ' '
      WRITE(*,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,'(A,/)')' PROTOCOL FOR SOLUTION OF ECKART EQUATIONS'
      WRITE(*,*)' INITIAL GUESSES'
      WRITE(*,*) (ANGLES(I),I=1,3)
Csolve for roots of Eckart equations with Newton-Raphson
C       NTRIAL=MTRIAL
        NROOT=3
        TOLX=TOL
        TOLF=TOL
        IF(NTRIAL.GT.0) CALL MNEWT(NTRIAL,ANGLES,NROOT,TOLX,TOLF,CP)
Ccheck results (test 1)
        THETA = ANGLES(1)
        PHI   = ANGLES(2)
        PSI   = ANGLES(3)
        IF (ABS(THETA).GT.360.) THEN
        IA = AINT(THETA/360.)
        THETA = THETA - FLOAT(IA)*360.
        ENDIF
        IF(THETA.LT.0.) THEN
         THETA=-THETA
         PHI=180.+PHI
         PSI=180.+PSI
        ENDIF
        IF(THETA.GT.180.) THEN
         THETA=360.-THETA
         PHI=180.+PHI
         PSI=180.+PSI
        ENDIF
        IF (ABS(PHI).GT.360.) THEN
        IA = AINT(PHI/360.)
        PHI = PHI - FLOAT(IA)*360.
        ENDIF
        IF(PHI.LT.0.) PHI=PHI+360.
        IF (ABS(PSI).GT.360.) THEN
        IA = AINT(PSI/360.)
        PSI = PSI - FLOAT(IA)*360.
        ENDIF
        IF(PSI.LT.0.) PSI=PSI+360.
        ANGLES(1)=THETA
        ANGLES(2)=PHI
        ANGLES(3)=PSI
        WRITE(*,*)' EULER ANGLES :'
        WRITE(*,*) THETA, PHI, PSI
        CALL ECKART(ANGLES,DF,F,CP)
        WRITE(*,*)' DIRECT TEST OF ECKART CONDITIONS :'
        WRITE(*,*) (F(I),I=1,3)
Ctest 2: rotate coordinates and check eckart conditions
 11   CONTINUE
      DO 300 J=0,(NATOM-1)
       XX(1)=Y(J*3+1)
       XX(2)=Y(J*3+2)
       XX(3)=Y(J*3+3)
       CALL REF(THETA,PHI,PSI,XX,YY)
       Y(J*3+1)=YY(1)
       Y(J*3+2)=YY(2)
       Y(J*3+3)=YY(3)
 300  CONTINUE
Ccalculate cross products of XE with rotated Y (test 2)
        DO 304 K=1,3
        DO 304 I=1,3
         SUM=0.0
         DO 305 J=0,(NATOM-1)
  305     SUM=SUM+U(J+1)*XE(3*J+I)*Y(3*J+K)
  304    CP(I,K)=SUM
        ANGLES(1) = 0.
        ANGLES(2) = 0.
        ANGLES(3) = 0.
        CALL ECKART(ANGLES,DF,F,CP)
      WRITE(*,*)' TEST OF ECKART CONDITIONS WITH ROTATED COORDINATES:'
      WRITE(*,*) (F(I),I=1,3)
      WRITE(*,*) ' ++++++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) ' '
      RETURN
      END
        SUBROUTINE ECKART(X,DF,F,PAR)
Cthis subroutine calculates the Eckart conditions suitable for use
Cin the Newton-Raphson routine.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        DIMENSION CHI(3,3),X(3),F(3),DF(3,3),PAR(3,3)
        DIMENSION CHITHE(3,3),CHIPHI(3,3),CHIPSI(3,3)
        ZCOS(R)=COS(R)
        ZSIN(R)=SIN(R)
        ZACOS(R)=ACOS(R)
        PI=ZACOS(-1.D0)
        URF=PI/180.
        THETA=X(1)
        PHI=X(2)
        PSI=X(3)
        APSI=PSI*URF
        APHI=PHI*URF
        ATHE=THETA*URF
        CPSI=ZCOS(APSI)
        CPHI=ZCOS(APHI)
        CTHE=ZCOS(ATHE)
        SPSI=ZSIN(APSI)
        SPHI=ZSIN(APHI)
        STHE=ZSIN(ATHE)
        CHI(1,1)= CTHE*CPHI*CPSI-SPHI*SPSI
        CHI(1,2)= CTHE*SPHI*CPSI+CPHI*SPSI
        CHI(1,3)=-STHE*CPSI
        CHI(2,1)=-CTHE*CPHI*SPSI-SPHI*CPSI
        CHI(2,2)=-CTHE*SPHI*SPSI+CPHI*CPSI
        CHI(2,3)= STHE*SPSI
        CHI(3,1)= STHE*CPHI
        CHI(3,2)= STHE*SPHI
        CHI(3,3)= CTHE
c        print '(3(x,g16.9))', (chi(1,k),k=1,3)
c        print '(3(x,g16.9))', (chi(2,k),k=1,3)
c        print '(3(x,g16.9))', (chi(3,k),k=1,3)
        CHITHE(1,1)=-STHE*CPHI*CPSI
        CHITHE(1,2)=-STHE*SPHI*CPSI
        CHITHE(1,3)=-CTHE*CPSI
        CHITHE(2,1)= STHE*CPHI*SPSI
        CHITHE(2,2)= STHE*SPHI*SPSI
        CHITHE(2,3)= CTHE*SPSI
        CHITHE(3,1)= CTHE*CPHI
        CHITHE(3,2)= CTHE*SPHI
        CHITHE(3,3)=-STHE
        CHIPHI(1,1)=-CHI(1,2)
        CHIPHI(1,2)= CHI(1,1)
        CHIPHI(1,3)= 0.0
        CHIPHI(2,1)=-CHI(2,2)
        CHIPHI(2,2)= CHI(2,1)
        CHIPHI(2,3)= 0.0
        CHIPHI(3,1)=-CHI(3,2)
        CHIPHI(3,2)= CHI(3,1)
        CHIPHI(3,3)= 0.0
        CHIPSI(1,1)= CHI(2,1)
        CHIPSI(1,2)= CHI(2,2)
        CHIPSI(1,3)= CHI(2,3)
        CHIPSI(2,1)=-CHI(1,1)
        CHIPSI(2,2)=-CHI(1,2)
        CHIPSI(2,3)=-CHI(1,3)
        CHIPSI(3,1)= 0.0
        CHIPSI(3,2)= 0.0
        CHIPSI(3,3)= 0.0
        DO 10 J=1,3
  10    F(J)=0.0
        DO 11 J=1,3
        DO 11 L=1,3
  11    DF(L,J)=0.0
        DO 1 J=1,3
         JN=J+1
         IF(J.EQ.3) JN=1
         DO 1111 I=1,3
          F(J)=F(J)+PAR(J,I)*CHI(JN,I)
c         print *, f(j)
          F(J)=F(J)-PAR(JN,I)*CHI(J,I)
c         print *, f(j)
 1111    continue
         F(J)=-F(J)
    1   continue
c    1     F(J)=-F(J)
        DO 2 J=1,3
         JN=J+1
         IF(J.EQ.3) JN=1
         DO 2 I=1,3
          DF(J,1)=DF(J,1)+PAR(J,I)*CHITHE(JN,I)*URF
    2     DF(J,1)=DF(J,1)-PAR(JN,I)*CHITHE(J,I)*URF
        DO 3 J=1,3
         JN=J+1
         IF(J.EQ.3) JN=1
         DO 3 I=1,3
          DF(J,2)=DF(J,2)+PAR(J,I)*CHIPHI(JN,I)*URF
    3     DF(J,2)=DF(J,2)-PAR(JN,I)*CHIPHI(J,I)*URF
        DO 4 J=1,3
         JN=J+1
         IF(J.EQ.3) JN=1
         DO 4 I=1,3
          DF(J,3)=DF(J,3)+PAR(J,I)*CHIPSI(JN,I)*URF
    4     DF(J,3)=DF(J,3)-PAR(JN,I)*CHIPSI(J,I)*URF
        RETURN
        END
      SUBROUTINE MNEWT(NTRIAL,X,N,TOLX,TOLF,PAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (NP=3)
      DIMENSION X(NP),ALPHA(NP,NP),BETA(NP),INDX(NP),PAR(3,3)
      DO 13  K=1,NTRIAL
        CALL ECKART(X,ALPHA,BETA,PAR)
        ERRF=0.
        DO 11 I=1,N
          ERRF=ERRF+ABS(BETA(I))
11      CONTINUE
c       PRINT *, ERRF
        IF(ERRF.LE.TOLF)RETURN
        CALL LUDCMP(ALPHA,N,NP,INDX,D)
        CALL LUBKSB(ALPHA,N,NP,INDX,BETA)
        ERRX=0.
        DO 12 I=1,N
          ERRX=ERRX+ABS(BETA(I))
          X(I)=X(I)+BETA(I)
12      CONTINUE
        IF(ERRX.LE.TOLX)RETURN
13    CONTINUE
      RETURN
      END
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (NMAX=100,TINY=1.0E-20)
      DIMENSION A(NP,NP),INDX(N),VV(NMAX)
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'SINGULAR MATRIX.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      DIMENSION A(NP,NP),INDX(N),B(N)
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
        SUBROUTINE REF(THETA,PHI,PSI,Y,X)
Cthis subroutine rotates the vector Y given in a molecule fixed frame
Cto  a vector X given in the Eckart frame.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        DIMENSION CHI(3,3),X(3),Y(3)
        ZCOS(R)=COS(R)
        ZSIN(R)=SIN(R)
        ZACOS(R)=ACOS(R)
        PI=ZACOS(-1.D0)
        URF=PI/180.
        APSI=PSI*URF
        APHI=PHI*URF
        ATHE=THETA*URF
        CPSI=ZCOS(APSI)
        CPHI=ZCOS(APHI)
        CTHE=ZCOS(ATHE)
        SPSI=ZSIN(APSI)
        SPHI=ZSIN(APHI)
        STHE=ZSIN(ATHE)
        CHI(1,1)= CTHE*CPHI*CPSI-SPHI*SPSI
        CHI(1,2)= CTHE*SPHI*CPSI+CPHI*SPSI
        CHI(1,3)=-STHE*CPSI
        CHI(2,1)=-CTHE*CPHI*SPSI-SPHI*CPSI
        CHI(2,2)=-CTHE*SPHI*SPSI+CPHI*CPSI
        CHI(2,3)= STHE*SPSI
        CHI(3,1)= STHE*CPHI
        CHI(3,2)= STHE*SPHI
        CHI(3,3)= CTHE
        DO 1 J=1,3
         SUM=0.0
         DO 2 I=1,3
   2      SUM=SUM+CHI(J,I)*Y(I)
   1     X(J)=SUM
        RETURN
        END
        SUBROUTINE REFT(THETA,PHI,PSI,Y,X)
Cthis subroutine rotates the vector Y given in the Eckart frame into
Ca vector X given in the moleule fixed frame
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
        DIMENSION CHI(3,3),X(3),Y(3)
        ZCOS(R)=COS(R)
        ZSIN(R)=SIN(R)
        ZACOS(R)=ACOS(R)
        PI=ZACOS(-1.D0)
        URF=PI/180.
        APSI=PSI*URF
        APHI=PHI*URF
        ATHE=THETA*URF
        CPSI=ZCOS(APSI)
        CPHI=ZCOS(APHI)
        CTHE=ZCOS(ATHE)
        SPSI=ZSIN(APSI)
        SPHI=ZSIN(APHI)
        STHE=ZSIN(ATHE)
        CHI(1,1)= CTHE*CPHI*CPSI-SPHI*SPSI
        CHI(2,1)= CTHE*SPHI*CPSI+CPHI*SPSI
        CHI(3,1)=-STHE*CPSI
        CHI(1,2)=-CTHE*CPHI*SPSI-SPHI*CPSI
        CHI(2,2)=-CTHE*SPHI*SPSI+CPHI*CPSI
        CHI(3,2)= STHE*SPSI
        CHI(1,3)= STHE*CPHI
        CHI(2,3)= STHE*SPHI
        CHI(3,3)= CTHE
        DO 1 J=1,3
         SUM=0.0
         DO 2 I=1,3
   2      SUM=SUM+CHI(J,I)*Y(I)
   1     X(J)=SUM
        RETURN
        END
c-------------------------------------------------------------------------------
      SUBROUTINE CGTRAF(NATOM,N,U,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(N),U(NATOM)
      XC=0.
      YC=0.
      ZC=0.
      UTOT=0.
      DO NA=1,NATOM
         XC=XC+U(NA)*X(3*NA-2)
         YC=YC+U(NA)*X(3*NA-1)
         ZC=ZC+U(NA)*X(3*NA-0)
         UTOT=UTOT+U(NA)
      ENDDO
      DO NA=1,NATOM
         X(3*NA-2)=X(3*NA-2)-XC/UTOT
         X(3*NA-1)=X(3*NA-1)-YC/UTOT
         X(3*NA-0)=X(3*NA-0)-ZC/UTOT
      ENDDO
      RETURN
      END
