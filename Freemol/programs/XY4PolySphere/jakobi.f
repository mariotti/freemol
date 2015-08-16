      subroutine jakobi (N,M,A,U,V,NROT)
C
C EIGENVALUES AND EIGENVECTORS OF A real SYMMETRIC MATRIX A OF logical
C SIZE A(N,N) STORED IN AN ARRAY OF PHYSICAL SIZE A(M,M), where
C M .ge. N.
C
C ON OUTPUT:
C (1) ELEMENTS OF A ABOVE THE DIAGONAL ARE DESTROYED.
C (2) THE ARRAY U(M) RETURNS THE EIGENVALUES OF A(N,N) IN ITS FIRST N
C     ELEMENTS.
C (3) THE COLUMNS OF THE MATRIX V(N,N), STORED IN THE ARRAY V(M,M),
C     CONTAIN THE EIGENVECTORS OF A.
C (4) THE VARIABLE NROT RETURNS THE NUMBER OF ITERATIONS OF JACOBI
C     ROTATION THAT WERE REQUIRED TO ANNIHILATE THE OFF-DIAGONAL
C     ELEMENTS OF A(N,N) TO MACHINE precision.
C
C FORTRAN CODE ADAPTED FROM WILLIAM H. PRESS, BRIAN P. FLANNERY, SAUL A.
C TEUKOLSKY, AND WILLIAM T. VETTERLING (1986).  NUMERICAL RECIPIES:  THE
C                                               --------- --------   ---
C ART OF SCIENTIFIC COMPUTING, PP. 335-349.  CAMBRIDGE, ENGLAND:
C --- -- ---------- ---------
C CAMBRIDGE UNIVERSITY PRESS.
C
      parameter (NMAX=999)
      dimension A(M,M),U(M),V(M,M),B(NMAX),Z(NMAX)
      double precision A,U,V,B,Z,T,THRESH,E,AII,AJJ,THETA,C,S,TAU,P,Q
C
C INITIALIZE EIGENVECTORS MATRIX TO AN IDENTITY MATRIX.
C
      do I=1,N
      do J=1,N
        V(I,J)=0
      end do
        V(I,I)=1
      end do
C
C INITIALIZE EIGENVALUES U AND WORK VECTOR B TO A(I,I) AND ZERO WORK
C VECTOR Z.
C
      do I=1,N
        T=A(I,I)
        U(I)=T
        B(I)=T
        Z(I)=0
      end do
C
C PERFORM UP TO 50 ITERATIONS OF UP TO N*(N - 1)/2 JACOBI ROTATIONS.
C
      NROT=0
      NCYCLE=1
      do while (NCYCLE.le.50)
C
C TEST FOR NORMAL return WHEN MAXIMUM MAGNITUDE OF AN OFF-DIAGONAL
C ELEMENT EQUALS ZERO TO MACHINE precision.  THE TEST PRESUMES THAT
C ARITHMETIC UNDERFLOW VALUES ARE SET TO ZERO.
C
        T=0
        do I=1,N-1
        do J=I+1,N
          T=T+abs(A(I,J))
        end do
        end do
        if (T.eq.0) return
C
C SET OFF-DIAGONAL THRESHOLD.
C
        if (NCYCLE.lt.4) then
          THRESH=T/(5*N**2)
        else
          THRESH=0
        end if
C
C ROTATE TO ANNIHILATE OFF-DIAGONAL ELEMENTS.
C
        do I=1,N-1
        do J=I+1,N
          T=abs(A(I,J))
          E=100*T
          AII=abs(U(I))
          AJJ=abs(U(J))
          if (NCYCLE.gt.4.and.AII+E.eq.AII.and.AJJ+E.eq.AJJ) then
C
C AFTER FOUR CYCLES, SKIP THE ROTATION if abs(A(I,J)) IS SMALL COMPARED
C TO BOTH abs(A(I,I)) AND abs(A(J,J)).
C
            A(I,J)=0
          else if (T.gt.THRESH) then
            T=abs(AJJ-AII)
            if (T+E.eq.T) then
C
C EFFECTIVELY, SET T = 1/(2*THETA)
C
              T=A(I,J)/(AJJ-AII)
            else
              THETA=(AJJ-AII)/(2*A(I,J))
              T=1/(abs(THETA)+sqrt(1+THETA**2))
              if (THETA.lt.0) T=-T
            end if
            C=1/sqrt(1+T**2)
            S=T*C
            TAU=S/(1+C)
C
C ADJUST EIGENVALUES U AND WORK VECTOR Z.
C
            E=T*A(I,J)
            Z(I)=Z(I)-E
            Z(J)=Z(J)+E
            U(I)=U(I)-E
            U(J)=U(J)+E
            A(I,J)=0
C
C ROTATIONS  1 .le. K .lt. I.
C
            do K=1,I-1
              P=A(K,I)
              Q=A(K,J)
              A(K,I)=P-S*(Q+P*TAU)
              A(K,J)=Q+S*(P-Q*TAU)
            end do
C
C ROTATIONS  I .lt. K .lt. J.
C
            do K=I+1,J-1
              P=A(I,K)
              Q=A(K,J)
              A(I,K)=P-S*(Q+P*TAU)
              A(K,J)=Q+S*(P-Q*TAU)
            end do
C
C ROTATIONS  J .lt. K .le. N.
C
            do K=J+1,N
              P=A(I,K)
              Q=A(J,K)
              A(I,K)=P-S*(Q+P*TAU)
              A(J,K)=Q+S*(P-Q*TAU)
            end do
C
C COMPUTE AND STORE EIGENVECTORS.
C
            do K=1,N
              P=V(K,I)
              Q=V(K,J)
              V(K,I)=P-S*(Q+P*TAU)
              V(K,J)=Q+S*(P-Q*TAU)
            end do
            NROT=NROT+1
          end if
        end do
        end do
C
C ADJUST EIGENVALUES AND WORK VECTORS.
C
        do I=1,N
          B(I)=B(I)+Z(I)
          U(I)=B(I)
          Z(I)=0
        end do
        NCYCLE=NCYCLE+1
      end do
      stop '50 CYCLES OF JACOBI ROTATIONS SHOULD NEVER BE NECESSARY.'
      end
