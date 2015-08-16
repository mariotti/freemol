      subroutine jakobi(N,M,A,U,V,NROT) 
!                                                                       
! EIGENVALUES AND EIGENVECTORS OF A real SYMMETRIC MATRIX A OF logical  
! SIZE A(N,N) STORED IN AN ARRAY OF PHYSICAL SIZE A(M,M), where         
! M .ge. N.                                                             
!                                                                       
! ON OUTPUT:                                                            
! (1) ELEMENTS OF A ABOVE THE DIAGONAL ARE DESTROYED.                   
! (2) THE ARRAY U(M) RETURNS THE EIGENVALUES OF A(N,N) IN ITS FIRST N   
!     ELEMENTS.                                                         
! (3) THE COLUMNS OF THE MATRIX V(N,N), STORED IN THE ARRAY V(M,M),     
!     CONTAIN THE EIGENVECTORS OF A.                                    
! (4) THE VARIABLE NROT RETURNS THE NUMBER OF ITERATIONS OF JACOBI      
!     ROTATION THAT WERE REQUIRED TO ANNIHILATE THE OFF-DIAGONAL        
!     ELEMENTS OF A(N,N) TO MACHINE precision.                          
!                                                                       
! FORTRAN CODE ADAPTED FROM WILLIAM H. PRESS, BRIAN P. FLANNERY, SAUL A.
! TEUKOLSKY, AND WILLIAM T. VETTERLING (1986).  NUMERICAL RECIPIES:  THE
!                                               --------- --------   ---
! ART OF SCIENTIFIC COMPUTING, PP. 335-349.  CAMBRIDGE, ENGLAND:        
! --- -- ---------- ---------                                           
! CAMBRIDGE UNIVERSITY PRESS.                                           
!                                                                       
      parameter (NMAX=999) 
      dimension A(M,M),U(M),V(M,M),B(NMAX),Z(NMAX) 
      double precision A,U,V,B,Z,T,THRESH,E,AII,AJJ,THETA,C,S,TAU,P,Q 
!                                                                       
! INITIALIZE EIGENVECTORS MATRIX TO AN IDENTITY MATRIX.                 
!                                                                       
      do I=1,N 
      do J=1,N 
        V(I,J)=0 
      end do 
        V(I,I)=1 
      end do 
!                                                                       
! INITIALIZE EIGENVALUES U AND WORK VECTOR B TO A(I,I) AND ZERO WORK    
! VECTOR Z.                                                             
!                                                                       
      do I=1,N 
        T=A(I,I) 
        U(I)=T 
        B(I)=T 
        Z(I)=0 
      end do 
!                                                                       
! PERFORM UP TO 50 ITERATIONS OF UP TO N*(N - 1)/2 JACOBI ROTATIONS.    
!                                                                       
      NROT=0 
      NCYCLE=1 
      do while (NCYCLE.le.50) 
!                                                                       
! TEST FOR NORMAL return WHEN MAXIMUM MAGNITUDE OF AN OFF-DIAGONAL      
! ELEMENT EQUALS ZERO TO MACHINE precision.  THE TEST PRESUMES THAT     
! ARITHMETIC UNDERFLOW VALUES ARE SET TO ZERO.                          
!                                                                       
        T=0 
        do I=1,N-1 
        do J=I+1,N 
          T=T+abs(A(I,J)) 
        end do 
        end do 
        if (T.eq.0) return 
!                                                                       
! SET OFF-DIAGONAL THRESHOLD.                                           
!                                                                       
        if (NCYCLE.lt.4) then 
          THRESH=T/(5*N**2) 
        else 
          THRESH=0 
        end if 
!                                                                       
! ROTATE TO ANNIHILATE OFF-DIAGONAL ELEMENTS.                           
!                                                                       
        do I=1,N-1 
        do J=I+1,N 
          T=abs(A(I,J)) 
          E=100*T 
          AII=abs(U(I)) 
          AJJ=abs(U(J)) 
          if (NCYCLE.gt.4.and.AII+E.eq.AII.and.AJJ+E.eq.AJJ) then 
!                                                                       
! AFTER FOUR CYCLES, SKIP THE ROTATION if abs(A(I,J)) IS SMALL COMPARED 
! TO BOTH abs(A(I,I)) AND abs(A(J,J)).                                  
!                                                                       
            A(I,J)=0 
          else if (T.gt.THRESH) then 
            T=abs(AJJ-AII) 
            if (T+E.eq.T) then 
!                                                                       
! EFFECTIVELY, SET T = 1/(2*THETA)                                      
!                                                                       
              T=A(I,J)/(AJJ-AII) 
            else 
              THETA=(AJJ-AII)/(2*A(I,J)) 
              T=1/(abs(THETA)+sqrt(1+THETA**2)) 
              if (THETA.lt.0) T=-T 
            end if 
            C=1/sqrt(1+T**2) 
            S=T*C 
            TAU=S/(1+C) 
!                                                                       
! ADJUST EIGENVALUES U AND WORK VECTOR Z.                               
!                                                                       
            E=T*A(I,J) 
            Z(I)=Z(I)-E 
            Z(J)=Z(J)+E 
            U(I)=U(I)-E 
            U(J)=U(J)+E 
            A(I,J)=0 
!                                                                       
! ROTATIONS  1 .le. K .lt. I.                                           
!                                                                       
            do K=1,I-1 
              P=A(K,I) 
              Q=A(K,J) 
              A(K,I)=P-S*(Q+P*TAU) 
              A(K,J)=Q+S*(P-Q*TAU) 
            end do 
!                                                                       
! ROTATIONS  I .lt. K .lt. J.                                           
!                                                                       
            do K=I+1,J-1 
              P=A(I,K) 
              Q=A(K,J) 
              A(I,K)=P-S*(Q+P*TAU) 
              A(K,J)=Q+S*(P-Q*TAU) 
            end do 
!                                                                       
! ROTATIONS  J .lt. K .le. N.                                           
!                                                                       
            do K=J+1,N 
              P=A(I,K) 
              Q=A(J,K) 
              A(I,K)=P-S*(Q+P*TAU) 
              A(J,K)=Q+S*(P-Q*TAU) 
            end do 
!                                                                       
! COMPUTE AND STORE EIGENVECTORS.                                       
!                                                                       
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
!                                                                       
! ADJUST EIGENVALUES AND WORK VECTORS.                                  
!                                                                       
        do I=1,N 
          B(I)=B(I)+Z(I) 
          U(I)=B(I) 
          Z(I)=0 
        end do 
        NCYCLE=NCYCLE+1 
      end do 
      stop '50 CYCLES OF JACOBI ROTATIONS SHOULD NEVER BE NECESSARY.' 
    end subroutine jakobi
