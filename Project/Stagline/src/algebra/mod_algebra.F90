!------------------------------------------------------------------------------!
!> This module contains subrouroutines and functions for dealing with albegra.
!! Subroutines are taken from Numerical Recipes.
  MODULE mod_algebra

#include"../config.h"

  IMPLICIT NONE

#ifdef GOTO_BLAS
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv_GOTO_BLAS
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: work_GOTO_BLAS
#endif
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dum_vec1, dum_vec2
  REAl(KIND=8), DIMENSION(:), ALLOCATABLE :: Yi_vec
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: alpha_i, invalpha_i, dum_mat
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Ai_mat, Bi_mat, Ci_mat, Gammai_mat
 
  ! Subroutine and functions
  CONTAINS

    !----------------------------------------------------!
    !> This subroutine initializes working arrays when using the 
    !! GOTO BLAS library 
#ifdef GOTO_BLAS
    SUBROUTINE initialize_GOTO_BLAS (r)

      INTEGER, INTENT(IN) :: r

     ALLOCATE(ipiv_GOTO_BLAS(r))
     ALLOCATE(work_GOTO_BLAS(r*r)) 

    END SUBROUTINE initialize_GOTO_BLAS 
#endif
    !----------------------------------------------------! 
    !> This subroutine initialized the arrays needed for the solution of 
    !! a block tridiagonal algebraic system of equations
    SUBROUTINE initialize_Thomas (r, nr)

      INTEGER, INTENT(IN) :: r, nr
 
      ! Array allocation
      ALLOCATE(Ai_mat(r,r*(nr + 1)))
      ALLOCATE(Bi_mat(r,r*(nr + 1)))
      ALLOCATE(Ci_mat(r,nr*r)) 
      ALLOCATE(Gammai_mat(r,(nr - 1)*r)) 
      ALLOCATE(Yi_vec(r*nr))
      ALLOCATE(dum_vec1(r))
      ALLOCATE(dum_vec2(r))
      ALLOCATE(invalpha_i(r,r))    
#ifndef GOTO_BLAS 
      ALLOCATE(alpha_i(r,r))
#endif
      ALLOCATE(dum_mat(r,r))

    END SUBROUTINE initialize_Thomas

    !----------------------------------------------------!
    !> This subroutine solves the tridiagonal block system by means of
    !! the generalized Thomas algotithm. Use is made of LU decomposition
    !! to speed up matrix inversion
    SUBROUTINE tridiag_block_solver (r, nr, d, x)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: r, nr
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: d
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: x
        
      INTEGER :: i, j, k, row, col
      INTEGER :: pos, pos1
      REAL(KIND=8) :: tmp

      ! alphai_i matrix
      DO j = 1,r
         DO i = 1,r
            alpha_i(i,j) = Ai_mat(i,j)
         ENDDO
      ENDDO

      ! Generalized Thomas algorithm for block tridiagonal matrix

      ! Inversion of alpha_i matrix
      CALL inv_matrix (r, alpha_i, invalpha_i)

      ! Matrix-Matrix product: Gammai_mat = invalpha_i*Ci_mat  
      DO col = 1,r
         DO row = 1,r
            Gammai_mat(row,col) = 0.d0
         ENDDO
      ENDDO

      DO col = 1,r
         DO row = 1,r
            tmp = 0.d0
            DO k = 1,r
               tmp = tmp + invalpha_i(row,k)*Ci_mat(k,col)
            ENDDO
            Gammai_mat(row,col) = tmp
         ENDDO
      ENDDO

      ! Matrix-vector product: Yi_vec = invalpha_i*d           
      Yi_vec(1:r) = 0.d0 
      DO col = 1,r
         tmp = d(col)
         DO row = 1,r
            Yi_vec(row) = Yi_vec(row) + invalpha_i(row,col)*tmp 
         ENDDO 
      ENDDO

      ! Loop: forward substitution to find YI
      DO i = 2,nr - 1

         ! Matrix-Matrix product: dum = Bi_mat*Gammai_mat      
         pos  = r*(i - 1)
         pos1 = r*(i - 2) 

         DO col = 1,r
            DO row = 1,r
               dum_mat(row,col) = 0.d0
            ENDDO
         ENDDO
           
         ! Matrix difference: alpha_i = Ai_mat - dum is added in the loop   
         DO col = 1,r
            DO row = 1,r
               tmp = 0.d0
               DO k = 1,r
                  !tmp = tmp + Bi_mat(row,pos1 + k)*Gammai_mat(k,pos1 + col)
                  tmp = tmp + Bi_mat(row,pos + k)*Gammai_mat(k,pos1 + col)
               ENDDO
               dum_mat(row,col)     = tmp
               alpha_i(row,col) = Ai_mat(row,pos + col) - tmp
            ENDDO
         ENDDO
           
         CALL inv_matrix (r, alpha_i, invalpha_i)

         ! Matrix-Matrix product: Gammai_mat = invalpha_i*Ci_mat        
         DO col = 1,r
            DO row = 1,r
               Gammai_mat(row,pos + col) = 0.d0
            ENDDO
         ENDDO

         DO col = 1,r
            DO row = 1,r
               tmp = 0.d0
               DO k = 1,r
                  tmp = tmp + invalpha_i(row,k)*Ci_mat(k,pos + col)
               ENDDO
               Gammai_mat(row,pos + col) = tmp
           ENDDO
         ENDDO

         ! Matrix-vector product: dum_vec1 = invalpha_i*d                      
         dum_vec1 = 0.d0
         DO col = 1,r
            tmp = d(pos + col)
            DO row = 1,r
               dum_vec1(row) = dum_vec1(row) + invalpha_i(row,col)*tmp 
            ENDDO 
         ENDDO 

         ! Matrix-Matrix product: dum = invalpha_i - Bi_mat            
         DO col = 1,r
            DO row = 1,r
               dum_mat(row,col) = 0.d0
            ENDDO
         ENDDO
         
         DO col = 1,r
            DO row = 1,r
               tmp = 0.d0
               DO k = 1,r
                  !tmp = tmp + invalpha_i(row,k)*Bi_mat(k,pos1 + col)
                  tmp = tmp + invalpha_i(row,k)*Bi_mat(k,pos + col)
               ENDDO
               dum_mat(row,col) = tmp
            ENDDO
         ENDDO

         ! Matrix-vector product: dum_vec2 = dum*Yi_vec                    
         dum_vec2 = 0.d0
         pos1 = r*(i - 2)
         DO col = 1,r
            tmp = Yi_vec(pos1 + col)
            DO row = 1,r
               dum_vec2(row) = dum_vec2(row) + dum_mat(row,col)*tmp 
            ENDDO 
         ENDDO 

         DO j = 1,r
            Yi_vec(pos + j) = dum_vec1(j) - dum_vec2(j) 
         ENDDO 
 
      ENDDO
      
      ! Last step of forward loop
      i = nr

      pos  = r*(i - 1)
      pos1 = pos - r  

      ! Matrix-Matrix product:                                       
      DO col = 1,r
         DO row = 1,r
            dum_mat(row,col) = 0.d0
         ENDDO
      ENDDO
           
      ! Matrix difference: alpha_i = Ai_mat - dum is added in the loop  
      DO col = 1,r
         DO row = 1,r
            tmp = 0.d0
            DO k = 1,r
               !tmp = tmp + Bi_mat(row,pos1 + k)*Gammai_mat(k,pos1 + col)
               tmp = tmp + Bi_mat(row,pos + k)*Gammai_mat(k,pos1 + col)
            ENDDO
            dum_mat(row,col)     = tmp
            alpha_i(row,col) = Ai_mat(row,pos + col) - tmp
         ENDDO
      ENDDO

      CALL inv_matrix (R, alpha_i, invalpha_i)                         
 
      ! Matrix-vector product: dum_vec1 = invalpha_i*d                     
      dum_vec1 = 0.d0
      DO col = 1,r
         tmp = d(pos + col)
         DO row = 1,r
            dum_vec1(row) = dum_vec1(row) + invalpha_i(row,col)*tmp 
         ENDDO 
      ENDDO  

      ! Matrix-Matrix product: dum = invalpha_i - Bi_mat               
      DO col = 1,r
         DO row = 1,r
            dum_mat(row,col) = 0.d0
         ENDDO
      ENDDO
          
      DO col = 1,r
         DO row = 1,r
            tmp = 0.d0
            DO k = 1,r
               !tmp = tmp + invalpha_i(row,k)*Bi_mat(k,pos1 + col)
               tmp = tmp + invalpha_i(row,k)*Bi_mat(k,pos + col)
            ENDDO
            dum_mat(row,col) = tmp
        ENDDO
      ENDDO 

      ! Matrix-vector product: dum_vec2 = dum*Yi_vec
      dum_vec2 = 0.d0
      DO col = 1,r
         tmp = Yi_vec(pos1 + col)
         DO row = 1,r
            dum_vec2(row) = dum_vec2(row) + dum_mat(row,col)*tmp 
         ENDDO 
      ENDDO 

      DO j = 1,r
         Yi_vec(pos + j) = dum_vec1(j) - dum_vec2(j) 
      ENDDO 

      ! Xn = Yn 
      DO i = 1,r
         pos    = r*(nr - 1) + i
         x(pos) = Yi_vec(pos)
      ENDDO

      ! Loop: backward substitution to find X (solution vector)      
      DO i = nr-1,1,-1

         ! Matrix vector product: dum_vec1 = Gammai_mat*x                
         dum_vec1 = 0.d0

         pos  = r*(i - 1)
         pos1 = pos + r
         DO col = 1,r
            tmp = x(pos1 + col)
            DO row = 1,r
               dum_vec1(row) = dum_vec1(row) + Gammai_mat(row,pos + col)*tmp 
            ENDDO 
         ENDDO 

         DO j = 1,r
            x(pos + j) = Yi_vec(pos + j) - dum_vec1(j)
         ENDDO 

      ENDDO        

    END SUBROUTINE tridiag_block_solver

    !------------------------------------------------------!
    !> This subroutine invrerts the matrix a. The output b is the inverted matrix.
    SUBROUTINE inv_matrix (r, a, b) 

      IMPLICIT NONE
      INTEGER  :: i, j
      REAL(KIND=8) :: d

      INTEGER, INTENT(IN) :: r
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: a
      REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: b

      INTEGER, DIMENSION(r) :: indx
 
      ! Matrix b initialization 
      DO j = 1,r 

         DO i = 1,j-1
            b(i,j) = 0.d0
         ENDDO

         b(j,j) = 1.d0

         DO i = j+1,r
            b(i,j) = 0.d0
         ENDDO

      ENDDO

      ! LU decomposition 
      CALL lu_decmp (a, indx, d)
      
      ! Columns of inverse matrix b = inv(a)
      DO j = 1,r
         CALL lu_solver (a, indx, b(1:r,j)) 
      ENDDO  

    END SUBROUTINE inv_matrix
    !--------------------------------------------------------!
    !> This subroutine performs the LU decomposition of the a matrix by means of the Crout's algorithm. 
    !! In order to save memory the decomposition is stored in the original matrix a.    
    SUBROUTINE lu_decmp (a, indx, d)

      IMPLICIT NONE
      INTEGER i, j, k, i_max, np
      REAL(KIND=8) :: toll, aamax, sum, dum
  
      REAL(KIND=8), INTENT(OUT) :: d
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: a
      INTEGER, DIMENSION(:), INTENT(OUT) :: indx
  
      REAL(KIND=8), DIMENSION(SIZE(a,1)) :: vv 
  
      ! Useful parameters
      d = 1.d0 ; toll = 1.d-20
        
      ! Loop over the rows to get the implicit scaling information
      np = SIZE(a,1) 
      DO i = 1,np 
         aamax = 0.d0
         DO j = 1,np
            IF (ABS(a(i,j)).GT.aamax) THEN 
               aamax = ABS(a(i,j))
            ENDIF
         ENDDO
         IF (aamax.EQ.0.d0) THEN
            PRINT*
            WRITE(*,10)'Sinlgular matrix in LU decomposition...'
            PRINT*,a
            PRINT*
            STOP
         ENDIF
         vv(i) = 1.d0/aamax     ! scaling is saved
      ENDDO 
  
      ! Loop over columns for Crout's method.
      DO j = 1,np
         DO i = 1,j - 1
            sum = a(i,j)
            DO k = 1,i - 1 
               sum = sum - a(i,k)*a(k,j)
            ENDDO
            a(i,j) = sum
         ENDDO 
         aamax = 0.d0
         DO i = j,np
            sum = a(i,j)
            DO k = 1,j - 1
               sum = sum - a(i,k)*a(k,j)
            ENDDO
            a(i,j) = sum 
            dum = vv(i)*ABS(sum)
            IF (dum.GE.aamax) THEN 
               i_max = i
               aamax = dum
            ENDIF  
         ENDDO
         IF (j.NE.i_max) THEN 
            DO k = 1,np 
               dum = a(i_max,k)
               a(i_max,k) = a(j,k)
               a(j,k) = dum
            ENDDO
            d = - d
            vv(i_max) = vv(j)
         ENDIF
         indx(j) = i_max
         IF (a(j,j).EQ.0.d0) THEN 
            a(j,j)= toll
         ENDIF
         IF (j.NE.np) THEN 
            dum = 1.d0/a(j,j)
            DO i = j + 1,np 
               a(i,j) = a(i,j)*dum
            ENDDO
         ENDIF
      ENDDO
  
10  FORMAT(A)

    END SUBROUTINE lu_decmp

    !-----------------------------------------------------------------------
    !> This subroutine computes the solution of the algebraic system of 
    !! equation once the LU decomposition of the A matrix has been 
    !! perfomed. The solution is found by means backward and forward 
    !! substitution on matrices L and U respectively
    SUBROUTINE lu_solver (a, indx, b)

    IMPLICIT NONE
 
    INTEGER i, j, ii, ll, np
    REAL(KIND=8) :: sum

    INTEGER, DIMENSION(:), INTENT(IN) :: indx(:)
    REAL(KIND=8), DIMENSION(:,:), INTENT(IN) :: a
    REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: b  

      ! Data inizialization 
      sum = 0.d0 
  
      ! Forward substitution
      np = SIZE(a,1)
      ii = 0   
      DO i = 1,np
         ll = indx(i)
         sum = b(ll)
         b(ll) = b(i) 
         IF (ii.NE.0) THEN 
            DO j = ii,i - 1 
               sum = sum - a(i,j)*b(j)
            ENDDO
         ELSEIF (sum.NE.0.d0) THEN 
            ii = i
         ENDIF
         b(i) = sum
      ENDDO

      ! Backward substitution
      DO i = np,1,-1
         sum = b(i)
         DO j = i + 1,np
            sum = sum - a(i,j)*b(j)
         ENDDO
         b(i) = sum/a(i,i)
      ENDDO

    END SUBROUTINE lu_solver   

#ifdef GOTO_BLAS
    !----------------------------------------------------!
    !> This subroutine solves the tridiagonal block system by means of
    !! the generalized Thomas algotithm. Use is made of the GOTO BLAS library 
    !! in order to speed the resolution of the algebraic system
    SUBROUTINE tridiag_block_solver_GOTO_BLAS (r, nr, d, x)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: r, nr
      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: d
      REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: x
        
      INTEGER :: i, j, k, row, col
      INTEGER :: pos, pos1, r2, info
      REAL(KIND=8) :: tmp

      ! Useful quantity 
      r2 = r*r

      ! alphai_i matrix
      DO j = 1,r
         DO i = 1,r
            invalpha_i(i,j) = Ai_mat(i,j)
         ENDDO
      ENDDO

      ! Generalized Thomas algorithm for block tridiagonal matrix

      ! Inversion of alpha_i matrix
      CALL DGETRF(r, r, invalpha_i, r, ipiv_GOTO_BLAS, info)
      CALL DGETRI(r, invalpha_i, r, ipiv_GOTO_BLAS, work_GOTO_BLAS, r2, info)

      ! Matrix-Matrix product: Gammai_mat = invalpha_i*Ci_mat  
      CALL DGEMM('n', 'n', r, r, r, 1.d0, invalpha_i, r, Ci_mat(1:r,1:r), r, 0.d0, Gammai_mat(1:r,1:r), r) 

      ! Matrix-vector product: Yi_vec = invalpha_i*d           
      CALL DGEMV ('n', r, r, 1.d0, invalpha_i, r, d(1:r), 1, 0.d0, Yi_vec(1:r), 1)

      ! Loop: forward substitution to find Yi_vec
      DO i = 2,nr - 1

         ! Matrix-Matrix product: dum = Bi_mat*Gammai_mat      
         pos  = r*(i - 1)
         pos1 = r*(i - 2) 

         ! Matrix difference: alpha_i = Ai_mat - dum is added in the loop   
         invalpha_i = Ai_mat(1:r,pos + 1:pos + r)
         CALL DGEMM('n', 'n', r, r, r, -1.d0, Bi_mat(1:r,pos + 1:pos + r), r, & 
        &           Gammai_mat(1:r,pos1 + 1:pos1 + r), r, 1.d0, invalpha_i, r)

         ! Inversion of alpha_i matrix
         CALL DGETRF(r, r, invalpha_i, r, ipiv_GOTO_BLAS, info)
         CALL DGETRI(r, invalpha_i, r, ipiv_GOTO_BLAS, work_GOTO_BLAS, r2, info)

         ! Matrix-Matrix product: Gammai_mat = invalpha_i*Ci_mat        !! GOTO_BLAS TO BE USED HERE !!
         CALL  DGEMM('n', 'n', r, r, r, 1.d0, invalpha_i, r, Ci_mat(1:r,pos + 1:pos + r), r, & 
        &            0.d0, Gammai_mat(1:r,pos + 1:pos + r), r)  

         ! Matrix-vector product: dum_vec1 = invalpha_i*d                   
         CALL DGEMV ('n', r, r, 1.d0, invalpha_i, r, d(pos + 1:pos + r), 1, 0.d0, dum_vec1, 1) 

         ! Matrix-Matrix product: dum = invalpha_i - Bi_mat            
         CALL DGEMM('n', 'n', r, r, r, 1.d0, invalpha_i, r, Bi_mat(1:r,pos + 1:pos + r), r, 0.d0, dum_mat, r)

         ! Matrix-vector product: dum_vec2 = dum*Yi_vec          
         pos1 = r*(i - 2)
         CALL DGEMV ('n', r, r, 1.d0, dum_mat, r, Yi_vec(pos1 + 1:pos1 + r), 1, 0.d0, dum_vec2, 1) 

         DO j = 1,r
            Yi_vec(pos + j) = dum_vec1(j) - dum_vec2(j) 
         ENDDO 
 
      ENDDO
      
      ! Last step of forward loop
      i = nr

      pos  = r*(i - 1)
      pos1 = pos - r  

      ! Matrix-Matrix product: alpha_i = Ai_mat - Bi_mat*Gammai_mat
      invalpha_i = Ai_mat(1:r,pos + 1:pos + r)
      CALL DGEMM('n', 'n', r, r, r, -1.d0, Bi_mat(1:r,pos + 1:pos + r), r, & 
        &           Gammai_mat(1:r,pos1 + 1:pos1 + r), r, 1.d0, invalpha_i, r)
                       
      ! Inversion of alpha_i matrix 
      CALL DGETRF(r, r, invalpha_i, r, ipiv_GOTO_BLAS, info)
      CALL DGETRI(r, invalpha_i, r, ipiv_GOTO_BLAS, work_GOTO_BLAS, r2, info)

      ! Matrix-vector product: dum_vec1 = invalpha_i*d                    
      CALL DGEMV ('n', r, r, 1.d0, invalpha_i, r, d(pos + 1:pos + r), 1, 0.d0, dum_vec1, 1)

      ! Matrix-Matrix product: dum = invalpha_i - Bi_mat               
      CALL DGEMM('n', 'n', r, r, r, 1.d0, invalpha_i, r, Bi_mat(1:r,pos + 1:pos + r), r, 0.d0, dum_mat, r)

      ! Matrix-vector product: dum_vec2 = dum*Yi_vec
      CALL DGEMV ('n', r, r, 1.d0, dum_mat, r, Yi_vec(pos1 + 1:pos1 + r), 1, 0.d0, dum_vec2, 1)

      DO j = 1,r
         Yi_vec(pos + j) = dum_vec1(j) - dum_vec2(j) 
      ENDDO 

      ! Xn = Yn 
      DO i = 1,r
         pos    = r*(nr - 1) + i
         x(pos) = Yi_vec(pos)
      ENDDO

      ! Loop: backward substitution to find X (solution vector)      
      DO i = nr-1,1,-1

         ! Matrix vector product: dum_vec1 = Gammai_mat*x               
         pos  = r*(i - 1)
         pos1 = pos + r

         CALL DGEMV ('n', r, r, 1.d0, Gammai_mat(1:r,pos + 1:pos + r), r, x(pos1 + 1:pos1 + r), 1, 0.d0, dum_vec1, 1)  

         DO j = 1,r
            x(pos + j) = Yi_vec(pos + j) - dum_vec1(j)
         ENDDO 

      ENDDO        

    END SUBROUTINE tridiag_block_solver_GOTO_BLAS

    !------------------------------------------------------!
    !> This subroutine inverts the matrix by using the GOTO BLAS library
    !! (the vectors ipiv and work need to be pre-allocated in order to 
    !! make use of this subroutine) 
    SUBROUTINE inv_matrix_GOTO_BLAS (n, a)

      IMPLICIT NONE

      INTEGER :: info

      INTEGER, INTENT(IN) :: n
      REAL(KIND=8), DIMENSION(:,:), INTENT(INOUT) :: a

      ! LU decomposition
      CALL DGETRF(n, n, a, n, ipiv_GOTO_BLAS, info)
    
      ! Inverse after performing LU factorization
      CALL DGETRI(n, a, n, ipiv_GOTO_BLAS, work_GOTO_BLAS, n*n, info)
      
    END SUBROUTINE inv_matrix_GOTO_BLAS 
#endif
 
  END MODULE mod_algebra
!------------------------------------------------------------------------------!
