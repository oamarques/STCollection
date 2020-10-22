SUBROUTINE SSTECHKRSLT( N, D, E, M, W, Z, LDZ, HBANDR, WORK, RESULT, TASK, INFO )
!
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: HBANDR, INFO, LDZ, M, N
!
!.. Array Arguments ..
LOGICAL :: TASK( 3 )
REAL( KIND=PREC ) :: D( * ), E( * ), RESULT( * ), Z( LDZ, * ), W( * ), WORK( * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTECHKRSLT checks the results of the LAPACK symmetric tridiagonal          !
!  eigensolvers. The following tests are performed:                            !
!  
!  If TASK(1)=.TRUE. then                                                      !
!     RESULT( 1 ) = | T - Z*W*Z' | / ( |T|*N*ULP ) if M = N                    !
!                 = | Z'*T*Z - W | / ( |T|*N*ULP ) otherwise                   !
!     RESULT( 2 ) = | I - Z*Z' | / ( N*ULP ) if M = N                          !
!                 = | I - Z'*Z | / ( N*ULP ) otherwise                         !
!  If TASK(2)=.TRUE. then                                                      !
!     RESULT( 6 ) = max | ||u|| - 1/sqrt(2) | / ULP                            !
!     RESULT( 7 ) = max | ||v|| - 1/sqrt(2) | / ULP                            !
!  If TASK(3)=.TRUE. then                                                      !
!     RESULT( 3 ) = | B - U*S*V' | / ( |B|*N*ULP ) if M=N                      !
!                 = | U'*B*V - S | / ( |B|*N*ULP ) otherwise                   !
!     RESULT( 4 ) = | I - U*U' | / ( N*ULP ) if M=N                            !
!                 = | I - U'*U | / ( N*ULP ) otherwise                         !
!     RESULT( 5 ) = | I - V*V' | / ( N*ULP ) if M=N                            !
!                 = | I - V'*V | / ( N*ULP ) otherwise                         !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          The dimension of the matrix.                                        !
!                                                                              !
!  D       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The N diagonal elements of the tridiagonal matrix.                  !
!                                                                              !
!  E       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The (N-1) off-diagonal elements of the tridiagonal matrix in        !
!          elements 1 to N-1, E(N) is not used.                                !
!                                                                              !
!  M       (input) INTEGER                                                     !
!          The total number of eigenvalues found, 0 <= M <= N.                 !
!                                                                              !
!  W       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The first M elements contain the selected eigenvalues in            !
!          ascending order.                                                    !
!                                                                              !
!  Z       (input) REAL( KIND=PREC ) array, dimension ( LDZ, M )               !
!          The first M columns of Z contain the orthonormal eigenvectors of    !
!          the tridiagonal matrix corresponding to the selected eigenvalues,   !
!          with the i-th column of Z holding the eigenvector associated        !
!          with W(i).                                                          !
!                                                                              !
!  LDZ     (input) INTEGER                                                     !
!          The leading dimension of the array Z.                               !
!                                                                              !
!  HBANDR  (input) INTEGER                                                     !
!          Sets the halfbandwidth of the matrices used in the tests            !
!          RESULT( 1 ) = | Z'*T*Z - W | / ( |T|*N*ULP )                        !
!          RESULT( 2 ) = | I - Z'*Z | / ( N*ULP )                              !
!          i.e. max(1,N*(HBANDR/100)) subdiagonals of Z'*T*Z and Z'*Z          !
!          are computed. This option should be used mainly for very            !
!          large problems in order to save time.                               !
!                                                                              !
!  WORK    (workspace) REAL( KIND=PREC ) array, dimension ( (N+1)*N )          !
!          Workspace. (TGK tests need ((N/2)^2)*4 + N/2 < (N+1)*N.)            !
!                                                                              !
!  RESULT  (output) REAL( KIND=PREC ) array, dimension ( 7 )                   !
!          The values computed by the two tests described above. The           !
!          values are currently limited to 1/ulp, to avoid overflow.           !
!                                                                              !
!  TASK    (input) LOGICAL array, dimension ( 3 )                              !
!          Activates the calculations as described in the leading comments,    !
!          TASK(1): residuals and orthogonality level for the eigenvalue       !
!                   problem T*Z = Z*W                                          !
!          TASK(2): norm of the vectors u and v obtained from the eigenvalue   !
!                   problem T*Z = Z*W (i.e. TGK mode)                          !
!          TASK(3): residuals and orthogonality level for the singular         !
!                   value decomposition B = U*S*V', where U, S and V           !
!                   are obtained from the eigenvalues and vectors of           !
!                   [ 0 B', B 0 ] (i.e. SVD mode)                              !
!                                                                              !
!  INFO    (input) INTEGER                                                     !
!          Exit status from routine used to compute W and Z.                   !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, IWRKU, IWRKV, IWRK1, IWRK2, J, K, L, MB, MINMN, NCOLS, NDIAG, NB
REAL( KIND=PREC ) :: AUKJ, BNORM, MXNRMU, MXNRMV, NRMU, NRMV, RNORM, &
                     RSQRT2, SQRTULP, SUM, TNORM, ULP, UNFL
!
!.. External Subroutines ..
LOGICAL, EXTERNAL           :: SISNAN
REAL( KIND=PREC ), EXTERNAL :: SLAMCH, SLANGE, SLANSB, SLANST, SLANSY, SNRM2 
EXTERNAL SCOPY, SGEMV, SLASET, SSCAL, SSYR, SSYRK
!
!.. Intrinsic Functions ..
INTRINSIC INT, MAX, MIN, REAL, SQRT
!
!.. Executable Statements ......................................................
!
RESULT( 1:7 ) = ZERO
!
! Early return.           
!
IF ( HBANDR == 0 .OR. INFO /= 0 .OR. M == 0 ) RETURN
!
NCOLS = MIN( M, N )
NDIAG = MAX( 1, MIN( M, INT( N*(HBANDR/HNDRD) ) ) )
MINMN = N  !** MIN( M, N ) can be pessimistic **!
!
! Get machine parameters. 
!
ULP = SLAMCH( 'Precision' )
UNFL = SLAMCH( 'Safe minimum' )
SQRTULP = SQRT( ULP )
!
IF ( TASK( 1 ) ) THEN
!
! Compute the norm of T (M or 1, set at installation time).
!
TNORM = SLANST( '1', N, D, E )
TNORM = MAX( TNORM, UNFL )
!
! Check the residuals.
!
IF      ( NDIAG /= N ) THEN
!
!       | Z'*T*Z - W | / ( |T|*N*ULP )         
!      
        L = N + 1
        DO I = 1, NCOLS
           K = MIN( NDIAG, NCOLS-I+1 ) 
           DO J = 1,N
              WORK( J ) = D( J )*Z( J,I )
           END DO
           DO J = 1,N-1
              WORK( J ) = WORK( J ) + E( J )*Z( J+1,I )
           END DO
           DO J = 2,N
              WORK( J ) = WORK( J ) + E( J-1 )*Z( J-1,I )
           END DO
           CALL SGEMV( 'T', N, K, ONE, Z( 1,I ), LDZ, WORK, 1, &
                       ZERO, WORK( L ), 1 )
           WORK( L ) = WORK( L ) - W( I )
           L = L + NDIAG
        END DO
        RNORM = SLANSB( 'M', 'L', NCOLS, NDIAG-1, WORK( N+1 ), NDIAG, WORK )
!
ELSE IF ( M == N ) THEN
!
!       | T - Z*W*Z' | / ( |T|*N*ULP )
!
        CALL SLASET( 'Full', N, N, ZERO, ZERO, WORK, N )
        DO I = 1, N-1
           J = ( N+1 )*( I-1 )
           WORK( J+1 ) = D( I )
           WORK( J+2 ) = E( I )
        END DO
        WORK( N**2 ) = D( N )
        DO I = 1, N
           CALL SSYR( 'L', N, -W( I ), Z( 1,I ), 1, WORK, N )
        END DO
        RNORM = SLANSY( 'M', 'L', N, WORK, N, WORK( N**2+1 ) )
!
ELSE
!
!       | Z'*T*Z - W | / ( |T|*N*ULP )
!
        CALL SLASET( 'Full', M, M, ZERO, ZERO, WORK, M )
        DO I = 1, M
           DO J = 1, I
              L = M*( J-1 ) + I
              DO K = 1, N
                 AUKJ = D( K )*Z( K, J )
                 IF ( K /= 1 ) AUKJ = AUKJ + E( K-1 )*Z( K-1, J )
                 IF ( K /= N ) AUKJ = AUKJ + E( K )*Z( K+1, J )
                 WORK( L ) = WORK( L ) + Z( K, I )*AUKJ
              END DO
           END DO
           L = M*( I-1 ) + I
           WORK( L ) = WORK( L ) - W( I )
        END DO
        RNORM = SLANSY( 'M', 'L', M, WORK, M, WORK( M**2+1 ) )
!
END IF
!
IF ( SISNAN( RNORM ) ) THEN
   RESULT( 1 ) = RNORM
ELSE
   IF ( TNORM > RNORM ) THEN
      RESULT( 1 ) = ( RNORM / TNORM ) / ( MINMN*ULP )
   ELSE
      IF ( TNORM < ONE ) THEN
         RESULT( 1 ) = ( MIN( RNORM, MINMN*TNORM ) / TNORM )/ ( MINMN*ULP )
      ELSE
         RESULT( 1 ) = MIN( RNORM / TNORM, REAL( MINMN,PREC ) ) / ( MINMN*ULP )
      END IF
   END IF
   RESULT( 1 ) = MAX( RESULT( 1 ), SQRTULP )
END IF
!
! Check the orthogonality of the computed eigenvectors.
!
IF      ( NDIAG /= N ) THEN
!
!       | I - Z'*Z | / ( N*ULP )        

        L = N + 1
        DO I = 1, NCOLS
           K = MIN( NDIAG, NCOLS-I+1 ) 
           CALL SGEMV( 'T', N, K, ONE, Z( 1,I ), LDZ, Z( 1,I ), 1, &
                       ZERO, WORK( L ), 1 )
           WORK( L ) = ONE - WORK( L )
           L = L + NDIAG
        END DO
        RNORM = SLANSB( 'M', 'L', NCOLS, NDIAG-1, WORK( N+1 ), NDIAG, WORK )
!
ELSE IF ( M == N ) THEN
!
!       | I - Z*Z' | / ( N*ULP ) 
!
        J = 1
        CALL SSYRK( 'L', 'N', N, N, ONE, Z, LDZ, ZERO, WORK, N )
        DO I = 1, N
           WORK( J ) = WORK( J ) - ONE
           J = J + N + 1
        END DO
        RNORM = SLANSY( 'M', 'L', N, WORK, N, WORK( N**2+1 ) )
!
ELSE
!
!       | I - Z'*Z | / ( N*ULP )
!
        J = I
        CALL SSYRK( 'L', 'T', M, N, ONE, Z, LDZ, ZERO, WORK, M )
        DO I = 1, M
           WORK( J ) = WORK( J ) - ONE
           J = J + M + 1
        END DO
        RNORM = SLANSY( 'M', 'L', M, WORK, M, WORK( M**2+1 ) )
!
END IF
!
IF ( SISNAN( RNORM ) ) THEN
   RESULT( 2 ) = RNORM
ELSE
   RESULT( 2 ) = MAX( MIN( REAL( MINMN,PREC ), RNORM ) / ( MINMN*ULP ), &
                      SQRTULP )
END IF
!
END IF !* TASK( 1 ) *!
!
!-------------------------------------------------------------------------------
! The following tests apply to the case TGK = [ 0 B, B' 0 ]. If (s,u,v) is a
! singular triplet of B with ||u||=||v||=1 then (+/-s,q) are eigenpairs of TGK 
! where ||q||=1 and q=P*(u' +/-v')/sqrt(2)=( v_1 u_1 v_2 u_2 ... v_n u_n )/ 
! sqrt(2). Therefore, if size(T)=N then size(B)=NB=N/2. 
!-------------------------------------------------------------------------------
!
NB = N/2
MB = MIN( M, NB )
!
IWRKU = 1               ! Copy of U (left singular vectors)
IWRKV = IWRKU + NB*MB   ! Copy of V (right singular vectors)
IWRK1 = IWRKV + NB*MB   ! Workspace (matrix products etc.)
IWRK2 = IWRK1 + MB*MB   ! Workspace (matrix products etc.)
MINMN = NB
!
! Copy U and V into work arrays. (Need loop because of LDZ and array stride.)
!
DO I = 0, MB-1 
   CALL SCOPY( NB, Z( 2, I+1 ), 2, WORK( IWRKU+NB*I ), 1 ) ! NB-by-MB (U)
   CALL SCOPY( NB, Z( 1, I+1 ), 2, WORK( IWRKV+NB*I ), 1 ) ! NB-by-MB (V) 
END DO
!
IF ( TASK( 2 ) ) THEN
!
! Normalize the singular vectors (mathematically, the norm should be sqrt(2)).
!
RSQRT2 = ONE / SQRT( TWO )
!
J = IWRKU
K = IWRKV
DO I = 1, MB
   NRMU = SNRM2( NB, WORK( J ), 1 )
   NRMV = SNRM2( NB, WORK( K ), 1 )
   MXNRMU = MAX( ABS( NRMU-RSQRT2 ), MXNRMU )
   MXNRMV = MAX( ABS( NRMV-RSQRT2 ), MXNRMV )
   IF ( NRMU == ZERO ) NRMU = ONE
   IF ( NRMV == ZERO ) NRMV = ONE
   CALL SSCAL( NB, ONE/NRMU, WORK( J ), 1 )
   CALL SSCAL( NB, ONE/NRMV, WORK( K ), 1 )   
   J = J + NB
   K = K + NB
END DO
!
RESULT( 6 ) = MXNRMU / ( NB*ULP )
RESULT( 7 ) = MXNRMV / ( NB*ULP )
!
END IF !* TASK( 2 ) *!
!
IF ( TASK( 3 ) ) THEN
!
! Compute the 1-norm of B.
!
BNORM = ABS( E( N-1 ) )
DO I = 1, N-2, 2
   SUM = ABS( E( I ) ) + ABS( E( I+1 ) )
   IF( BNORM<SUM ) BNORM = SUM
END DO
BNORM = MAX( BNORM, UNFL )
!
! Check the residuals.
!
IF ( MB == NB ) THEN
!
!  | B - U*S*V' | / ( |B|*NB*ULP )
!
!  Create B from the entries of ( E ).
   CALL SLASET( 'Full', NB, NB, ZERO, ZERO, WORK( IWRK1 ), NB )
   K = IWRK1
   DO I = 1, NB*2-1, 2
      WORK( K ) = E( I )
      K = K + NB
      WORK( K ) = E( I+1 )
      K = K + 1
   END DO
!  Compute U*S (i.e. scale the columns of U).
   K = IWRKU
   L = IWRK2
   DO I = 1, NB
      CALL SCOPY( NB, WORK( K ), 1, WORK( L ), 1 )
      CALL SSCAL( NB, W( I ), WORK( L ), 1 )
      K = K + NB
      L = L + NB
   END DO
!  Residual matrix.
   CALL SGEMM( 'N', 'T', NB, NB, NB, -ONE, WORK( IWRK2 ), NB, &
               WORK( IWRKV), NB, ONE, WORK( IWRK1 ), NB ) 
!  Compute the norm of the residual matrix; it is OK to overwrite WORK( IWRK2 ).
   RNORM = SLANGE( 'M', NB, NB, WORK( IWRK1 ), NB, WORK( IWRK2 ) )
!
ELSE
!
!  | U'*B*V - S | / ( |B|*N*ULP )
!
!  Compute C = B*V.
   K = IWRKV
   L = IWRK2
   DO I = 1, MB
      DO J = 1, NB*2-1, 2
         WORK( L ) = E( J )*WORK( K ) + E( J+1 )*WORK( K+1 )
         K = K + 1
         L = L + 1
      END DO
      WORK( L ) = E( NB*2-1 )*WORK( K )
   END DO
!  Compute U'*C and subtract S.
   CALL SGEMM( 'T', 'N', MB, MB, NB, ONE, WORK( IWRKU ), NB, &
               WORK( IWRK2 ), NB, ZERO, WORK( IWRK1 ), MB )
   J = IWRK1
   DO I = 1, MB
      WORK( J ) = WORK( J ) - W( I )
      J = J + MB + 1
   END DO     
!  Compute the norm of the residual matrix; it is OK to overwrite WORK( IWRK2 ).
   RNORM = SLANGE( 'M', MB, MB, WORK( IWRK1 ), MB, WORK( IWRK2 ) )
!
END IF
IF ( BNORM > RNORM ) THEN
   RESULT( 3 ) = ( RNORM / BNORM ) / ( MINMN*ULP )
ELSE
   IF ( BNORM < ONE ) THEN
      RESULT( 3 ) = ( MIN( RNORM, MINMN*BNORM ) / BNORM )/ ( MINMN*ULP )
   ELSE
      RESULT( 3 ) = MIN( RNORM / BNORM, REAL( MINMN,PREC ) ) / ( MINMN*ULP )
   END IF
END IF
RESULT( 3 ) = MAX( RESULT( 3 ), SQRTULP )
!
! Check the orthogonality of the left singular vectors.
!
IF ( MB == NB ) THEN
!
!  | I - U*U' | / ( NB*ULP ) 
!
   CALL SSYRK( 'L', 'N', NB, NB, ONE, WORK( IWRKU ), NB, ZERO, WORK( IWRK1 ), NB )
   J = IWRK1
   DO I = 1, NB
      WORK( J ) = WORK( J ) - ONE
      J = J + NB + 1
   END DO
   RNORM = SLANSY( 'M', 'L', NB, WORK( IWRK1 ), NB, WORK( IWRK2 ) )
!
ELSE
!
!  | I - U'*U | / ( NB*ULP )
!
   CALL SSYRK( 'L', 'T', MB, NB, ONE, WORK( IWRKU ), NB, ZERO, WORK( IWRK1 ), MB )
   J = IWRK1
   DO I = 1, MB
      WORK( J ) = WORK( J ) - ONE
      J = J + MB + 1
   END DO
   RNORM = SLANSY( 'M', 'L', MB, WORK( IWRK1 ), MB, WORK( IWRK2 ) )
!
END IF
!
RESULT( 4 ) = MIN( REAL( MINMN,PREC ), RNORM ) / ( MINMN*ULP )
RESULT( 4 ) = MAX( RESULT( 6 ), SQRTULP )
!
! Check the orthogonality of the right singular vectors.
!
IF ( MB == NB ) THEN
!
!  | I - V*V' | / ( NB*ULP ) 
!
   CALL SSYRK( 'L', 'N', NB, NB, ONE, WORK( IWRKV ), NB, ZERO, WORK( IWRK1 ), NB )
   J = IWRK1
   DO I = 1, NB
      WORK( J ) = WORK( J ) - ONE
      J = J + NB + 1
   END DO
   RNORM = SLANSY( 'M', 'L', NB, WORK( IWRK1 ), NB, WORK( IWRK2 ) )
!
ELSE
!
!  | I - V'*V | / ( NB*ULP )
!
   CALL SSYRK( 'L', 'T', MB, NB, ONE, WORK( IWRKV ), NB, ZERO, WORK( IWRK1 ), MB )
   J = IWRK1
   DO I = 1, MB
      WORK( J ) = WORK( J ) - ONE
      J = J + MB + 1
   END DO
   RNORM = SLANSY( 'M', 'L', MB, WORK( IWRK1 ), MB, WORK( IWRK2 ) )
!
END IF
!
RESULT( 5 ) = MIN( REAL( MINMN,PREC ), RNORM ) / ( MINMN*ULP )
RESULT( 5 ) = MAX( RESULT( 5 ), SQRTULP )
!
END IF !* TASK( 3 ) *!
!
END SUBROUTINE SSTECHKRSLT
