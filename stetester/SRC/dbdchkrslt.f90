SUBROUTINE DBDCHKRSLT( N, D, E, M, W, Z, LDZ, HBANDR, WORK, RESULT, INFO )
!
USE DSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: HBANDR, INFO, LDZ, M, N
!
!.. Array Arguments ..
REAL( KIND=PREC ) :: D( * ), E( * ), RESULT( * ), Z( LDZ, * ), W( * ), WORK( * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DBDCHKRSLT checks the results of the LAPACK bidiagonal SVD. The             !
!  following tests are performed:                                              !
!                                                                              !
!  RESULT( 1 ) = | B - U*S*V' | / ( |B|*N*ULP ) if M=N                         !
!              = | U'*B*V - S | / ( |B|*N*ULP ) otherwise                      !
!  RESULT( 2 ) = | I - U*U' | / ( N*ULP ) if M=N                               !
!              = | I - U'*U | / ( N*ULP ) otherwise                            !
!  RESULT( 3 ) = | I - V*V' | / ( N*ULP ) if M=N                               !
!              = | I - V'*V | / ( N*ULP ) otherwise                            !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          The dimension of the bidiagonal matrix.                             !
!                                                                              !
!  D       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The N diagonal elements of the bidiagonal matrix.                   !
!                                                                              !
!  E       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The (N-1) off-diagonal elements of the bidiagonal matrix in         !
!          elements 1 to N-1, E(N) is not used.                                !
!                                                                              !
!  M       (input) INTEGER                                                     !
!          The total number of singular values found, 0 <= M <= N.             !
!                                                                              !
!  W       (input) REAL( KIND=PREC ) array, dimension ( N )                    !
!          The first M elements contain the selected singular values in        !
!          descending order.                                                   !
!                                                                              !
!  Z       (input) REAL( KIND=PREC ) array, dimension ( LDZ, M )               !
!          The first M columns of Z contain the singular vectors of the        !
!          bidiagonal matrix corresponding to the selected singular values,    !
!          with U in rows 1 to N and V in rows N+1 to N*2, i.e.                !
!          Z = [ U ]                                                           !
!              [ V ]                                                           !
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
!  RESULT  (output) REAL( KIND=PREC ) array, dimension ( 3 )                   !
!          The values computed by the two tests described above. The           !
!          values are currently limited to 1/ulp, to avoid overflow.           !
!                                                                              !
!  INFO    (input) INTEGER                                                     !
!          Exit status from routine used to compute W and Z.                   !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
INTEGER :: I, IWRKU, IWRKV, IWRK1, IWRK2, J, K, L, MB, MINMN, NB
REAL( KIND=PREC ) :: BNORM, RNORM, SQRTULP, SUM, ULP, UNFL
!
!.. External Subroutines ..
REAL( KIND=PREC ), EXTERNAL :: DLAMCH, DLANGE, DLANSB, DLANST, DLANSY, DNRM2 
EXTERNAL DCOPY, DGEMV, DLASET, DSCAL, DSYR, DSYRK
!
!.. Intrinsic Functions ..
INTRINSIC INT, MAX, MIN, REAL, SQRT
!
!.. Executable Statements ......................................................
!
RESULT( 1:3 ) = ZERO
!
! Early return.           
!
IF ( HBANDR == 0 .OR. INFO /= 0 .OR. M == 0 ) RETURN
!
! Get machine parameters. 
!
ULP = DLAMCH( 'Precision' )
UNFL = DLAMCH( 'Safe minimum' )
SQRTULP = SQRT( ULP )
!
NB = N
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
   CALL DCOPY( NB, Z( 1   , I+1 ), 1, WORK( IWRKU+NB*I ), 1 ) ! NB-by-MB (U)
   CALL DCOPY( NB, Z( 1+NB, I+1 ), 1, WORK( IWRKV+NB*I ), 1 ) ! NB-by-MB (V)
END DO
!
! Compute the 1-norm of B.
!
BNORM = ABS( D( NB ) )
DO I = 1, NB-1
   SUM = ABS( D( I ) ) + ABS( E( I ) )
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
   CALL DLASET( 'Full', NB, NB, ZERO, ZERO, WORK( IWRK1 ), NB )
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
      CALL DCOPY( NB, WORK( K ), 1, WORK( L ), 1 )
      CALL DSCAL( NB, W( I ), WORK( L ), 1 )
      K = K + NB
      L = L + NB
   END DO
!  Residual matrix.
   CALL DGEMM( 'N', 'T', NB, NB, NB, -ONE, WORK( IWRK2 ), NB, &
               WORK( IWRKV), NB, ONE, WORK( IWRK1 ), NB ) 
!  Compute the norm of the residual matrix; it is OK to overwrite WORK( IWRK2 ).
   RNORM = DLANGE( 'M', NB, NB, WORK( IWRK1 ), NB, WORK( IWRK2 ) )
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
   CALL DGEMM( 'T', 'N', MB, MB, NB, ONE, WORK( IWRKU ), NB, &
               WORK( IWRK2 ), NB, ZERO, WORK( IWRK1 ), MB )
   J = IWRK1
   DO I = 1, MB
      WORK( J ) = WORK( J ) - W( I )
      J = J + MB + 1
   END DO     
!  Compute the norm of the residual matrix; it is OK to overwrite WORK( IWRK2 ).
   RNORM = DLANGE( 'M', MB, MB, WORK( IWRK1 ), MB, WORK( IWRK2 ) )
!
END IF
IF ( BNORM > RNORM ) THEN
   RESULT( 1 ) = ( RNORM / BNORM ) / ( MINMN*ULP )
ELSE
   IF ( BNORM < ONE ) THEN
      RESULT( 1 ) = ( MIN( RNORM, MINMN*BNORM ) / BNORM )/ ( MINMN*ULP )
   ELSE
      RESULT( 1 ) = MIN( RNORM / BNORM, REAL( MINMN,PREC ) ) / ( MINMN*ULP )
   END IF
END IF
RESULT( 1 ) = MAX( RESULT( 1 ), SQRTULP )
!
! Check the orthogonality of the left singular vectors.
!
IF ( MB == NB ) THEN
!
!  | I - U*U' | / ( NB*ULP ) 
!
   CALL DSYRK( 'L', 'N', NB, NB, ONE, WORK( IWRKU ), NB, ZERO, WORK( IWRK1 ), NB )
   J = IWRK1
   DO I = 1, NB
      WORK( J ) = WORK( J ) - ONE
      J = J + NB + 1
   END DO
   RNORM = DLANSY( 'M', 'L', NB, WORK( IWRK1 ), NB, WORK( IWRK2 ) )
!
ELSE
!
!  | I - U'*U | / ( NB*ULP )
!
   CALL DSYRK( 'L', 'T', MB, NB, ONE, WORK( IWRKU ), NB, ZERO, WORK( IWRK1 ), MB )
   J = IWRK1
   DO I = 1, MB
      WORK( J ) = WORK( J ) - ONE
      J = J + MB + 1
   END DO
   RNORM = DLANSY( 'M', 'L', MB, WORK( IWRK1 ), MB, WORK( IWRK2 ) )
!
END IF
!
RESULT( 2 ) = MIN( REAL( MINMN,PREC ), RNORM ) / ( MINMN*ULP )
RESULT( 2 ) = MAX( RESULT( 2 ), SQRTULP )
!
! Check the orthogonality of the right singular vectors.
!
IF ( MB == NB ) THEN
!
!  | I - V*V' | / ( NB*ULP ) 
!
   CALL DSYRK( 'L', 'N', NB, NB, ONE, WORK( IWRKV ), NB, ZERO, WORK( IWRK1 ), NB )
   J = IWRK1
   DO I = 1, NB
      WORK( J ) = WORK( J ) - ONE
      J = J + NB + 1
   END DO
   RNORM = DLANSY( 'M', 'L', NB, WORK( IWRK1 ), NB, WORK( IWRK2 ) )
!
ELSE
!
!  | I - V'*V | / ( NB*ULP )
!
   CALL DSYRK( 'L', 'T', MB, NB, ONE, WORK( IWRKV ), NB, ZERO, WORK( IWRK1 ), MB )
   J = IWRK1
   DO I = 1, MB
      WORK( J ) = WORK( J ) - ONE
      J = J + MB + 1
   END DO
   RNORM = DLANSY( 'M', 'L', MB, WORK( IWRK1 ), MB, WORK( IWRK2 ) )
!
END IF
!
RESULT( 3 ) = MIN( REAL( MINMN,PREC ), RNORM ) / ( MINMN*ULP )
RESULT( 3 ) = MAX( RESULT( 3 ), SQRTULP )
!
END SUBROUTINE DBDCHKRSLT
