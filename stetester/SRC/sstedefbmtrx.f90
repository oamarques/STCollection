SUBROUTINE SSTEDEFBMTRX( MTYPE, N, D, E, ISEED, WORK )
!
USE SSTEDEFINITIONS
!
!.. Scalar Arguments ..
INTEGER :: MTYPE, N
!
!.. Array Arguments ..
INTEGER :: ISEED( 4 )
REAL( KIND=PREC ) :: D( N ), E( N ), WORK( * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEDEFBMTRX generates a bidiagonal matrix as indicated by MTYPE.           !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  MTYPE  (input) INTEGER                                                      !
!         Sets the matrix type:                                                !
!         = 0, zero matrix                                                     !
!         = 1, identity matrix                                                 !
!         = 2, random bidiagonal matrix                                        !
!         = 3, random bidiagonal matrix                                        !
!                                                                              !
!  D      (output) REAL( KIND=PREC ) array, dimension ( N )                    !
!         The N diagonal elements of the bidiagonal matrix.                    !
!                                                                              !
!  E      (output) REAL( KIND=PREC ) array, dimension ( N )                    !
!         The (N-1) off-diagonal elements of the bidiagonal                    !
!         matrix in elements 1 to N-1, E(N) is set to 0.                       !
!                                                                              !
!  N      (input) INTEGER                                                      !
!         The order of the matrix. N > 0.                                      !
!                                                                              !
!  ISEED  (input/output) INTEGER array, dimension ( 4 )                        !
!         Seed for the random number generator. Each entry of ISEED should     !
!         lie between 0 and 4095 inclusive and ISEED(4) should be odd.         !
!                                                                              !
!  WORK   (workspace) REAL( KIND=PREC ) array, dimension ( N )                 !
!         Workspace.                                                           !
!                                                                              !
!==============================================================================!
!
!.. Local Scalar ..
INTEGER :: I
REAL( KIND=PREC ) :: GAMMA, ULP
!
!.. Local Array ..
INTEGER :: JSEED( 4 )
!
!.. External Functions ..
REAL( KIND=PREC ), EXTERNAL :: SLAMCH, SLARND, SNRM2
!
!.. Executable Statements ......................................................
!
SELECT CASE ( MTYPE )
!
CASE ( 0 )
!
!    Zero matrix
!
     D = ZERO
     E = ZERO
!
CASE ( 1 )
!
!    Identity matrix
!
     D = ONE
     E = ZERO
!
CASE ( 2 )
!
!    Random matrix 
!    a(j) = nrm2(rand(1:n-j+1,1))
!    b(j) = nrm2(rand(1:n-j,1))                    
!
     DO I = 1, N
        JSEED( 1 ) = ISEED( 1 ) ! I
        JSEED( 2 ) = ISEED( 2 ) ! I + 1
        JSEED( 3 ) = ISEED( 3 ) ! I + 2
        JSEED( 4 ) = 2*I + 3
        CALL SLARNV( 3, JSEED, N+1-I, WORK )
        D( I ) = SNRM2( N+1-I, WORK, 1 )
        JSEED( 1 ) = ISEED( 1 ) ! I
        JSEED( 2 ) = ISEED( 2 ) ! I + 1
        JSEED( 3 ) = ISEED( 3 ) ! I + 2
        JSEED( 4 ) = 2*I + 11
        CALL SLARNV( 3, JSEED, N-I, WORK )
        E( I ) = SNRM2( N-I, WORK, 1 )
     END DO
     E( N ) = ZERO
!
CASE ( 3 )
!
!    Random matrix                     
!    a(j) = e^x, x = rand[ 2*log(ulp), -2*log(ulp) ]
!    b(j) = e^x, x = rand[ 2*log(ulp), -2*log(ulp) ]
!
     ULP = SLAMCH( 'Precision' )
     GAMMA = -TWO*LOG( ULP )
     DO I = 1, N
        D( I ) = EXP( GAMMA*SLARND( 2, ISEED ) )
        E( I ) = EXP( GAMMA*SLARND( 2, ISEED ) )
     END DO
     E( N ) = ZERO
!
END SELECT
! 
END SUBROUTINE SSTEDEFBMTRX
