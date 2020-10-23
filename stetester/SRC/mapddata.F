!
!*******************************************************************************
!
MODULE MAPDDATA
    CONTAINS
!-------------------------------------------------------------------------------
    SUBROUTINE DB_TO_DT( N, A, B, ALPHA, BETA )
!
!   Map bidiagonal representation into a tridiagonal representation.
!
    USE DSTEDEFINITIONS
    INTEGER :: I, N
    REAL( KIND=PREC ) :: A( N ), ALPHA( N ), B( N ), BETA( N )
    ALPHA( 1 ) = A( 1 )**2
    DO I = 1,N-1
       BETA( I ) = A( I )*B( I )
       ALPHA( I+1 ) = A( I+1 )**2 + B( I )**2
    END DO
    BETA( N ) = ZERO
    END SUBROUTINE DB_TO_DT
!-------------------------------------------------------------------------------
    SUBROUTINE DB_TO_DZ( N, A, B, Q, E )
!
!   Map bidiagonal representation into q's and e's.
!
!   Let
!       | a_1  b_1           |
!   B = |      a_2  b_2      |
!       |           a_2  b_3 |
!       |                a_4 |
!   then B'*B = 
!       | (a_1)^2  a_1*b_1                                           |
!   T = | a_1*b_1  (a_2)^2+(b_1)^2  a_2*b_2                          |
!       |          a_2*b_2          (a_3)^2+(b_2)^2  a_3*b_3         |
!       |                           a_3*b_3          (a_4)^2+(b_3)^2 |
!   which is diagonally similar to (see the characteristic polynomial)
!        | (a_1)^2      1                                                 |
!   T1 = | (a_1*b_1)^2  (a_2)^2+(b_1)^2  1                                |
!        |              (a_2*b_2)^2      (a_3)^2+(b_2)^2  1               |
!        |                               (a_3*b_3)^2      (a_4)^2+(b_3)^2 |
!   and then the mapping into q's and e's implemented below can be easily 
!   obtained, see the mapping T_TO_Z.
!
    USE DSTEDEFINITIONS
    INTEGER :: I, N
    REAL( KIND=PREC ) :: A( N ), B( N ), E( N ), Q( N )
    DO I = 1,N-1
       Q( I ) = A( I )**2
       E( I ) = B( I )**2
    END DO
    Q( N ) = A( N )**2
    E( N ) = ZERO
    END SUBROUTINE DB_TO_DZ
!-------------------------------------------------------------------------------
    SUBROUTINE DT_TO_DZ( N, ALPHA, BETA, Q, E )
!
!   Map tridiagonal representation into q's and e's.
!
!   Let
!       | alpha_1  beta_1                    |
!   T = | beta_1   alpha_2  beta_2           |
!       |          beta_2   alpha_3  beta_3  |
!       |                   beta_3   alpha_4 |
!   which is diagonally similar to (see the characteristic polynomial)
!        | alpha_1     1                               |
!   T1 = | (beta_1)^2  alpha_2     1                   |
!        |             (beta_2)^2  alpha_3     1       |
!        |                         (beta_3)^2  alpha_4 |
!   Let (this is the Z notation)
!       | q_1 1           |         | 1             |
!   U = |     q_2 1       | and U = | e_1 1         |
!       |         q_3 1   |         |     e_2 1     |
!       |             q_4 |         |         e_3 1 |
!   then
!         | q_1      1                         |
!   L*U = | e_1*q_1  e_1+q_2  1                |
!         |          e_2*q_2  e_2+q_3  1       |
!         |                   e_3*q_3  e_3+q_4 |
!   By writting T1=L*U we obtain the mapping implemented below.
!
    USE DSTEDEFINITIONS
    INTEGER :: I, N
    REAL( KIND=PREC ) :: ALPHA( N ), BETA( N ), E( N ), Q( N )
    Q( 1 ) = ALPHA( 1 )
    DO I = 1,N-1
       E( I ) = ( BETA( I ) / Q( I ) ) * BETA( I )
       Q( I+1 ) = ALPHA( I+1 ) - E( I )
       Q( I+1 ) = MAX( MICRO, Q( I+1 ) )
    END DO
    E( N ) = ZERO
    END SUBROUTINE DT_TO_DZ
!-------------------------------------------------------------------------------
    SUBROUTINE DZ_TO_DB( N, A, B, Q, E )
!
!   Map q's and e's into bidiagonal representation.
!
    USE DSTEDEFINITIONS
    INTEGER :: I, N
    REAL( KIND=PREC ) :: A( N ), B( N ), E( N ), Q( N )
    A( 1 ) = SQRT( Q( 1 ) )
    DO I = 1,N-1
       B( I ) = SQRT( E( I ) )
       A( I+1 ) = SQRT( Q( I+1 ) )
    END DO
    B( N ) = ZERO
    END SUBROUTINE DZ_TO_DB
!-------------------------------------------------------------------------------
    SUBROUTINE DSQTSHIFT( N, ALPHA, BETA )
!
!   Shift the matrix (using Gershgorin bounds) if necessary.
!
    USE DSTEDEFINITIONS
    INTEGER :: I, N
    REAL( KIND=PREC ) :: ALPHA( N ), BETA( N ), GL, GR
    GL = ALPHA( 1 ) - ABS( BETA(1) )
    GR = ALPHA( 1 ) + ABS( BETA(1) )
    DO I = 2,N
       GL = MIN( GL, ALPHA( I ) - ( ABS( BETA( I-1 ) ) + ABS( BETA( I ) ) ) )
       GR = MAX( GR, ALPHA( I ) + ( ABS( BETA( I-1 ) ) + ABS( BETA( I ) ) ) )
    END DO
    IF ( GL < ZERO ) THEN
       IF ( ABS( GR ) < ABS( GL ) ) THEN
          DO I = 1,N
             ALPHA( I ) = GR - ALPHA( I )
          END DO
       ELSE
          GL = MIN( ZERO, GL )
          DO I = 1,N
             ALPHA( I ) = ALPHA( I ) - GL
          END DO
       END IF
    END IF
    END SUBROUTINE DSQTSHIFT
!-------------------------------------------------------------------------------
END MODULE MAPDDATA
