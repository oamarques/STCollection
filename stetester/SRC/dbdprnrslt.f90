SUBROUTINE DBDPRNRSLT( N, M, S, Z, RESULT, INFO, ITEST, ICASE, DUMP, &
                       IINTO, ILO, IUO, VLO, VUO )
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: INFO, ITEST, ICASE, M, N
INTEGER, INTENT( IN ), OPTIONAL :: IINTO, ILO, IUO
REAL( KIND=PREC ), INTENT( IN ), OPTIONAL :: VLO, VUO
!
!.. Array Arguments ..
LOGICAL :: DUMP( 8 )
REAL( KIND=PREC ) :: RESULT( * ), S( * ), Z( * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DBDPRNRSLT prints/checks the results of DBDSVDX SVD tests.                  !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N        (input) INTEGER                                                    !
!           The dimension of the bidiagonal matrix.                            !
!                                                                              !
!  M        (input) INTEGER                                                    !
!           The total number of singular values/vectors found, 0 <= M <= N.    !
!                                                                              !
!  S        (input) REAL( KIND=PREC ) array, dimension ( N )                   !
!           The first M elements contain the selected singular values in       !
!           descending order.                                                  !
!                                                                              !
!  Z        (input) REAL( KIND=PREC ) array, dimension ( N*2, M )              !
!           The first M columns of Z contain the singular vectors of the       !
!           corresponding to S.                                                !
!                                                                              !
!  RESULT   (input) REAL( KIND=PREC ) array, dimension ( 4 )                   !
!           The results of the tests computed in DBDCHKRSLT.                   !
!                                                                              !
!  INFO     (input) INTEGER                                                    !
!           Exit status from routine used to compute S and Z.                  !
!                                                                              !
!  ITEST    (input) INTEGER                                                    !
!           Identifies the routine used to compute the SVD.                    !
!                                                                              !
!  ICASE    (input) INTEGER                                                    !
!           Identifies the case associated to S and Z.                         !
!                                                                              !
!  DUMP     (input) LOGICAL array, dimension ( 8 )                             !
!           Defines data to be written into files,                             !
!           DUMP( 1 ) : tridiagonal matrix (i,d_i,e_i)                         !
!           DUMP( 2 ) : eigenvalues                                            !
!           DUMP( 3 ) : eigenvectors                                           !
!           DUMP( 4 ) : timing, residuals, orthogonality                       !
!           DUMP( 5 ) : tridiagonal matrix (i,d_i,e_i) in Matlab format        !
!           DUMP( 6 ) : eigenvalues in Matlab format                           !
!           DUMP( 7 ) : eigenvectors in Matlab format                          !
!           DUMP( 8 ) : singular values and vectors in Matlab format           !
!                       (see TGK mode)                                         !
!                                                                              !
!  IINTO    (input,optional) INTEGER                                           !
!           Identifies the interval associated to S and Z.                     !
!                                                                              !
!  ILO      (input,optional) INTEGER                                           !
!           Index of the smallest computed singular value.                     !
!                                                                              !
!  IUO      (input,optional) INTEGER                                           !
!           Index of the largest computed singular value.                      !
!                                                                              !
!  VLO      (input,optional) REAL( KIND=PREC )                                 !
!           Lower bound of the interval searched for singular values.          !
!                                                                              !
!  VUO      (input,optional) REAL( KIND=PREC )                                 !
!           Upper bound of the interval searched for singular values.          !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
CHARACTER( LEN=6 ) :: TAG
INTEGER :: J
!
!.. Local Array ..
CHARACTER( LEN=8 ), DIMENSION( 5 ) :: STEST  = (/ 'DBDSQR  ', &
                                                  'DBDSDC  ', &
                                                  'DBDSVDXA', &
                                                  'DBDSVDXI', &
                                                  'DBDSVDXV' /)
!
!.. Intrinsic Function ..
INTRINSIC TRIM
!
!.. Executable Statements ......................................................
!
IF ( INFO /= 0 ) RESULT( 1 ) = -INFO
!
SELECT CASE ( ITEST )
!
CASE ( 1 )
!
!    DBDSQR
!
     IF ( DUMP( 4 ) ) WRITE( UNIT=FUDUMP( 4 ), &
                      FMT='(3X,2A,1P,E11.4,56X,2(A,E12.4E3),A)' ) &
                      STEST( ITEST ), ': TIME=', RESULT( 1 ), &
                      '(s_1=', S( 1 ), ', s_n=', S( M ), ')'
!
CASE ( 2 )
!
!    DBDSDC
!
     IF ( DUMP( 4 ) ) WRITE( UNIT=FUDUMP( 4 ), &
                      FMT='(3X,2A,1P,E11.4,56X,2(A,E12.4E3),A)' ) &
                      STEST( ITEST ), ': TIME=', RESULT( 1 ), &
                      '(s_1=', S( 1 ), ', s_n=', S( M ), ')'
!  
CASE ( 3 )
!
!    DBDSVDX( RANGE='A' )
!
     WRITE( UNIT=FUOUT, FMT='(3X,A,1P,4(A,E11.4))' ) &
            STEST( ITEST ), ': TIME=', RESULT( 1 ), &
            ', RESD=', RESULT( 2), ', ORTU=', RESULT( 3 ), &
            ', ORTV=', RESULT( 4 )
     IF ( DUMP( 4 ) ) WRITE( UNIT=FUDUMP( 4 ), &
                      FMT='(3X,A,1P,4(A,E11.4),2X,2(A,E12.4E3 ),A)' ) &
                      STEST( ITEST ), ': TIME=', RESULT( 1 ), &
                      ', RESD=', RESULT( 2), ', ORTU=', RESULT( 3 ), &
                      ', ORTV=', RESULT( 4 ), '(s_1=', S( 1 ), &
                      ', s_n=', S( M ), ')'
!
CASE ( 4 )
!
!    DBDSVDX( RANGE='I' )
!
     WRITE( UNIT=FUOUT, FMT='(3X,A,1P,4(A,E11.4),3(A,I5))' ) &
            STEST( ITEST ), ': TIME=', RESULT( 1 ), &
            ', RESD=', RESULT( 2 ), ', ORTU=', RESULT( 3 ), &
            ', ORTV=', RESULT( 4 ), ', M=', M, &
            ', IL=', ILO, ', IU=', IUO
     IF ( DUMP( 4 ) ) WRITE( UNIT=FUDUMP( 4 ), &
                      FMT='(3X,A,1P,4(A,E11.4),3(A,I5))' ) &
                      STEST( ITEST ), ': TIME=', RESULT( 1 ), &
                      ', RESD=', RESULT( 2 ), ', ORTU=', RESULT( 3 ), &
                      ', ORTV=', RESULT( 4 ), ', M=', M, &
                      ', IL=', ILO, ', IU=', IUO
!
CASE ( 5 )
!
!    DBDSVDX( RANGE='V' )
!
     WRITE( UNIT=FUOUT, FMT='(3X,A,1P,4(A,E11.4),A,I5,2(A,E12.4E3))' ) &
            STEST( ITEST ), ': TIME=', RESULT( 1 ), &
            ', RESD=', RESULT( 2 ), ', ORTU=', RESULT( 3 ), &
            ', ORTV=', RESULT( 4 ), ', M=', M, &
            ', VL=', VLO, ', VU=', VUO
     IF ( DUMP( 4 ) ) WRITE( UNIT=FUDUMP( 4 ), &
                      FMT='(3X,A,1P,4(A,E11.4),A,I5,2(A,E23.15E3))' ) &
                      STEST( ITEST ), ': TIME=', RESULT( 1 ), &
                      ', RESD=', RESULT( 2 ), ', ORTU=', RESULT( 3 ), &
                      ', ORTV=', RESULT( 4 ), ', M=', M, &
                      ', VL=', VLO, ', VU=', VUO
!     
END SELECT
!
IF ( ITEST.GE.3 .AND. DUMP( 8 ) ) THEN
   IF ( ITEST.EQ.3 ) THEN
      WRITE( TAG, FMT='(I3.3)' ) ICASE
   ELSE
      WRITE( TAG, FMT='(I3.3,A,I2.2)' ) ICASE, '_', IINTO
   END IF
   WRITE( UNIT=FUDUMP( 5 ), FMT='(4A)' ) &
          '% ', STEST( ITEST ), '_', TAG
   WRITE( UNIT=FUDUMP( 5 ), FMT='(2(A,I3.3),A,I5,A)' ) &
          'B = B_', ICASE, &
          '; NB = NB_', ICASE, &
          '; NS =', M, ';'
   WRITE( UNIT=FUDUMP( 5 ), FMT='(5A)' ) &
          'NS_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = NS;'
      WRITE( UNIT=FUDUMP( 5 ), FMT='(5A)' ) &
             'fprintf(''', STEST( ITEST ), '_', TAG, ', N =%4i:\n'',NB)'
   IF ( M==0 ) THEN
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) &
             'fprintf(''   NS = 0\n'')'
      WRITE( UNIT=FUDUMP( 5 ), FMT='(5A)' ) &
             'S_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = [];', &
             'U_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = [];', &
             'V_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = [];', &
             'clear B NB NS'
   ELSE
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) &
             'S = zeros(NS,1);'
      WRITE( UNIT=FUDUMP( 5 ), FMT='(2(A,I6,A,1P,E23.15E3,A))' ) &
             ( 'S(', J, ')=', S( J ), '; ', J = 1,M ) 
      WRITE( UNIT=FUDUMP( 5 ), FMT='(2(A,I6,A,1P,E23.15E3,A))' ) &
             ( 'Z_SVDX(', J, ')=', Z( J ), '; ', J = 1,N*2*M ) 
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) &
             'S = diag(S);', &
             'Z_SVDX = reshape(Z_SVDX,NB*2,NS);', &
             'U = Z_SVDX(1   :NB  ,1:NS);', &
             'V = Z_SVDX(1+NB:NB*2,1:NS);', &
             'normB = norm(B); if (normB==0), normB=1.0; end', &
             'normR = norm(U''*B*V-S)/(normB*NB*eps);', &
             'normU = norm(eye(NS)-U''*U)/(NB*eps);', &
             'normV = norm(eye(NS)-V''*V)/(NB*eps);', &
             'fprintf(''   ||U''''*B*V-S||/(norm(B)*N*eps) =%11.4e\n'',normR)', &
             'fprintf(''   ||I-U''''*U||/(N*eps) =%11.4e\n'',normU)', &
             'fprintf(''   ||I-V''''*V||/(N*eps) =%11.4e\n'',normV)'
      WRITE( UNIT=FUDUMP( 5 ), FMT='(5A)' ) &
             'S_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = diag(S);', &
             'U_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = U;', &
             'V_', TRIM( TAG ), '_', STEST( ITEST )(4:8), ' = V;', &
             'clear B NB NS S U V Z_SVDX normB normR normU normV'
   END IF
END IF
!
END SUBROUTINE DBDPRNRSLT
