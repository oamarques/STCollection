SUBROUTINE SSTEPRNRSLT( N, M, W, Z, LDZ, RESULT, INFO, ITEST, ICASE, &
                        DUMP, IINTO, ILO, IUO, VLO, VUO, GAPMINO )
!
USE GSTEDEFINITIONS
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
INTEGER :: ICASE, INFO, ITEST, LDZ, M, N
INTEGER, INTENT( IN ), OPTIONAL :: IINTO, ILO, IUO
REAL( KIND=PREC ), INTENT( IN ), OPTIONAL :: VLO, VUO, GAPMINO
!
!.. Array Arguments ..
LOGICAL :: DUMP( 8 )
REAL( KIND=PREC ) :: RESULT( 8 ), W( * ), Z( LDZ, * )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEPRNRSLT prints the results of the LAPACK symmetric tridiagonal          !
!  eigensolver tests.                                                          !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N        (input) INTEGER                                                    !
!           The dimension of the matrix.                                       !
!                                                                              !
!  M        (input) INTEGER                                                    !
!           The total number of eigenvalues found, 0 <= M <= N.                !
!                                                                              !
!  W        (input) REAL( KIND=PREC ) array, dimension ( N )                   !
!           The first M elements contain the selected eigenvalues in           !
!           ascending order.                                                   !
!                                                                              !
!  Z        (input) REAL( KIND=PREC ) array, dimension ( LDZ, M )              !
!           The first M columns of Z contain the orthonormal eigenvectors of   !
!           the tridiagonal matrix corresponding to the selected eigenvalues,  !
!           with the i-th column of Z holding the eigenvector associated       !
!           with W(i).                                                         !
!                                                                              !
!  LDZ      (input) INTEGER                                                    !
!           The leading dimension of the array Z.                              !
!                                                                              !
!  RESULT   (input) REAL( KIND=PREC ) array, dimension ( 3 )                   !
!           The results of the tests computed in SSTECHKRSLT.                  !
!                                                                              !
!  INFO     (input) INTEGER                                                    !
!           Exit status from routine used to compute S and Z.                  !
!                                                                              !
!  ITEST    (input) INTEGER                                                    !
!           Identifies the routine used to compute W and Z.                    !
!                                                                              !
!  ICASE    (input) INTEGER                                                    !
!           Identifies the case associated to W and Z.                         !
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
!           Identifies the interval associated to W and Z.                     !
!                                                                              !
!  ILO      (input,optional) INTEGER                                           !
!           Index of the smallest computed eigenvalue.                         !
!                                                                              !
!  IUO      (input,optional) INTEGER                                           !
!           Index of the largest computed eigenvalue.                          !
!                                                                              !
!  VLO      (input,optional) REAL( KIND=PREC )                                 !
!           Lower bound of the interval searched for eigenvalues.              !
!                                                                              !
!  VUO      (input,optional) REAL( KIND=PREC )                                 !
!           Upper bound of the interval searched for eigenvalues.              !
!                                                                              !
!  GAPMINO  (input,optional) REAL( KIND=PREC )                                 !
!           Minimum gap between | W(ILO-1)-W(ILO) | and | W(IUO-W(IUO+1) |     !
!           or between | W(right_of_VLO)-W(left_of_VLO) | and                  !
!           | W(right_of_VUO)-W(left_of_VUO) |.                                !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
CHARACTER( LEN=3 ) :: STEST
INTEGER :: FIZZLE = 543210, I, J
!
!.. Executable Statements ......................................................
!
! Standard output.
!
IF      ( PRESENT( ILO ) .AND. PRESENT( IUO ) ) THEN
        IF ( INFO == 0 ) THEN
           WRITE( UNIT=FUOUT, FMT='(3X,1P,A,3(A,E9.2),3(A,I5),A,E9.2)' ) &
                  IDTEST( ITEST ), ': TIME=', RESULT( 1 ), &
                  ', RESD=', RESULT( 2 ), ', ORTH=', RESULT( 3 ), &
                  ', M=', M, ', IL=', ILO, ', IU=', IUO, &
                  ', GAPMIN=', GAPMINO
        ELSE
           WRITE( UNIT=FUOUT, FMT='(3X,1P,2A,I9,2(A,I5))' ) &
                  IDTEST( ITEST ), ': INFO=', INFO, &
                  ', IL=', ILO, ', IU=', IUO
        END IF
ELSE IF ( PRESENT( VLO ) .AND. PRESENT( VUO ) ) THEN
        IF ( INFO == 0 ) THEN
           WRITE( UNIT=FUOUT, FMT='(3X,1P,A,3(A,E9.2),A,I5,3(A,E11.4))' ) &
                  IDTEST( ITEST ), ': TIME=', RESULT( 1 ), &
                  ', RESD=', RESULT( 2 ), ', ORTH=', RESULT( 3 ), &
                  ', M=', M, ', VL=', VLO, ', VU=', VUO, &
                  ', GAPMIN=', GAPMINO
        ELSE
           WRITE( UNIT=FUOUT, FMT='(3X,1P,2A,I9,2(A,E11.4))' ) &
                  IDTEST( ITEST ), ': INFO=', INFO, &
                  ', VL=', VLO, ', VU=', VUO
        END IF
ELSE
        IF ( INFO == 0 ) THEN
           WRITE( UNIT=FUOUT, FMT='(3X,1P,A,3(A,E9.2))' ) &
                  IDTEST( ITEST ), ': TIME=', RESULT( 1 ), &
                  ', RESD=', RESULT( 2 ), ', ORTH=', RESULT( 3 )
        ELSE
           WRITE( UNIT=FUOUT, FMT='(3X,1P,2A,I9)' ) &
                  IDTEST( ITEST ), ': INFO=', INFO
        END IF
END IF
!
! Print eigenvalues in appropriate file.
!
IF ( DUMP( 2 ) ) THEN
   IF      ( PRESENT( ILO ) .AND. PRESENT( IUO ) ) THEN
           WRITE( UNIT=FUDUMP( 2 ), FMT='(5X,A,2(A,I5))' ) &
                  IDTEST( ITEST ), ', IL=', ILO, ', IU=', IUO 
   ELSE IF ( PRESENT( VLO ) .AND. PRESENT( VUO ) ) THEN
           WRITE( UNIT=FUDUMP( 2 ), FMT='(5X,A,1P,2(A,E15.7E3))' ) &
                  IDTEST( ITEST ), ', VL=', VLO, ', VU=', VUO 
   ELSE
           WRITE( UNIT=FUDUMP( 2 ), FMT='(5X,A)' ) &
                  IDTEST( ITEST )
   END IF
   IF ( M == 0 ) THEN
      WRITE( UNIT=FUDUMP( 2 ), FMT='(8X,A)' ) &
             'M = 0, W = [ ]'
   ELSE
      WRITE( UNIT=FUDUMP( 2 ), FMT='(8X,A,I5,A,/,(3X,I5,1P,E17.7))' ) &
             'W( 1:', M, ' ) =', ( I, W( I ), I = 1,M )
   END IF
END IF
!
! Print eigenvectors in appropriate file.
!
IF ( DUMP( 3 ) ) THEN
   IF      ( PRESENT( ILO ) .AND. PRESENT( IUO ) ) THEN
           WRITE( UNIT=FUDUMP( 3 ), FMT='(5X,A,2(A,I5))' ) &
                  IDTEST( ITEST ), ', IL=', ILO, ', IU=', IUO 
   ELSE IF ( PRESENT( VLO ) .AND. PRESENT( VUO ) ) THEN
           WRITE( UNIT=FUDUMP( 3 ), FMT='(5X,A,1P,2(A,E11.4))' ) &
                  IDTEST( ITEST ), ', VL=', VLO, ', VU=', VUO 
   ELSE
           WRITE( UNIT=FUDUMP( 3 ), FMT='(5X,A)' ) IDTEST( ITEST )
   END IF
   IF ( M == 0 ) THEN
      WRITE( UNIT=FUDUMP( 3 ), FMT='(7X,A)' ) &
             'M = 0, Z = [ ]'
   ELSE
      DO I = 1, M
         WRITE( UNIT=FUDUMP( 3 ), FMT='(8X,A,I5,A,1P,/,(3X,E15.7E3))' ) &
                'Z( :,', I, ' ) =', ( Z( J, I ), J = 1,N )
      END DO
   END IF
END IF
!
! Print timing, residuals and orthogonality.
!
IF ( DUMP( 4 ) ) THEN
   IF      ( PRESENT( ILO ) .AND. PRESENT( IUO ) ) THEN
           IF ( INFO == 0 ) THEN
              WRITE( UNIT=FUDUMP( 4 ), FMT='(3X,1P,A,3(A,E11.4),3(A,I5),A,E11.4)' ) &
                     IDTEST( ITEST ), ': TIME=', RESULT( 1 ), &
                     ', RESD=', RESULT( 2 ), ', ORTH=', RESULT( 3 ), &
                     ', M=', M, ', IL=', ILO, ', IU=', IUO, &
                     ', GAPMIN=', GAPMINO
           ELSE
              WRITE( UNIT=FUDUMP( 4 ), FMT='(3X,1P,A,3(A,I11),3(A,I5),A)' ) &
                     IDTEST( ITEST ), ': TIME=', -INFO, ', RESD=', FIZZLE, &
                     ', ORTH=', FIZZLE, ', M=', 0, ', IL=', ILO, ', IU=', IUO, &
                     '  * INFO > 0 *'
           END IF
   ELSE IF ( PRESENT( VLO ) .AND. PRESENT( VUO ) ) THEN
           IF ( INFO == 0 ) THEN
              WRITE( UNIT=FUDUMP( 4 ), FMT='(3X,1P,A,3(A,E11.4),A,I5,2(A,E15.7E3),A,E11.4)' ) &
                     IDTEST( ITEST ), ': TIME=', RESULT( 1 ), &
                     ', RESD=', RESULT( 2 ), ', ORTH=', RESULT( 3 ), &
                     ', M=', M, ', VL= ', VLO, ', VU= ', VUO, &
                     ', GAPMIN=', GAPMINO
           ELSE
              WRITE( UNIT=FUDUMP( 4 ), FMT='(3X,1P,A,3(A,I11),A,I5,2(A,E15.7E3),A)' ) &
                     IDTEST( ITEST ), ': TIME=', -INFO, ', RESD=', FIZZLE, &
                     ', ORTH=', FIZZLE, ', M=', 0, ', VL= ', VLO, ', VU= ', VUO, &
                     '  * INFO > 0 *'
           END IF
   ELSE
           IF ( INFO == 0 ) THEN
              WRITE( UNIT=FUDUMP( 4 ), FMT='(3X,1P,A,3(A,E11.4))' ) &
                     IDTEST( ITEST ), ': TIME=', RESULT( 1 ), &
                     ', RESD=', RESULT( 2 ), ', ORTH=', RESULT( 3 )
           ELSE
              WRITE( UNIT=FUDUMP( 4 ), FMT='(3X,1P,A,3(A,I11),A)' ) &
                     IDTEST( ITEST ), ': TIME=', -INFO, ', RESD=', FIZZLE, &
                     ', ORTH=', FIZZLE, '  * INFO > 0 *'
           END IF
   END IF
END IF
IF ( DUMP( 4 ) .AND. DUMP( 8 ) .AND. RESULT( 4 ).GT.0 ) THEN
     WRITE( UNIT=FUDUMP( 4 ), FMT='(6X,A,I1,A,1P,5(A,E11.4))' ) &
            'TGK(', ITEST, ')', & 
            ': RESD=', RESULT( 4 ), ', ORTU=', RESULT( 5 ), &
            ', ORTV=', RESULT( 6 ), ', NRMU=', RESULT( 7 ), &
            ', NRMV=', RESULT( 8 )
END IF
!
! Print eigenvalues and eigenvectors in Matlab format.
!
STEST = IDTEST(ITEST)(5:6) // IDTEST(ITEST)(8:8)
!
IF ( DUMP( 6 ) .OR. DUMP( 7 ) .OR. DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP( 5 ), FMT='(''% '',A,1X,59(''=''))' ) IDTEST( ITEST )
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I5,A)' ) 'M =', M, ';'
END IF
IF ( DUMP( 6 ) .OR. DUMP( 8 ) ) THEN
   IF ( M == 0 ) THEN
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'W = [ ];'
   ELSE
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'W = zeros(M,1);'
      WRITE( UNIT=FUDUMP( 5 ), FMT='(2(A,I5,A,1P,E15.7E3,A))' ) &
             ( 'W(', I, ')=', W( I ), '; ', I = 1,M )
   END IF
   IF ( PRESENT( IINTO ) ) THEN
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A,I2.2,A)' ) &
             'W_', ICASE, '_', STEST, '_', IINTO, ' = W;'
   ELSE
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A)' ) &
             'W_', ICASE, '_', STEST, ' = W;'
   END IF
END IF
IF ( DUMP( 7 ) .OR. DUMP( 8 ) ) THEN
   IF ( M == 0 ) THEN
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'Z = [ ];'
   ELSE
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'Z = zeros(N,M);'
      WRITE( UNIT=FUDUMP( 5 ), FMT='(2(2(A,I5),A,1P,E15.7E3,A))' ) &
             ( ( 'Z(', J, ',', I, ')=', Z( J, I ), '; ', J = 1,N ), I = 1,M ) 
   END IF
   IF ( PRESENT( IINTO ) ) THEN
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A,I2.2,A)' ) &
             'Z_', ICASE, '_', STEST, '_', IINTO, ' = Z;'
   ELSE
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A)' ) &
             'Z_', ICASE, '_', STEST, ' = Z;'
   END IF
END IF
IF ( DUMP( 6 ) .OR. DUMP( 7 ) ) THEN
   IF ( PRESENT( IINTO ) ) THEN
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A,I2.2,A)' ) &
             'M_', ICASE, '_', STEST, '_', IINTO, ' = M;'
   ELSE
      WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A)' ) &
             'M_', ICASE, '_', STEST, ' = M;'
   END IF
END IF
IF ( DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) '% T = P*[ 0 B, B'' 0 ]*P'''
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A)' ) &
          'S_', ICASE, '_', STEST, ' = W( 1:N/2 ) ;'
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) &
          'V_temp = Z( [1:2:N],1:N/2 );', &
          'U_temp = Z( [2:2:N],1:N/2 );', &
          'for i=1:size(V_temp,2), V_temp(:,i)=V_temp(:,i)/norm(V_temp(:,i)); end', &
          'for i=1:size(U_temp,2), U_temp(:,i)=U_temp(:,i)/norm(U_temp(:,i)); end'
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A,I3.3,3A)' ) &
          'V_', ICASE, '_', STEST, ' = V_temp;', &
          'U_', ICASE, '_', STEST, ' = U_temp;'
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A,4(I3.3,3A))' ) &
          'B_', ICASE, '_', STEST, ' = U_', ICASE, '_', STEST, &
          '*diag(S_', ICASE, '_', STEST, ')*V_', ICASE, '_', STEST, ''';'
   WRITE( UNIT=FUDUMP( 5 ), FMT='(2(A,I3.3,3A))' ) &
          'S_', ICASE, '_', STEST, ' = ', 'abs(S_', ICASE, '_', STEST, ');', &
          'V_', ICASE, '_', STEST, ' = ', '-V_', ICASE, '_', STEST, ';'
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'clear i V_temp U_temp;'
END IF
IF ( DUMP( 6 ) .OR. DUMP( 7 ) .OR. DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'clear M Z W;'
END IF
!
END SUBROUTINE SSTEPRNRSLT
