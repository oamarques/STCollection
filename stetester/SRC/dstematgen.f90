SUBROUTINE DSTEMATGEN( N, D, E, HBANDA, SETTO0, ISEED, WORK, LWORK, &
                       VALS, MTRX, M, ICASE, DUMP, TGK )
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
USE MAPDDATA
!
!.. Scalar Arguments ..
LOGICAL :: TGK
INTEGER :: HBANDA, ICASE, LWORK, N, SETTO0
!
!.. Array Arguments ..
LOGICAL :: DUMP( 8 )
INTEGER :: ISEED( 4 )
REAL( KIND=PREC ) :: D( * ), E( * ), WORK( * )
!
!.. Derived Data Type Argument ..
TYPE( VALS_LIST ), POINTER :: VALS
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( M_LIST ), POINTER :: M
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEMATGEN generates various kinds of symmetric tridiagonal and             !
!  bidiagonal matrices.                                                        !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  N        (output) INTEGER                                                   !
!           The dimension of the matrix.                                       !
!                                                                              !
!  D        (output) REAL( KIND=PREC ) array, dimension ( N )                  !
!           The N diagonal elements of the tridiagonal matrix.                 !
!                                                                              !
!  E        (output) REAL( KIND=PREC ) array, dimension ( N )                  !
!           The N-1 off-diagonal elements of the tridiagonal matrix in         !
!           elements 1 to N-1, E(N) is set to zero.                            !
!                                                                              !
!  HBANDA   (input) INTEGER                                                    !
!           Sets the halfbandwidth of the symmetric matrix to be generated and !
!           then tridiagonalized, i.e. a matrix with max(1,N*(HBANDA/100))     !
!           subdiagonals is generated.                                         !
!                                                                              !
!  SETTO0   (input) INTEGER                                                    !
!           Specifies the number of entries of the generated matrix that       !
!           will be randomnly set to 0, i.e. (N*2-1)*(SETTO0/100) entries      !
!           will be set to 0.                                                  !
!                                                                              !
!  ISEED    (input/output) INTEGER array, dimension ( 4 )                      !
!           Seed for the random number generator. Each entry of ISEED should   !
!           lie between 0 and 4095 inclusive and ISEED(4) should be odd.       !
!                                                                              !
!  WORK     (workspace) REAL( KIND=PREC ) array, dimension ( LWORK )           !
!           Workspace.                                                         !
!                                                                              !
!  LWORK    (input) INTEGER                                                    !
!           Dimension of WORK.                                                 !
!                                                                              !
!  VALS     (input) VALS_LIST (derived data type)                              !
!           Values read from files.                                            !
!                                                                              !
!  MTRX     (input) MTRX_LIST (derived data type)                              !
!           Tridiagonal matrices read from files.                              !
!                                                                              !
!  M        (input) M_LIST (derived data type)                                 !
!           Properties of the matrices to be used in the tests.                !
!                                                                              !
!  ICASE    (input) INTEGER                                                    !
!           Case number.                                                       !
!                                                                              !
!  DUMP     (input) LOGICAL, dimension ( 8 )                                   !
!           Defines data to be written into files,                             !
!           DUMP( 1 ) : tridiagonal matrix (i,d_i,e_i)                         !
!           DUMP( 2 ) : eigenvalues                                            !
!           DUMP( 3 ) : eigenvectors                                           !
!           DUMP( 4 ) : timing, residuals, orthogonality                       !
!           DUMP( 5 ) : tridiagonal matrix (i,d_i,e_i) in Matlab format        !
!           DUMP( 6 ) : eigenvalues in Matlab format                           !
!           DUMP( 7 ) : eigenvectors in Matlab format                          !
!           DUMP( 8 ) : singular values and vectors in Matlab format           !
!                                                                              !
!  TGK      (input) LOGICAL                                                    !
!           Activates the generation of an augmented tridiagonal as a          !
!           permutation of [ 0 B', B 0 ] where B is a bidiagonal               !
!           matrix. When set, it applies to all test matrices.                 !
!                                                                              !
!==============================================================================!
!  
!.. Local Scalars ..
LOGICAL :: BDIAG, TDIAG
CHARACTER( LEN=1 ) :: PACK
INTEGER :: ECOND, EDIST, ESIGN, IA, IBA, IBB, INFO, IW, IWORK, &
           IZE, IZQ, J, K, MFORM, MSIZE, MTYPE, NDIAG
REAL( KIND=PREC ) :: GAMMA
!
!.. External Subroutines ..
EXTERNAL HANDLER, DLATMS, DSBTRD, DSYTRD
!
!.. External Function ..
REAL( KIND=PREC ), EXTERNAL :: DLARND
!
!.. Executable Statements ......................................................
!
WRITE( UNIT=FUOUT, FMT='(A4,I4)' ) 'Case', ICASE
!
IF ( DUMP( 1 ) ) WRITE( UNIT=FUDUMP( 1 ), FMT='(A,I4)' ) 'Case', ICASE
IF ( DUMP( 2 ) ) WRITE( UNIT=FUDUMP( 2 ), FMT='(A,I4)' ) 'Case', ICASE
IF ( DUMP( 3 ) ) WRITE( UNIT=FUDUMP( 3 ), FMT='(A,I4)' ) 'Case', ICASE
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) .OR. DUMP( 8 ) ) &
   WRITE( UNIT=FUDUMP( 5 ), FMT='(''% '',A,I5,1X,64(''#''))' ) 'Case', ICASE
!
N = 0
BDIAG = .FALSE.
TDIAG = .FALSE.
!
DO
!
!  Get data from M while GAMMA != 0
!
!  form    type    size      en   id
!  -------------------------------------------------------------------------
!    1    ETYPE   ESIZE   GAMMA   built-in eigenvalue distribution
!    2    MTYPE   MSIZE   GAMMA   built-in tridiagonal matrix
!    3        1    NEIG       0   eigenvalue distribution read from file
!    4        1     NDE       0   tridiagonal matrix read from file
!    5    STYPE   SSIZE   GAMMA   built-in singular value distribution
!    6    MTYPE   MSIZE   GAMMA   built-in bidiagonal matrix
!    7        1    NSVL       0   singular value distribution read from file
!    8        1     NDE       0   bidiagonal matrix read from file
!    9        -       -       -   reserved for future use
!  -------------------------------------------------------------------------
!
   MFORM = M%DATA%FORM
   MTYPE = M%DATA%TYPE
   MSIZE = M%DATA%SIZE
   ECOND = M%DATA%COND
   EDIST = M%DATA%DIST
   ESIGN = M%DATA%SIGN
   ISEED = M%DATA%SEED
   GAMMA = M%DATA%EN
!
   IF ( GAMMA == ZERO ) THEN
      WRITE( UNIT=FUOUT, FMT='(3X,A,I2,A,I2,A,I5)' ) &
             'Matrix data: form=', MFORM, ', type=', MTYPE, &
             ', size=', MSIZE
   ELSE
      WRITE( UNIT=FUOUT, FMT='(3X,A,I2,A,I2,A,I5,A,1P,E11.4)' ) &
             'Matrix data: form=', MFORM, ', type=', MTYPE, &
             ', size=', MSIZE, ', glue=', GAMMA 
   END IF
   IF ( DUMP( 1 ) ) THEN
      IF ( GAMMA == ZERO ) THEN
         WRITE( UNIT=FUDUMP(1), &
                FMT='(3X,A,I2,A,I2,A,I5,/,17X,A,4I5)' ) &
                'Matrix data:  form =', MFORM, ', type =', MTYPE, &
                ', size =', MSIZE, 'SEED =', ISEED
      ELSE
         WRITE( UNIT=FUDUMP(1), &
                FMT='(3X,A,I2,A,I2,A,I5,A,1P,E11.4,/,17X,A,4I5)' ) &
                'Matrix data:  form =', MFORM, ', type =', MTYPE, &
                ', size =', MSIZE, ', glue =', GAMMA, &
                'SEED =', ISEED
      END IF
   END IF
!
!  Generate matrix or read it from file.
!
   SELECT CASE ( MFORM )
!
   CASE ( 1, 5 )
!
!       Built-in eigenvalue / singular value distribution.
!
        CALL DSTEDEFVALS( ECOND, EDIST, ESIGN, MTYPE, ISEED, MSIZE, WORK( 1 ) )
!
   CASE ( 2 )
!
!       Built-in tridiagonal matrix.
!
        CALL DSTEDEFTMTRX( MTYPE, MSIZE, D( N+1 ), E( N+1 ) )
!
   CASE ( 3, 7 )
!
!       Eigenvalue / singular value distribution read from file.
!
        WORK( 1:MSIZE ) = VALS%S( 1:MSIZE )
        VALS => VALS%NEXT
!
   CASE ( 4, 8 )
!
!       Tridiagonal matrix read from file.
!
        D( 1:MSIZE ) = MTRX%D( 1:MSIZE )
        E( 1:MSIZE ) = MTRX%E( 1:MSIZE )
        MTRX => MTRX%NEXT
!
   CASE ( 6 )
!
!       Built-in bidiagonal matrix.
!
        CALL DSTEDEFBMTRX( MTYPE, MSIZE, D( N+1 ), E( N+1 ), ISEED, WORK )
!
   CASE DEFAULT
!
        CALL HANDLER( 7, 'Matrix form not implemented!' )
!
   END SELECT
!
   IF ( MFORM.GE.1 .AND. MFORM.LE.4 ) TDIAG = .TRUE.
   IF ( MFORM.GE.5 .AND. MFORM.LE.8 ) BDIAG = .TRUE.
!
!  Tridiagonal or bidiagonal form.

   IW = 1
   IA = IW + MSIZE
   IWORK = IA + MSIZE*MSIZE
!
   IF      ( MFORM == 1 .OR. MFORM == 3 ) THEN
!
!          Generate symmetric matrix, halfbandwidth max(1,MSIZE*(HBANDA/100)).
!
           IF ( HBANDA > 50 ) THEN
              PACK = 'N'  ! No packing.
           ELSE
              PACK = 'B'  ! Store the lower triangle.
           END IF
           NDIAG = MAX( 1, INT( MSIZE*(HBANDA/HNDRD) ) )
           CALL DLATMS( MSIZE, MSIZE, DIST( EDIST ), ISEED, SYMM( ESIGN ), &
                        WORK( IW ), 0, ZERO, ONE, NDIAG, NDIAG, PACK, &
                        WORK( IA ), MSIZE, WORK( IWORK ), INFO )  
           IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'DLATMS' )
!
!          Tridiagonalize matrix.
!
           IF ( PACK == 'N' ) THEN
              CALL DSYTRD( 'L', MSIZE, WORK( IA ), MSIZE, D( N+1 ), E( N+1 ), &
                           WORK( IW ), WORK( IWORK ), LWORK-IWORK+1, INFO )
              IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'DSYTRD' )
           ELSE
              CALL DSBTRD( 'N', 'L', MSIZE, NDIAG, WORK( IA ), MSIZE, D( N+1 ), &
                           E( N+1 ), WORK( IW ), 1, WORK( IW ), INFO )
              IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'DSBTRD' )
           END IF
!
   ELSE IF ( MFORM == 5 .OR. MFORM == 7 ) THEN
!
!          Generate a bidiagonal matrix B for TGK = [ 0 B', B 0 ]
!
           CALL DLATMS( MSIZE, MSIZE, DIST( EDIST ), ISEED, &
                       'N', WORK( IW ), 0, ZERO, ONE, 0, 1, 'N', &
                       WORK( IA ), MSIZE, WORK( IWORK ), INFO )  
           IF ( INFO /= 0 ) CALL HANDLER ( INFO, 'DLATMS' )
           CALL DCOPY( MSIZE, WORK( IA ), MSIZE+1, D( N+1), 1 )
           CALL DCOPY( MSIZE-1, WORK( IA+MSIZE ), MSIZE+1, E( N+1), 1 ) 
!
   END IF
!
   N = N + MSIZE
   E( N ) = GAMMA
   IF ( GAMMA == ZERO ) EXIT
   M => M%NEXT
!
END DO
!
IF ( SETTO0.GT.0 ) THEN
!
!  Set entries to zero randomnly.
!
   DO J = 1, MAX( 1, INT( (N*SETTO0)/HNDRD ) )
      K = 1 + INT( ( 2*N-1 )*DLARND( 1, ISEED ) )
      IF ( K > N ) THEN
         E( K-N ) = ZERO
      ELSE
         D( K ) = ZERO
      END IF
   END DO 
!
END IF
!
IZQ = 1
IZE = IZQ + N
IBA = IZE + N
IBB = IBA + N
!
IF      ( TGK ) THEN
!
!       Augmented tridiagonal mode: TGK = [ 0 B', B 0 ].
!   
        IF ( TDIAG ) THEN
!
!          B is obtained from T: 
!            1) Shift T if necessary; 
!            2) Map T into Z representation; 
!            3) Map Z into B.
!
           CALL DSQTSHIFT( N, D, E )
           CALL DT_TO_DZ( N, D, E, WORK( IZQ ), WORK( IZE ) )
           CALL DZ_TO_DB( N, WORK( IBA ), WORK( IBB ), WORK( IZQ ), WORK( IZE ) )
!
        ELSE
!
!          Matrix B is already bidiagonal.
!
           CALL DCOPY( N, D, 1, WORK( IBA ), 1 )
           CALL DCOPY( N, E, 1, WORK( IBB ), 1 )
!
        END IF
!
!       Copy the entries of B into D and E.
!
        D( 1:N*2 ) = ZERO
        DO J = 1, N-1
           E( 2*J-1 ) = WORK( IBA+J-1 )
           E( 2*J ) = WORK( IBB+J-1 )
        END DO
        E( 2*N-1 ) = WORK( IBA+N-1 )
        N = N*2
!
ELSE IF ( BDIAG ) THEN
!
!       T is obtained from B.
! 
        CALL DCOPY( N, D, 1, WORK( IBA ), 1 )
        CALL DCOPY( N, E, 1, WORK( IBB ), 1 )
        CALL DB_TO_DT( N, WORK( IBA ), WORK( IBB ), D, E )
!
END IF
!
E( N ) = ZERO
!
! Write matrix data into file.
!
IF ( DUMP( 1 ) ) THEN
   WRITE( UNIT=FUDUMP(1), FMT='(5X,A,6X,A,19X,A,/,(I6,1P,2E25.15))' ) &
          'J', 'D( J )', 'E( J )', ( J, D( J ), E( J ), J = 1,N )
END IF
IF ( DUMP( 4 ) ) THEN
   IF ( MFORM .EQ. 9 ) THEN
      WRITE( UNIT=FUDUMP( 4 ), FMT='(A,I4,A,I5,2X,5(A,I2),A)' ) &
             'Case', ICASE, ': N= ', N, '(form=', 9, ', type=', 9, &
             ', cond=', ECOND, ', dist=', EDIST, ', sign=', ESIGN, &
             ', glued matrix)'
   ELSE     
      WRITE( UNIT=FUDUMP( 4 ), FMT='(A,I4,A,I5,2X,5(A,I2),A)' ) &
             'Case', ICASE, ': N= ', N, '(form=', MFORM, ', type=', MTYPE, &
             ', cond=', ECOND, ', dist=', EDIST, ', sign=', ESIGN, ')'
   END IF
END IF
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) .OR. DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP(5), FMT='(A,I5,A)' ) 'N =', N, ';'
   WRITE( UNIT=FUDUMP(5), FMT='(A,I3.3,A)' ) 'N_', ICASE, ' = N;'
END IF
IF ( DUMP( 5 ) .OR. DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP(5), FMT='(A)' ) 'D = zeros(N,1); E = zeros(N,1);'
   WRITE( UNIT=FUDUMP(5), FMT='(2(A,I5,A,1P,E23.15E3,A))' ) &
          ( 'D(', J, ')=', D( J ), '; ', &
            'E(', J, ')=', E( J ), '; ', &
            J = 1, N )
   WRITE( UNIT=FUDUMP(5), FMT='(2(A,I3.3,A))' ) &
          'D_', ICASE, ' = D; ', 'E_', ICASE, ' = E;'
END IF
IF ( DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP(5), FMT='(A,I3.3,A)' ) &   
          'B_', ICASE, ' = diag(E(1:2:N)) + diag(E(2:2:N-1),1);', &
          'NB_', ICASE, ' = N/2;'
END IF
!
M => M%NEXT
!
END SUBROUTINE DSTEMATGEN
