SUBROUTINE SSTEGETDATA( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, &
                        MINTT, TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, &
                        ILIU, NILIU, VLVU, NVLVU, VALS, MTRX, M )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1
USE SSTEDEFINITIONS
USE SSTEINTERFACES1
USE SSTEINTERFACES2
! 
!.. Scalar Arguments ..
LOGICAL :: SVD, TGK
INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU, SETTO0
REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
!
!.. Array Arguments ..
LOGICAL :: CALLST( 8 ), DUMP( 8 )
!
!.. Derived Data Type Arguments ..
TYPE( VALS_LIST ), POINTER :: VALS
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( I_LIST ), POINTER :: ILIU
TYPE( V_LIST ), POINTER :: VLVU
TYPE( M_LIST ), POINTER :: M
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEGETDATA interprets the input file.                                      !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  HBANDA   (output) INTEGER                                                   !
!           Sets the halfbandwidth of the symmetric matrix to be generated and !
!           then tridiagonalized, i.e. a matrix with max(1,N*(HBANDA/100))     !
!           subdiagonals is generated.                                         !
!                                                                              !
!  HBANDR   (output) INTEGER                                                   !
!           Sets the halfbandwidth of the matrices Z'*T*Z and Z'*Z used in the !
!           tests, i.e. max(1,N*(HBANDR/100)) subdiagonals are computed. If    !
!           0 then the tests are not performed.                                !
!                                                                              !
!  SETTO0   (output) INTEGER                                                   !
!           Specifies the number of entries that will be randomnly set to      !
!           0 in the test matrix; N*(SETTO0/100) entries will be set to 0.     !
!                                                                              !
!  MAXN     (output) INTEGER                                                   !
!           Maximum dimension of the matrix to be tested.                      !
!                                                                              !
!  NCASE    (output) INTEGER                                                   !
!           Number of tridiagonal matrices defined in T.                       !
!                                                                              !
!  NRILIU   (output) INTEGER                                                   !
!           Number of random pair of indexes of the smallest and largest       !
!           eigenvalues to be computed.                                        !
!                                                                              !
!  NRVLVU   (output) INTEGER                                                   !
!           Number of random intervals to be searched for eigenvalues.         !
!                                                                              !
!  MINTT    (output) REAL( KIND=PREC )                                         !
!           Minimum total timing for each subroutine call, i.e.                !
!           time_per_call * number_of_calls > MINTT                            !
!                                                                              !
!  TOLEBZ   (output) REAL( KIND=PREC )                                         !
!           Absolute error tolerance for the eigenvalues to be computed by     !
!           SSTEBZ (see the leading comments of SSTEBZ for details).           !
!                                                                              !
!  TOLEGR   (output) REAL( KIND=PREC )                                         !
!           Absolute error tolerance for the eigenvalues to be computed by     !
!           SSTEGR (see the leading comments of SSTEGR for details).           !
!                                                                              !
!  DUMP     (output) LOGICAL array, dimension ( 8 )                            !
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
!  CALLST   (output) LOGICAL array, dimension ( 8 )                            !
!           Flags the subroutines to be tested,                                !
!           CALLST( 1 ) : call SSTEQR with COMPZ='V'                           !
!           CALLST( 2 ) : call SSTEBZ with RANGE='A'                           !
!           CALLST( 3 ) : call SSTEBZ with RANGE='I'                           !
!           CALLST( 4 ) : call SSTEBZ with RANGE='V'                           !
!           CALLST( 5 ) : call SSTEDC with COMPZ='I'                           !
!           CALLST( 6 ) : call SSTEGR with RANGE='A'                           !
!           CALLST( 7 ) : call SSTEGR with RANGE='I'                           !
!           CALLST( 8 ) : call SSTEGR with RANGE='V'                           !
!                                                                              !
!  SVD      (output) LOGICAL                                                   !
!           Computes SVD(B) in the tridiagonal derived from [ 0 B, B' 0 ].     ! 
!                                                                              !
!  TGK      (output) LOGICAL                                                   !
!           Activates the generation of an augmented tridiagonal as a          !
!           permutation of [ 0 B', B 0 ] where B is a bidiagonal               !
!           matrix. When set, it applies to all test matrices.                 !
!                                                                              !
!  ISEED    (input/output) INTEGER array, dimension ( 4 )                      !
!           Initial seed of the random number generator. Each entry of ISEED   !
!           should lie between 0 and 4095 and ISEED(4) should be odd.          !
!                                                                              !
!  ILIU     (output) I_LIST (derived data type)                                !
!           Indices (in ascending order) of the smallest and largest           !
!           eigenvalues to be computed, used only when RANGE='I'.              !
!                                                                              !
!  NILIU    (output) INTEGER                                                   !
!           Number of pair of indices stored in ILIU.                          !
!                                                                              !
!  VLVU     (output) V_LIST (derived data type)                                !
!           Lower and upper bounds of intervals to be searched for             !
!           eigenvalues, used only when RANGE='V'.                             !
!                                                                              !
!  NVLVU    (output) INTEGER                                                   !
!           Number of intervals stored in VLVU.                                !
!                                                                              !
!  VALS     (output) VALS_LIST (derived data type)                             !
!           Values read from files.                                            !
!                                                                              !
!  MTRX     (output) MTRX_LIST (derived data type)                             !
!           Tridiagonal matrices read from files.                              !
!                                                                              !
!  M        (output) M_LIST (derived data type)                                !
!           Properties of the matrices to be used in the tests.                !
!                                                                              !
!==============================================================================!
!                                                                              !
!  Notes:                                                                      !
!                                                                              !
!  ECOND is the condition number for some built-in eigenvalue distributions    !
!        = 1, COND = 1 / SQRT(ULP), default,                                   !
!        = 2, COND = 1 / (N*SQRT(ULP)),                                        !
!        = 3, COND = 1 / (10*N*SQRT(ULP)),                                     !
!        = 4, COND = 1 / ULP,                                                  !
!        = 5, COND = 1 / (N*ULP),                                              !
!        = 6, COND = 1 / (10*N*ULP),                                           !
!                                                                              !
!  EDIST sets the type of the distribution to be used for random eigenvaules   !
!        = 1, uniform distribution (-1,1), default                             !
!        = 2, uniform distribution ( 0,1)                                      !
!        = 3, normal distribution  ( 0,1)                                      !
!                                                                              !
!  ESIGN assigns signs to the eigenvalues                                      !
!        = 0, the eigenvalues will not be negative, default                    !
!        = 1, the eigenvalues may be positive, negative, or zero               !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO 
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: ECOND = 1, EDIST = 1, ESIGN = 0, NTEST = 0
!
!.. Local Static Arrays ..
INTEGER :: ISEED(4), ITEMP( 1 )
REAL( KIND=PREC ) :: RTEMP( 1 )
! 
!.. Local Derived Data Types ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
TYPE( VALS_LIST ), POINTER :: FIRST_VALS
TYPE( MTRX_LIST ), POINTER :: FIRST_MTRX
TYPE( I_LIST ), POINTER :: FIRST_ILIU
TYPE( V_LIST ), POINTER :: FIRST_VLVU
TYPE( M_LIST ), POINTER :: FIRST_M
!
!.. External Subroutine ..
EXTERNAL HANDLER 
!
!.. External Function ..
CHARACTER( LEN=RECORD_LENGTH ), EXTERNAL :: GETRECORD
!
!.. Executable statements ......................................................
!
ISEED = ISEED_INIT
!
WRITE( UNIT=FUOUT, FMT='(A,A)' ) '========== Reading data =================='
!
ALLOCATE( FIRST_VALS ); NULLIFY( FIRST_VALS%NEXT ); FIRST_VALS => VALS
ALLOCATE( FIRST_MTRX ); NULLIFY( FIRST_MTRX%NEXT ); FIRST_MTRX => MTRX
ALLOCATE( FIRST_ILIU ); NULLIFY( FIRST_ILIU%NEXT ); FIRST_ILIU => ILIU
ALLOCATE( FIRST_VLVU ); NULLIFY( FIRST_VLVU%NEXT ); FIRST_VLVU => VLVU
ALLOCATE( FIRST_M ); NULLIFY( FIRST_M%NEXT ); FIRST_M => M
!
! Parse input file.
!
INPUT: DO
!
   RECORD = GETRECORD( )
   LIST => PARSER( RECORD )
   MACRO = LIST%FIELD
   LIST => LIST%NEXT
   WRITE( UNIT=FUOUT, FMT='(A,A)' ) '** Interpreting key word: ', MACRO
!
   SELECT CASE ( MACRO )
!
   CASE ( 'CALLST', 'callst' ) 
!
!       Read list of subroutines to be tested.
!
        CALL SSTEMCCALLST( NTEST, MACRO, LIST, CALLST )
!
   CASE ( 'DUMP', 'dump' ) 
!
!       Set variables to be saved.
!
        CALL SSTEMCDUMP( MACRO, LIST, DUMP )

   CASE ( 'ECOND', 'econd' ) 
!
!       Set condition number.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); ECOND = ITEMP( 1 )
        IF ( ECOND<1 .OR. ECOND>6 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'EDIST', 'edist' ) 
!
!       Set distribution for random eigenvalues.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); EDIST = ITEMP( 1 )
        IF ( EDIST<1 .OR. EDIST>3 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'ESIGN', 'esign' ) 
!
!       Set eigenvalue signs.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); ESIGN = ITEMP( 1 )
        IF ( ESIGN<0 .OR. ESIGN>1 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'GLUED', 'glued' ) 
!
!       Glued tridiagonal matrices.
!
        CALL SSTEMCGLUED( ECOND, EDIST, ESIGN, ISEED, MAXN, NCASE, MACRO, M )
!
   CASE ( 'HBANDA', 'hbanda' )
!
!       Set halfbandwidth of the symmetric matrix to be generated.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); HBANDA = ITEMP( 1 )
        IF ( HBANDA<1 .OR. HBANDA>100 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'HBANDR', 'hbandr' )
!
!       Set halfbandwidth of the the matrices Z'*T*Z and Z'*Z.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); HBANDR = ITEMP( 1 )
        IF ( HBANDR<0 .OR. HBANDR>100 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'ISEED', 'iseed' ) 
!
!       Set seed for random number generator.
!
        ISEED = GETINTGR( MACRO, LIST, 4 )
        IF ( MINVAL( ISEED ) < 1 ) CALL HANDLER( 1, MACRO )

   CASE ( 'MATRIX', 'matrix' ) 
!
!       Built-in tridiagonal matrices.
!
        CALL SSTEMCMATRIX( MAXN, NCASE, MACRO, LIST, M )
!
   CASE ( 'MATRIXF', 'matrixf' ) 
!
!       Read tridiagonal matrix stored in file.
!
        CALL SSTEMCMATRIXF( MAXN, NCASE, MACRO, LIST, MTRX, M )
!
   CASE ( 'MINTT', 'mintt' ) 
!
!       Set minimum total timing for each subroutine call
!
        RTEMP = GETSREAL( MACRO, LIST, 1 ); MINTT = RTEMP( 1 )
!
   CASE ( 'NRILIU', 'nriliu' )
!
!       Set number of random intervals.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); NRILIU = ITEMP( 1 )
        IF ( NRILIU<0 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'NRVLVU', 'nrvlvu' )
!
!       Set number of random intervals.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); NRVLVU = ITEMP( 1 )
        IF ( NRVLVU<0 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'SETTO0', 'setto0' )
!
!       Entries of the matrix (%) that will be randomnly set to 0.
!
        ITEMP = GETINTGR( MACRO, LIST, 1 ); SETTO0 = ITEMP( 1 )
        IF ( SETTO0<0 .OR. SETTO0>100 ) CALL HANDLER( 1, MACRO )
!
   CASE ( 'SNGVAL', 'sngval' )
!
!       Built-in singular value distributions.
!
        CALL SSTEMCVALUES( ECOND, EDIST, ESIGN, ISEED, MAXN, NCASE, MACRO, 5, M)
!
   CASE ( 'SNGVALF', 'sngvalf' )
!
!       Read singular value distribution stored in file.
!
        CALL SSTEMCVALUESF( ISEED, MAXN, NCASE, MACRO, LIST, VALS, 7, M )
!
   CASE ( 'SVD', 'svd' ) 
!
!       Computes SVD(B) in the tridiagonal derived from [ 0 B, B' 0 ]. 
!       Used only if TGK = .TRUE. (see below)
!
        SVD = .TRUE.
!
   CASE ( 'TGK', 'tgk' ) 
!
!       Activates the augmented tridiagonal mode.
!
        TGK = .TRUE.
!
   CASE ( 'TOLEBZ', 'tolebz' )
!
!       Set absolute error for SSTEBZ.
!
        RTEMP = GETSREAL( MACRO, LIST, 1 ); TOLEBZ = RTEMP( 1 )
!
   CASE ( 'TOLEGR', 'tolegr' )
!
!       Set absolute error for SSTEGR.
!
        RTEMP = GETSREAL( MACRO, LIST, 1 ); TOLEGR = RTEMP( 1 )
!
   CASE ( 'VALSI', 'valsi' ) 
!
!       Read indices of the smallest and largest eigenvalues to be computed.
!
        CALL SSTEMCVALSI( MACRO, ILIU, NILIU )
!
   CASE ( 'VALSV', 'valsv' ) 
!
!       Read lower and upper bounds of intervals to be searched for eigenvalues.
!
        CALL SSTEMCVALSV( MACRO, VLVU, NVLVU )
!
   CASE ( 'VALUES', 'values' ) 
!
!       Built-in eigenvalue distributions.
!
        CALL SSTEMCVALUES( ECOND, EDIST, ESIGN, ISEED, MAXN, NCASE, MACRO, 1, M )
!
   CASE ( 'VALUESF', 'valuesf' ) 
!
!       Read eigenvalue distribution stored in file.
!
        CALL SSTEMCVALUESF( ISEED, MAXN, NCASE, MACRO, LIST, VALS, 3, M )
!
   CASE ( 'END', 'end' ) 
!
!       End of input file.
!
        EXIT INPUT
!
   CASE DEFAULT
!
!       Macro not implemented.
!
        CALL HANDLER( 5, MACRO )
!
   END SELECT
!
END DO INPUT
!
IF ( NCASE == 0 ) CALL HANDLER( 7, 'NCASE = 0 !' )
IF ( NTEST == 0 .AND. TGK ) CALL HANDLER( 7, 'NTEST = 0 and TGK mode on!' )
IF ( NTEST == 0 .AND. .NOT.SVD ) CALL HANDLER( 7, 'NTEST = 0 and SVD mode off !' )
!
VALS => FIRST_VALS
MTRX => FIRST_MTRX
ILIU => FIRST_ILIU
VLVU => FIRST_VLVU
M => FIRST_M
!
IF ( SVD ) WRITE( UNIT=FUOUT, FMT='(A)' ) '** SSTETESTER warning: SVD mode on'
IF ( TGK ) WRITE( UNIT=FUOUT, FMT='(A)' ) '** SSTETESTER warning: TGK mode on'
IF ( SVD .OR. TGK ) MAXN = MAXN*2
!
END SUBROUTINE SSTEGETDATA
