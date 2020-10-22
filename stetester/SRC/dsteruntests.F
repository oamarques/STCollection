SUBROUTINE DSTERUNTESTS( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, &
                         MINTT, TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, &
                         NILIU, VLVU, NVLVU, VALS, MTRX, M )
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
USE DSTEINTERFACES2, ONLY : DBDPRNRSLT, DSTEBNDGAP, DSTEMATGEN, DSTEPRNRSLT 
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
!  DSTERUNTESTS runs the LAPACK symmetric tridiagonal eigensolver subroutines  !
!  and checks the results.                                                     !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  HBANDA   (input) INTEGER                                                    !
!           Sets the halfbandwidth of the symmetric matrix to be generated and !
!           then tridiagonalized, i.e. a matrix with max(1,N*(HBANDA/100))     !
!           subdiagonals is generated.                                         !
!                                                                              !
!  HBANDR   (input) INTEGER                                                    !
!           Sets the halfbandwidth of the matrices Z'*T*Z and Z'*Z used in the !
!           tests, i.e. max(1,N*(HBANDR/100)) subdiagonals are computed. If    !
!           0 then the tests are not performed.                                !
!                                                                              !
!  SETTO0   (input) INTEGER                                                    !
!           Specifies the number of entries that will be randomnly set to      !
!           0 in the test matrix; N*(SETTO0/100) entries will be set to 0.     !
!                                                                              !
!  MAXN     (input) INTEGER                                                    !
!           Maximum dimension of the matrix to be tested.                      !
!                                                                              !
!  NCASE    (input) INTEGER                                                    !
!           Number of matrices defined in M.                                   !
!                                                                              !
!  NRILIU   (input) INTEGER                                                    !
!           Number of pairs of random indices of the smallest and largest      !
!           eigenvalues to be returned, used only when RANGE='I'.              !
!                                                                              !
!  NRVLVU   (input) INTEGER                                                    !
!           Number of random intervals to be searched for eigenvalues,         !
!           used only when RANGE='V'.                                          !
!                                                                              !
!  MINTT    (input) REAL( KIND=PREC )                                          !
!           Minimum total timing for each subroutine call, i.e.                !
!           time_per_call * number_of_calls > MINTT                            !
!                                                                              !
!  TOLEBZ   (input) REAL( KIND=PREC )                                          !
!           Absolute error tolerance for the eigenvalues to be computed by     !
!           DSTEVX (see the leading comments of DSTEVX                         !
!           for details).                                                      !
!                                                                              !
!  TOLEGR   (input) REAL( KIND=PREC )                                          !
!           Absolute error tolerance for the eigenvalues to be computed by     !
!           DSTEGR (see the leading comments of DSTEGR for details).           !
!                                                                              !
!  CALLST   (input) LOGICAL array, dimension ( 8 )                             !
!           Flags the subroutines to be tested,                                !
!           CALLST( 1 ) : call DSTEQR with COMPZ='V'                           !
!           CALLST( 2 ) : call DSTEVX with RANGE='A'                           !
!           CALLST( 3 ) : call DSTEVX with RANGE='I'                           !
!           CALLST( 4 ) : call DSTEVX with RANGE='V'                           !
!           CALLST( 5 ) : call DSTEDC with COMPZ='I'                           !
!           CALLST( 6 ) : call DSTEGR with RANGE='A'                           !
!           CALLST( 7 ) : call DSTEGR with RANGE='I'                           !
!           CALLST( 8 ) : call DSTEGR with RANGE='V'                           !
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
!           DUMP( 8 ) : compute B = U*S*V' through DBDSQR and DBDSDC           !
!                       in TGK mode for comparison purposes                    !
!                                                                              !
!  SVD      (input) LOGICAL                                                    !
!           Computes SVD(B) in the tridiagonal derived from [ 0 B, B' 0 ].     ! 
!                                                                              !
!  TGK      (input) LOGICAL                                                    !
!           Triggers the generation of an augmented tridiagonal as a           !
!           permutation of [ 0 B, B' 0 ] where B is a bidiagonal               !
!           representation of the input T. When set, it applies                !
!           to all test matrices.                                              !
!                                                                              !
!  ILIU     (input) I_LIST (derived data type)                                 !
!           Indices (in ascending order) of the smallest and largest           !
!           eigenvalues to be returned, used only when RANGE='I'.              !
!                                                                              !
!  NILIU    (input) INTEGER                                                    !
!           Number of pairs of indices in ILIU.                                !
!                                                                              !
!  VLVU     (input) V_LIST (derived data type)                                 !
!           Lower and upper bounds of intervals to be searched for             !
!           eigenvalues, used only when RANGE='V'.                             !
!                                                                              !
!  NVLVU    (input) INTEGER                                                    !
!           Number of intervals in VLVU.                                       !
!                                                                              !
!  VALS     (input) VALS_LIST (derived data type)                              !
!           Values read from files.                                            !
!                                                                              !
!  MTRX     (input) MTRX_LIST (derived data type)                              !
!           Matrices read from files.                                          !
!                                                                              !
!  M        (input) M_LIST (derived data type)                                 !
!           Properties of the matrices to be used in the tests.                !
!                                                                              !
!==============================================================================!
!
!.. Local Scalars ..
LOGICAL :: NOVLVU, REPEAT, TESTILIU, TESTVLVU, TGKMODE
INTEGER :: IBLOCK, ICASE, ID, IE, IFAIL, IINT, IIWORK, IL, INFO, IRWORK, &
           IRWRKU, IRWRKV, ISPLIT, ISUPPZ, IU, IW, IWEBZ, J, K, LDZ, &
           LIWORK, LRWORK, N, NBD, NEIG, NRUN, NRUNMX = 123, NS
REAL( KIND=PREC ) :: FLOPS, GAPMIN, TIME0, TIME1, TOL, TTIME, VL, VU
!
!.. Static Arrays ..
LOGICAL :: CHKRSLT( 3 )
INTEGER :: ISEED( 4 )
REAL( KIND=PREC ) :: RESULT( 8 )
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: IWORK( : ), ILIU_A( :,: )
REAL( KIND=PREC ), ALLOCATABLE :: D( : ), E( : ), RWORK( : ), &
                                  VLVU_A( :,: ), Z( :,: )
!
!.. Local Derived Data Types ..
TYPE( VALS_LIST ), POINTER :: FIRST_VALS
TYPE( MTRX_LIST ), POINTER :: FIRST_MTRX
TYPE( I_LIST ), POINTER :: FIRST_ILIU 
TYPE( V_LIST ), POINTER :: FIRST_VLVU
TYPE( M_LIST ), POINTER :: FIRST_M
!
!.. External Functions ..
REAL( KIND=PREC ), EXTERNAL :: DLAMCH, DSTEFLOPS, DSTETIMER
!
!.. External Subroutines ..
EXTERNAL DSTEDC, DSTEGR, DSTEQR, DSTEVX
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
WRITE( UNIT=FUOUT, FMT='(A,A)' ) '========== Running test cases ============'
!
IF ( DUMP( 1 ) ) OPEN( UNIT=FUDUMP( 1 ), FILE=FDUMP( 1 ) )
IF ( DUMP( 2 ) ) OPEN( UNIT=FUDUMP( 2 ), FILE=FDUMP( 2 ) )
IF ( DUMP( 3 ) ) OPEN( UNIT=FUDUMP( 3 ), FILE=FDUMP( 3 ) )
IF ( DUMP( 4 ) ) OPEN( UNIT=FUDUMP( 4 ), FILE=FDUMP( 4 ) )
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) .OR. DUMP( 8 ) ) THEN
   OPEN( UNIT=FUDUMP( 5 ), FILE=FDUMP( 5 ) )
END IF
!
ALLOCATE( FIRST_VALS ); NULLIFY( FIRST_VALS%NEXT ); FIRST_VALS => VALS
ALLOCATE( FIRST_MTRX ); NULLIFY( FIRST_MTRX%NEXT ); FIRST_MTRX => MTRX
ALLOCATE( FIRST_ILIU ); NULLIFY( FIRST_ILIU%NEXT ); FIRST_ILIU => ILIU
ALLOCATE( FIRST_VLVU ); NULLIFY( FIRST_VLVU%NEXT ); FIRST_VLVU => VLVU
ALLOCATE( FIRST_M ); NULLIFY( FIRST_M%NEXT ); FIRST_M => M
ALLOCATE( ILIU_A( 2,MAX(1,NILIU+NRILIU+NVLVU+NRVLVU) ) )
ALLOCATE( VLVU_A( 2,MAX(1,NILIU+NRILIU+NVLVU+NRVLVU) ) )
!
N = 0
IL = 0; IU = 0
VL = ZERO; VU = ZERO
TESTILIU = CALLST( 3 ) .OR. CALLST( 7 ) .OR. SVD
TESTVLVU = CALLST( 4 ) .OR. CALLST( 8 ) .OR. SVD
CHKRSLT( 1 ) = .TRUE.; CHKRSLT( 2:3 ) = TGK
TGKMODE = SVD .OR. TGK
ISEED = ISEED_INIT
FLOPS = ZERO
!
DO IINT = 1, NILIU
   ILIU_A( 1,IINT ) = ILIU%IL
   ILIU_A( 2,IINT ) = ILIU%IU
   ILIU => ILIU%NEXT
END DO
DO IINT = 1, NVLVU
   VLVU_A( 1,IINT ) = VLVU%VL
   VLVU_A( 2,IINT ) = VLVU%VU
   VLVU => VLVU%NEXT
END DO
!
LIWORK = MAX( 3 + 5*MAXN, & !* Work space required by DSTEDC *!
              10*MAXN     ) !* Work space required by DSTEGR *!
LRWORK = MAX( 1 + 4*MAXN + MAXN**2, & !* Work space required by DSTEDC *!
              18*MAXN,              & !* Work space required by DSTEGR *!
              2*MAXN-2,             & !* Work space required by DSTEQR *!
              5*MAXN,               & !* Work space required by DSTEVX *!      
              (MAXN+1)*MAXN         ) !* Work space required by DSTECHKRSLT *! 
!
! SVD mode:
! DBDSQR      needs U(NBD,NBD), VT(NBD,NBD), WORK( 4*NBD ) 
! DBDSDC      needs U(NBD,NBD), VT(NBD,NBD), WORK( 3*NBD**2+4*NBD )
! DBDSVDX     needs WORK( MAXN*(NBD+1)+4*NBD )
! DSTECHKRSLT needs WORK( MAXN*(NBD+1)*2 ), which includes the space 
!             needed in DBDSVDX
!
NBD = MAXN/2
LIWORK = LIWORK + NBD*20 + NRVLVU*2
LRWORK = LRWORK + MAX( NBD*NBD*5 + 4*NBD, MAXN*(NBD+1)*2+4*NBD )
!
LRWORK = LRWORK + MAXN*4
!
ALLOCATE( D( MAXN ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'D (in DSTERUNTESTS)' )
ALLOCATE( E( MAXN ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'E (in DSTERUNTESTS)' ) 
ALLOCATE( IWORK( LIWORK ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'IWORK (in DSTERUNTESTS)' ) 
ALLOCATE( RWORK( LRWORK ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'RWORK (in DSTERUNTESTS)' ) 
ALLOCATE( Z( MAXN,MAXN ), STAT=INFO )
IF ( INFO /= 0 ) CALL HANDLER( 2, 'Z (in DSTERUNTESTS)' ) 
!
NOVLVU = MAX( NVLVU, NRVLVU ).EQ.0
!
! If PAPI is used, initialize it.
!
TIME0 = DSTETIMER( 0 )
CALL DSTEPRNFLOP( 0, FLOPS, 0, 0 )
!
! Run test cases.
!
DO ICASE = 1, NCASE
!
   INFO = 0
   RESULT = 0
!
!  Generate the test matrix.
!
   CALL DSTEMATGEN( N, D, E, HBANDA, SETTO0, ISEED, RWORK, LRWORK, &
                    VALS, MTRX, M, ICASE, DUMP, TGKMODE )
!
   LDZ = N
!
!  Set pointers for the real and integer work arrays.
!
   ID = 1 
   IE = ID + N 
   IW = IE + N 
   IWEBZ = IW + N 
   IRWORK = IWEBZ + N 
!
   IFAIL = 1
   IBLOCK = 1
   ISPLIT = IBLOCK + N
   ISUPPZ = 1
   IIWORK = ISUPPZ + N*2
!
!  Compute pairs IL,IU and VL,VU for cases RANGE="V"/"I".
!
   IF ( TESTILIU .OR. TESTVLVU ) THEN
!
!     Compute the eigenvalues of the matrix using bisection, by calling
!     DSTEVX with JOBZ='N'. This information will be used later to 
!     compare the results. (DSTEVX scaling of the matrix may be
!     important in some cases.)
!
      TOL = 2*DLAMCH('Safe minimum')
      CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
      CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
      CALL DSTEVX( 'N', 'A', N, RWORK( ID ), RWORK( IE ), ZERO, ZERO, 0, 0, &
                   TOLEBZ, NEIG, RWORK( IWEBZ ), Z, LDZ, RWORK( IRWORK ), &
                   IWORK( IIWORK ), IWORK( IFAIL ), INFO)
!
!     Choose random values for IL and IU. In SVD or TGK mode the 
!     values must be in [1,N/2].
!
      IF ( TGKMODE ) THEN
         K = N/2
         J = IWEBZ + K
      ELSE
         K = N
         J = IWEBZ
      END IF
!
      CALL DSTERNDILIU( K, ISEED, NRILIU, ILIU_A( 1,NILIU+1 ) )
!
!     Choose values for VL and VU. If NVLVU=NRVLVU=0 then the values of VL 
!     and VU are computed based on the NRILIU values of IL and IU obtained 
!     above, otherwise the values of VL and VU are chosen randomly. 
!
      IF ( NOVLVU ) THEN
         CALL DSTERNDVLVU( K, RWORK( J ), NRILIU, VLVU_A( 1,1 ), &
                           ILIU_A( 1,1 ) )
         NRVLVU = NRILIU
      ELSE
         CALL DSTERNDILIU( K, ISEED, NRVLVU, IWORK )
         CALL DSTERNDVLVU( K, RWORK( J ), NRVLVU, VLVU_A( 1,NVLVU+1 ), &
                           IWORK )
      END IF
!
      IF ( TGKMODE ) THEN
         DO IINT = 1, NVLVU+NRVLVU
            VLVU_A( 1,IINT ) = MAX( ZERO, VLVU_A( 1,IINT ) )
         END DO
      END IF
!
   END IF
!
!  Test DSTEQR( COMPZ='I' ) ....................................................
!
   IF ( CALLST( 1 ) ) THEN
      NRUN = 0
      TTIME = ZERO
      REPEAT = .TRUE.
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   ) 
         TIME0 = DSTETIMER( 1 )
         CALL DSTEQR( 'I', N, RWORK( ID ), RWORK( IE ), &
                      Z, LDZ, RWORK( IRWORK ), INFO )
         TIME1 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME1 < TIME0 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 TTIME = TTIME + ( TIME1 - TIME0 )
                 REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
         END IF
      END DO
      RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
      CALL DSTECHKRSLT( N, D, E, N, RWORK( ID ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
      CALL DSTEPRNRSLT( N, N, RWORK( ID: ), Z, LDZ, RESULT, &
                        INFO, 1, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 1, ICASE )
   END IF
!
!  Test DSTEVX( RANGE='A' ) ....................................................
!
   IF ( CALLST( 2 ) ) THEN
      NRUN = 0
      TTIME = ZERO
      REPEAT = .TRUE.
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )      
         TIME0 = DSTETIMER( 1 )
         CALL DSTEVX( 'V', 'A', N, RWORK( ID ), RWORK( IE ), VL, VU, IL, IU, &
                      TOLEBZ, NEIG, RWORK( IW ), Z, LDZ, RWORK( IRWORK ), &
                      IWORK( IIWORK ), IWORK( IFAIL ), INFO)
         TIME1 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME1 < TIME0 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 TTIME = TTIME + ( TIME1 - TIME0 )
                 REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
         END IF
      END DO
      RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
      CALL DSTECHKRSLT( N, D, E, NEIG, RWORK( IW ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
      CALL DSTEPRNRSLT( N, NEIG, RWORK( IW: ), Z, LDZ, RESULT, &
                        INFO, 2, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 2, ICASE )
   END IF
!
!  Test DSTEVX( RANGE='I' ) ....................................................
!
   IF ( CALLST( 3 ) ) THEN
      DO IINT = 1, NILIU+NRILIU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         IL = ILIU_A( 1,IINT )
         IU = ILIU_A( 2,IINT )         
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME0 = DSTETIMER( 1 ) 
            CALL DSTEVX( 'V', 'I', N, RWORK( ID ), RWORK( IE ), VL, VU, IL, IU, &
                         TOLEBZ, NEIG, RWORK( IW ), Z, LDZ, RWORK( IRWORK ), &
                         IWORK( IIWORK ), IWORK( IFAIL ), INFO)
            TIME1 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME1 < TIME0 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    TTIME = TTIME + ( TIME1 - TIME0 )
                    REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
            END IF
         END DO
         RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, ILO=IL, IUO=IU )
         CALL DSTECHKRSLT( N, D, E, NEIG, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
         CALL DSTEPRNRSLT( N, NEIG, RWORK( IW: ), Z, LDZ, RESULT, INFO, 3, &
                           ICASE, DUMP, IINTO=IINT, ILO=IL, IUO=IU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 3, ICASE )
      END DO
   END IF
!
!  Test DSTEVX( RANGE='V' ) ....................................................
!
   IF ( CALLST( 4 ) ) THEN
      DO IINT = 1, NVLVU+NRVLVU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         VL = VLVU_A( 1,IINT )
         VU = VLVU_A( 2,IINT )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )         
            TIME0 = DSTETIMER( 1 )
            CALL DSTEVX( 'V', 'V', N, RWORK( ID ), RWORK( IE ), VL, VU, IL, IU, &
                         TOLEBZ, NEIG, RWORK( IW ), Z, LDZ, RWORK( IRWORK ), &
                         IWORK( IIWORK ), IWORK( IFAIL ), INFO)
            TIME1 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME1 < TIME0 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    TTIME = TTIME + ( TIME1 - TIME0 )
                    REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
            END IF
         END DO
         RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, VLO=VL, VUO=VU )
         CALL DSTECHKRSLT( N, D, E, NEIG, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
         CALL DSTEPRNRSLT( N, NEIG, RWORK( IW: ), Z, LDZ, RESULT, INFO, 4, &
                           ICASE, DUMP, IINTO=IINT, VLO=VL, VUO=VU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 4, ICASE )
      END DO
   END IF
!
!  Test DSTEDC( COMPZ='I' ) ....................................................
!
   IF ( CALLST( 5 ) ) THEN
      NRUN = 0
      TTIME = ZERO
      REPEAT = .TRUE.
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )
         TIME0 = DSTETIMER( 1 )
         CALL DSTEDC( 'I', N, RWORK( ID ), RWORK( IE ), Z, LDZ, RWORK( IRWORK ), &
                      LRWORK-N*4, IWORK, LIWORK, INFO )
         TIME1 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME1 < TIME0 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 TTIME = TTIME + ( TIME1 - TIME0 )
                 REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
         END IF
      END DO
      RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
      CALL DSTECHKRSLT( N, D, E, N, RWORK( ID ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
      CALL DSTEPRNRSLT( N, N, RWORK( ID: ), Z, LDZ, RESULT, INFO, &
                        5, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 5, ICASE )
   END IF
!    
!  Test DSTEGR( RANGE='A' ) ....................................................
!
   IF ( CALLST( 6 ) ) THEN
      NRUN = 0       
      TTIME = ZERO
      REPEAT = .TRUE.
      DO WHILE ( REPEAT )
         CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
         CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
         FLOPS = DSTEFLOPS(   )
         TIME0 = DSTETIMER( 1 )   
         CALL DSTEGR( 'V', 'A', N, RWORK( ID ), RWORK( IE ), VL, VU, &
                      IL, IU, TOLEGR, NEIG, RWORK( IW ), Z, LDZ, IWORK( ISUPPZ ), &
                      RWORK( IRWORK ), LRWORK-N*4, IWORK( IIWORK ), &
                      LIWORK-N*2, INFO )
         TIME1 = DSTETIMER( 1 )
         FLOPS = DSTEFLOPS(   )
         IF      ( INFO /= 0 ) THEN
                 REPEAT = .FALSE.
         ELSE IF ( TIME1 < TIME0 ) THEN
                 REPEAT = .TRUE.
         ELSE
                 NRUN = NRUN + 1
                 TTIME = TTIME + ( TIME1 - TIME0 )
                 REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
         END IF
      END DO
      RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
      CALL DSTECHKRSLT( N, D, E, NEIG, RWORK( IW ), Z, LDZ, HBANDR, &
                        RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
      CALL DSTEPRNRSLT( N, NEIG, RWORK( IW: ), Z, LDZ, RESULT, INFO, &
                        6, ICASE, DUMP )
      CALL DSTEPRNFLOP( N, FLOPS, 6, ICASE )
   END IF
!
!  Test DSTEGR( RANGE='I' ) ....................................................
!
   IF ( CALLST( 7 ) ) THEN
      DO IINT = 1, NILIU+NRILIU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         IL = ILIU_A( 1,IINT )
         IU = ILIU_A( 2,IINT )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME0 = DSTETIMER( 1 )
            CALL DSTEGR( 'V', 'I', N, RWORK( ID ), RWORK( IE ), VL, VU, &
                         IL, IU, TOLEGR, NEIG, RWORK( IW ), Z, LDZ, IWORK( ISUPPZ ), &
                         RWORK( IRWORK ), LRWORK-N*4, IWORK( IIWORK ), &
                         LIWORK-N*2, INFO )
            TIME1 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME1 < TIME0 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    TTIME = TTIME + ( TIME1 - TIME0 )
                    REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
            END IF
         END DO
         RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, ILO=IL, IUO=IU )
         CALL DSTECHKRSLT( N, D, E, NEIG, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
         CALL DSTEPRNRSLT( N, NEIG, RWORK( IW: ), Z, LDZ, RESULT, INFO, 7, &
                           ICASE, DUMP, IINTO=IINT, ILO=IL, IUO=IU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 7, ICASE )
      END DO
   END IF
!
!  Test DSTEGR( RANGE='V' ) ....................................................
!
   IF ( CALLST( 8 ) ) THEN
      DO IINT = 1, NVLVU+NRVLVU
         NRUN = 0
         TTIME = ZERO
         REPEAT = .TRUE.
         VL = VLVU_A( 1,IINT )
         VU = VLVU_A( 2,IINT )
         DO WHILE ( REPEAT )
            CALL DCOPY( N, D, 1, RWORK( ID ), 1 )
            CALL DCOPY( N, E, 1, RWORK( IE ), 1 )
            FLOPS = DSTEFLOPS(   )
            TIME0 = DSTETIMER( 1 )
            CALL DSTEGR( 'V', 'V', N, RWORK( ID ), RWORK( IE ), VL, VU, &
                         IL, IU, TOLEGR, NEIG, RWORK( IW ), Z, LDZ, IWORK( ISUPPZ ), &
                         RWORK( IRWORK ), LRWORK-N*4, IWORK( IIWORK ), &
                         LIWORK-N*2, INFO )
            TIME1 = DSTETIMER( 1 )
            FLOPS = DSTEFLOPS(   )
            IF      ( INFO /= 0 ) THEN
                    REPEAT = .FALSE.
            ELSE IF ( TIME1 < TIME0 ) THEN
                    REPEAT = .TRUE.
            ELSE
                    NRUN = NRUN + 1
                    TTIME = TTIME + ( TIME1 - TIME0 )
                    REPEAT = ( TTIME>ZERO .AND. TTIME<MINTT .AND. NRUN<NRUNMX )
            END IF
         END DO
         RESULT( 1 ) = MAX( TTIME / MAX( 1, NRUN ), MICRO )
         CALL DSTEBNDGAP( N, RWORK( IWEBZ: ), GAPMIN, VLO=VL, VUO=VU )
         CALL DSTECHKRSLT( N, D, E, NEIG, RWORK( IW ), Z, LDZ, HBANDR, &
                           RWORK( IRWORK ), RESULT( 2 ), CHKRSLT, INFO )
         CALL DSTEPRNRSLT( N, NEIG, RWORK( IW: ), Z, LDZ, RESULT, INFO, 8, &
                           ICASE, DUMP, IINTO=IINT, VLO=VL, VUO=VU, &
                           GAPMINO=GAPMIN )
         CALL DSTEPRNFLOP( N, FLOPS, 8, ICASE )
      END DO
   END IF
!
!
END DO
! 
IF ( DUMP( 1 ) ) CLOSE( UNIT=FUDUMP( 1 ) )
IF ( DUMP( 2 ) ) CLOSE( UNIT=FUDUMP( 2 ) )
IF ( DUMP( 3 ) ) CLOSE( UNIT=FUDUMP( 3 ) )
IF ( DUMP( 4 ) ) CLOSE( UNIT=FUDUMP( 4 ) )
IF ( DUMP( 5 ) .OR. DUMP( 6 ) .OR. DUMP( 7 ) .OR. DUMP( 8 ) ) THEN
   WRITE( UNIT=FUDUMP( 5 ), FMT='(A)' ) 'clear D E N;'
   CLOSE( UNIT=FUDUMP( 5 ) )
END IF
CALL DSTEPRNFLOP( 0, FLOPS, 0, NCASE )
!
VALS => FIRST_VALS
MTRX => FIRST_MTRX
ILIU => FIRST_ILIU
VLVU => FIRST_VLVU
M => FIRST_M
!
DEALLOCATE( ILIU_A )
DEALLOCATE( VLVU_A )
DEALLOCATE( D )
DEALLOCATE( E )
DEALLOCATE( IWORK )
DEALLOCATE( RWORK )
DEALLOCATE( Z )
!
END SUBROUTINE DSTERUNTESTS
