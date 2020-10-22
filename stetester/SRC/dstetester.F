PROGRAM DSTETESTER
!
USE GSTEDEFINITIONS
USE DSTEDEFINITIONS
USE DSTEINTERFACES3
! 
!==============================================================================!
!                                                                              !
! Version 02/12/2015, by Osni Marques (LBNL)                                   !
!                                                                              !
! Purpose:                                                                     !
! =======                                                                      !
!                                                                              !
! DSTETESTER is a program to test the symmetric tridiagonal eigenvalue solvers !
! implemented in LAPACK.                                                       !
!                                                                              !
! The test matrices used by DSTETESTER can be generated automatically (as      !
! described below) or generated from data read from files. In this case, the   !
! the data input a symmetric tridiagonal matrix or an eigenvalue distribution, !
! from which a symmetric tridiagonal matrix will then be generated.            !
!                                                                              !
! Input data for DSTETESTER:                                                   !
! =========================                                                    !
!                                                                              !
! Input data for DSTETESTER are specified by means of macro commands (or key   !
! key words). These macros are groups of characters that uniquely define a     !
! specific subset of the input data, such as matrix types and dimensions,      !
! matrices to be read from files, etc. The input data are usually read from    !
! a file.                                                                      !
!                                                                              !
! Macros can be specified in any order, either in lower or upper case, with    !
! the corresponding subset of data they define. Data can be separated by       !
! blanks or commas, and the character % is interpreted as the beginning        !
! of a comment. Therefore, a line beginning with % is ignored.                 !
!                                                                              !
! In the directory INPUT the user should find three examples of input files:   !
!                                                                              !
!    input.test.easy                                                           !
!    input.test.mild                                                           !
!    input.test.hard                                                           !
!                                                                              !
! which exercise various features of DSTETESTER. After a successful            !
! compilation, the user should do, for example,                                !
!                                                                              !
!    dstetester.x < input.easy                                                 !
!                                                                              !
!******************************************************************************!
! The macros currently supported by DSTETESTER are described in                !
! in the file README.input                                                     !
!******************************************************************************!
!                                                                              !
! Data stored in derived data type M_DATA, module DSTEDEFINITIONS:             !
! ---------------------------------------------------------------              !
!                                                                              !
! form    type    size      en   id                                            !
!                                                                              !
!   1    ETYPE   ESIZE   GAMMA   built-in eigenvalue distribution              !
!   2    MTYPE   MSIZE   GAMMA   built-in tridiagonal matrix                   !
!   3        1    NEIG       0   eigenvalue distribution read from file        !
!   4        1     NDE       0   tridiagonal matrix read from file             !
!   5    STYPE   SSIZE   GAMMA   built-in singular value distribution          !
!   6    MTYPE   MSIZE   GAMMA   built-in bidiagonal matrix                    !
!   7        1    NSVL       0   singular value distribution read from file    !
!   8        1     NDE       0   bidiagonal matrix read from file              !
!   9        -       -       -   reserved for future use                       !
!                                                                              !
!==============================================================================!
!                                                                              !
! The user is welcome to incorporate additional test matrices in the routines  !
!    DSTEDEFVALS (built-in eigenvalue / singular value distributions),         !
!    DSTEDEFTMTRX (built-in tridiagonal matrices) and                          !
!    DSTEDEFBMTRX (built-in bidiagonal matrices)                               !
! in which case the user is encouraged to send a note to OAMarques@lbl.gov so  !
! we can consider incorporating other test matrices in our distribution.       !
!                                                                              !
!==============================================================================!
!
!.. Scalars ..
LOGICAL :: IEEE = .TRUE., SVD = .FALSE., TGK = .FALSE.
INTEGER :: ILAENV_10, ILAENV_11
INTEGER :: HBANDA = 100, HBANDR = 100, MAXN   = 0, NCASE = 0, &
           NILIU  =   0, NRILIU =   0, NRVLVU = 0, NVLVU = 0, &
           SETTO0 =   0 
REAL( KIND=PREC ) :: MINTT = 0.05, TOLEBZ, TOLEGR
!
!.. Static Arrays ..
LOGICAL, DIMENSION( 8 ) :: CALLST = .FALSE. 
LOGICAL, DIMENSION( 8 ) :: DUMP = .FALSE.
INTEGER, DIMENSION( 4 ) :: ISEED = 0
!
!.. Derived Data Types ..
TYPE( VALS_LIST ), POINTER :: VALS
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( I_LIST ), POINTER :: ILIU
TYPE( V_LIST ), POINTER :: VLVU
TYPE( M_LIST ), POINTER :: M
!
!.. External Functions ..
INTEGER, EXTERNAL :: ILAENV
REAL( KIND=PREC ), EXTERNAL :: DLAMCH
!
!.. Intrinsic Function ..
INTRINSIC SQRT
!
!.. Executable Statements ......................................................
!
ISEED = ISEED_INIT
!
! Set TOLEBZ, TOLEGR and check for IEEE arithmetic support.
!
TOLEBZ = 2*DLAMCH('Safe minimum')
TOLEGR = 2*DLAMCH('Safe minimum')
ILAENV_10 = ILAENV( 10, 'DSTETESTER', 'N', 1, 2, 3, 4 )
ILAENV_11 = ILAENV( 11, 'DSTETESTER', 'N', 1, 2, 3, 4 )
IEEE = ILAENV_10.EQ.1 .AND. ILAENV_11.EQ.1      
IF ( .NOT.IEEE ) THEN
   WRITE( UNIT=FUOUT, FMT='(A)' ) '** DSTETESTER warning: IEEE = .FALSE.'
   WRITE( UNIT=FUOUT, FMT='(A,I2)' ) &
   '** IEEE NaN arithmetic can be trusted not to trap =', ILAENV_10, &
   '** Infinity arithmetic can be trusted not to trap =', ILAENV_11
END IF
!
! Allocate derived data types.
!
ALLOCATE( VALS ); NULLIFY( VALS%NEXT )
ALLOCATE( MTRX ); NULLIFY( MTRX%NEXT )
ALLOCATE( ILIU ); NULLIFY( ILIU%NEXT )
ALLOCATE( VLVU ); NULLIFY( VLVU%NEXT )
ALLOCATE( M ); NULLIFY( M%NEXT )
!
! Parse the data file.
!
CALL DSTEGETDATA( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, &
                  MINTT, TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, &
                  NILIU, VLVU, NVLVU, VALS, MTRX, M )
!
! Perform the tests.
!
CALL DSTERUNTESTS( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, &
                   MINTT, TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, &
                   NILIU, VLVU, NVLVU, VALS, MTRX, M )
!
! Deallocate derived data types.
!
DEALLOCATE( VALS )
DEALLOCATE( MTRX )
DEALLOCATE( ILIU )
DEALLOCATE( VLVU )
DEALLOCATE( M )
!
END PROGRAM DSTETESTER
