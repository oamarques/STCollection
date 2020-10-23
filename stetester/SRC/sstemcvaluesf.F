SUBROUTINE SSTEMCVALUESF( ISEED, MAXN, NCASE, MACRO, LIST, VALS, KEY, M )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETSTRNG
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: KEY, MAXN, NCASE
!
!.. Array Argument ..
INTEGER :: ISEED( 4 )
!
!.. Derived Data Type Arguments ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
TYPE( VALS_LIST ), POINTER :: VALS
TYPE( M_LIST ), POINTER :: M
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCVALUESF deals with the macro that reads an eigenvalue distribution    !
!  from a file.                                                                !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  ISEED  (input/output) INTEGER array, dimension ( 4 )                        !
!         Seed for the random number generator. Each entry of ISEED should     !
!         lie between 0 and 4095 inclusive and ISEED(4) should be odd.         !
!                                                                              !
!  MAXN   (input/output) INTEGER                                               !
!         Maximum dimension of the matrix to be tested.                        !
!                                                                              !
!  NCASE  (input/output) INTEGER                                               !
!         Number of matrices defined in M.                                     !
!                                                                              !
!  MACRO  (input) CHARACTER                                                    !
!         Macro definition.                                                    !
!                                                                              !
!  LIST   (input) DATA_FROM_RECORD (derived data type)                         !
!         String (file name).                                                  !
!                                                                              !
!  VALS   (input/output) VALS_LIST (derived data type)                         !
!         Values read from files.                                              !
!                                                                              !
!  KEY    (input) INTEGER                                                      !
!         =3, eigenvalues; = 7, singular values                                !
!                                                                              !
!  M      (input/output) M_LIST (derived data type)                            !
!         Properties of the matrices to be used in the tests.                  !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
LOGICAL :: FILE_EXIST
INTEGER :: I, IERR, N
REAL( KIND=PREC ) :: TEMP
! 
!.. Static Array ..
CHARACTER( LEN=FILE_NAME_LENGTH ) :: FILE_NAME( 1 )
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
FILE_NAME = GETSTRNG( MACRO, LIST, 1 )
INQUIRE ( FILE=FILE_NAME( 1 ), EXIST=FILE_EXIST )
!
IF ( FILE_EXIST ) THEN
   OPEN( UNIT=FUEXT, FILE=FILE_NAME( 1 ), IOSTAT=IERR )
   IF ( IERR /= 0 ) CALL HANDLER( 4, FILE_NAME( 1 ) )
   READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) N
   IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_NAME( 1 ) )
   ALLOCATE( VALS%S( N ) )
   DO I = 1, N
      READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) TEMP
      IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_NAME( 1 ) )
      VALS%S( I ) = TEMP
   END DO
   CLOSE( UNIT=FUEXT )
   M%DATA = M_DATA( KEY, 1, N, 1, 1, 0, ISEED( 1:4 ), ZERO )
   ALLOCATE( M%NEXT );  M => M%NEXT; NULLIFY( M%NEXT )
   ALLOCATE( VALS%NEXT );  VALS => VALS%NEXT; NULLIFY( VALS%NEXT )
ELSE
   CALL HANDLER( 4, FILE_NAME( 1 ) )
END IF
!
MAXN = MAX( MAXN, N )
NCASE = NCASE + 1
!
END SUBROUTINE SSTEMCVALUESF
