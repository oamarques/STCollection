SUBROUTINE SSTEMCMATRIXF( MAXN, NCASE, MACRO, LIST, MTRX, M )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETSTRNG
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: MAXN, NCASE
!
!.. Derived Data Type Arguments ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
TYPE( MTRX_LIST ), POINTER :: MTRX
TYPE( M_LIST ), POINTER :: M
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCMATRIXF deals with the macro that reads a matrix from a file.         !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
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
!  MTRX   (input/output) MTRX_LIST (derived data type)                         !
!         Matrices read from files.                                            !
!                                                                              !
!  M      (input/output) M_LIST (derived data type)                            !
!         Properties of the matrices to be used in the tests.                  !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
LOGICAL :: FILE_EXIST
INTEGER :: I, IERR, J, MFORM, N
REAL( KIND=PREC ) :: TEMP1, TEMP2
!
!.. Static Array ..
CHARACTER( LEN=FILE_NAME_LENGTH ) :: FILE_DATA( 2 )
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Intrinsic Functions ..
INTRINSIC MAX, TRIM
!
!.. Executable Statements ......................................................
!
FILE_DATA = GETSTRNG( MACRO, LIST, 2 )
!
SELECT CASE ( TRIM( FILE_DATA( 1 ) ) )
CASE ( 'T', 't' )
     MFORM = 4
CASE ( 'B', 'b' )
     MFORM = 8
CASE ( 'Z', 'z' )
     MFORM = 9
CASE DEFAULT
     CALL HANDLER( 1, MACRO )
END SELECT
!
INQUIRE ( FILE=FILE_DATA( 2 ), EXIST=FILE_EXIST )
!
IF ( FILE_EXIST ) THEN
   OPEN( UNIT=FUEXT, FILE=FILE_DATA( 2 ), IOSTAT=IERR )
   IF ( IERR /= 0 ) CALL HANDLER( 4, FILE_DATA( 2 ) )
   READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) N
   IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_DATA( 2 ) )
   ALLOCATE( MTRX%D( N ), MTRX%E( N ), STAT=IERR )
   DO I = 1, N
      READ ( UNIT=FUEXT, FMT=*, IOSTAT=IERR ) J, TEMP1, TEMP2
      IF ( IERR /= 0 ) CALL HANDLER( 1, FILE_DATA( 2 ) )
      MTRX%D( J ) = TEMP1; MTRX%E( J ) = TEMP2
   END DO
   CLOSE( UNIT=FUEXT )
   M%DATA = M_DATA( MFORM, 1, N, 1, 1, 0, ISEED_INIT( 1:4 ), ZERO )
   ALLOCATE( M%NEXT);  M => M%NEXT; NULLIFY( M%NEXT )
   ALLOCATE( MTRX%NEXT );  MTRX => MTRX%NEXT; NULLIFY( MTRX%NEXT )
ELSE
   CALL HANDLER( 4, FILE_DATA( 2 ) )
END IF
!
MAXN = MAX( MAXN, N )
NCASE = NCASE + 1
!
END SUBROUTINE SSTEMCMATRIXF
