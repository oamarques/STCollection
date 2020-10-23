SUBROUTINE SSTEMCMATRIX( MAXN, NCASE, MACRO, LIST, M )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETINTGR, GETSTRNG, PARSER
USE GSTEINTERFACES2, ONLY : PARSERLIST
USE SSTEDEFINITIONS
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: MAXN, NCASE
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
TYPE( M_LIST ), POINTER :: M
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  SSTEMCMATRIX deals with the macro that sets matrices to be generated.       !
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
!  M      (input/output) M_LIST (derived data type)                            !
!         Properties of the matrices to be used in the tests.                  !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: I, IERR, J, MFORM, NMSIZE, NMTYPE
!
!.. Static Arrays ..
CHARACTER( LEN=1 ) :: MATRIX_TYPE( 1 )
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: MSIZE( : ), MTYPE( : )
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. External Functions ..
CHARACTER( LEN=RECORD_LENGTH ), EXTERNAL :: GETRECORD
INTEGER, EXTERNAL :: LISTLENGTH
!
!.. Intrinsic Function ..
INTRINSIC MAX
!
!.. Executable Statements ......................................................
!
MATRIX_TYPE = GETSTRNG( MACRO, LIST, 1 )
!
SELECT CASE ( MATRIX_TYPE( 1 ) )
CASE ( 'T', 't' )
     MFORM = 2
CASE ( 'B', 'b' )
     MFORM = 6
CASE DEFAULT
     CALL HANDLER( 1, MACRO )
END SELECT

! Read matrix types.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NMTYPE = LISTLENGTH( LIST )
ALLOCATE( MTYPE( NMTYPE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'MTYPE (subroutine SSTEMCMATRIX)' )
MTYPE = GETINTGR( 'MTYPE', LIST, NMTYPE )
!
! Read matrix sizes.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD )
LIST => PARSERLIST( MACRO // '(invalid matrix size)', LIST )
NMSIZE = LISTLENGTH( LIST )
ALLOCATE( MSIZE( NMSIZE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'MSIZE (subroutine SSTEMCMATRIX)' )
MSIZE = GETINTGR( 'MSIZE', LIST, NMSIZE )
!
! Check data consistency.
!
IF      ( MINVAL( MTYPE ) < 0 ) THEN
        CALL HANDLER( 1, MACRO // '(matrix type < 0)' )
ELSE IF ( MINVAL( MSIZE ) < 0 ) THEN
        CALL HANDLER( 1, MACRO // '(matrix size < 0' )
END IF
!
! Store information in M_DATA.
!
DO I = 1, NMTYPE
   DO J = 1, NMSIZE
      M%DATA = M_DATA( MFORM, MTYPE( I ), MSIZE( J ), 1, 1, 0, &
                       ISEED_INIT( 1:4 ), ZERO )
      ALLOCATE( M%NEXT);  M => M%NEXT; NULLIFY( M%NEXT )
      MAXN = MAX( MAXN, MSIZE( J ) )
   END DO
END DO
!
NCASE = NCASE + NMTYPE*NMSIZE
!
! Deallocate arrays.
!
DEALLOCATE( MTYPE, MSIZE )
!
END SUBROUTINE SSTEMCMATRIX
