SUBROUTINE DSTEMCGLUED( ECOND, EDIST, ESIGN, ISEED, MAXN, NCASE, MACRO, M )
!
USE GSTEDEFINITIONS
USE GSTEINTERFACES1, ONLY : GETINTGR, PARSER
USE DSTEDEFINITIONS
USE DSTEINTERFACES1, ONLY : GETDREAL
! 
!.. Scalar Arguments ..
CHARACTER( LEN=MACRO_NAME_LENGTH ) :: MACRO
INTEGER :: ECOND, EDIST, ESIGN, MAXN, NCASE
!
!.. Array Argument ..
INTEGER :: ISEED( 4 )
!
!.. Derived Data Type Argument ..
TYPE( M_LIST ), POINTER :: M
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  DSTEMCGLUED deals with the macro that defines glued matrices.               !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  ECOND  (input) INTEGER                                                      !
!         Sets the condition number for ETYPE = 1:5,                           !
!         = 1, COND = 1 / SQRT(ULP), default                                   !
!         = 2, COND = 1 / (N*SQRT(ULP)),                                       !
!         = 3, COND = 1 / (10*N*SQRT(ULP)),                                    !
!         = 4, COND = 1 / ULP,                                                 !
!         = 5, COND = 1 / (N*ULP),                                             !
!         = 6, COND = 1 / (10*N*ULP),                                          !
!                                                                              !
!  EDIST  (input) INTEGER                                                      !
!         Specifies the type of the distribution to be used in the generation  !
!         of random eigenvalues                                                !
!         = 1, UNIFORM(  0, 1 )                                                !
!         = 2, UNIFORM( -1, 1 )                                                !
!         = 3, NORMAL( 0, 1 )                                                  !
!                                                                              !
!  ESIGN  (input) INTEGER                                                      !
!         Attributes signs to the entries in W                                 !
!         = 0, the entries of W are unchanged                                  !
!         = 1, multiplies each entry of W by +/- 1 with probability 0.5        !
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
!  M      (input/output) M_LIST (derived data type)                            !
!         Properties of the matrices to be used in the tests.                  !
!                                                                              !
!==============================================================================!
! 
!.. Local Scalars ..
CHARACTER( LEN=RECORD_LENGTH ) :: RECORD
INTEGER :: I, IERR, NGAMMA, NGFORM, NGSIZE, NGTYPE
REAL( KIND=PREC ) :: TEMP 
!
!.. Allocatable Arrays ..
INTEGER, ALLOCATABLE :: GFORM( : ), GSIZE( : ), GTYPE( : )
REAL( KIND=PREC ), ALLOCATABLE :: GAMMA( : )
!
!.. Derived Data Types ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. External Functions ..
CHARACTER( LEN=RECORD_LENGTH ), EXTERNAL :: GETRECORD
INTEGER, EXTERNAL :: LISTLENGTH
REAL( KIND=PREC ), EXTERNAL :: DLARAN
!
!.. Intrinsic Function ..
INTRINSIC MAX, MIN
!
!.. Executable Statements ......................................................
!
! Read matrix modes.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NGFORM = LISTLENGTH( LIST )
ALLOCATE( GFORM( NGFORM ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'GFORM (subroutine DSTEMCGLUED)' ) 
GFORM = GETINTGR( MACRO, LIST, NGFORM )
!
! Read matrix types.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NGTYPE = LISTLENGTH( LIST )
ALLOCATE( GTYPE( NGTYPE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'GTYPE (subroutine DSTEMCGLUED)' ) 
GTYPE = GETINTGR( MACRO, LIST, NGTYPE )
!
! Read matrix sizes.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NGSIZE = LISTLENGTH( LIST )
ALLOCATE( GSIZE( NGSIZE ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'GSIZE (subroutine DSTEMCGLUED)' ) 
GSIZE = GETINTGR( MACRO, LIST, NGSIZE )
!
! Read glue factors.
!
RECORD = GETRECORD( )
LIST => PARSER( RECORD ); NGAMMA = LISTLENGTH( LIST )
ALLOCATE( GAMMA( NGAMMA+1 ), STAT=IERR )
IF ( IERR /= 0 ) CALL HANDLER( 2, 'GAMMA (subroutine DSTEMCGLUED)' ) 
GAMMA( 1:NGAMMA ) = GETDREAL( MACRO, LIST, NGAMMA )
GAMMA( NGAMMA+1 ) = ZERO
!
! Check data consistency.
!
IF      ( MINVAL( GFORM ) < 1 ) THEN
        CALL HANDLER( 1, MACRO // '(glued matrix form < 1)' )
ELSE IF ( MAXVAL( GFORM ) /= 1 .AND. &
          MAXVAL( GFORM ) /= 2 .AND. &
          MAXVAL( GFORM ) /= 5 .AND. &
          MAXVAL( GFORM ) /= 6 ) THEN
        CALL HANDLER( 1, MACRO // '(glued matrix form not allowed)' )
ELSE IF ( MINVAL( GTYPE ) < -9 ) THEN
        CALL HANDLER( 1, MACRO // '(glued matrix type not implemented)' )
ELSE IF ( MINVAL( GSIZE ) < 0 ) THEN
        CALL HANDLER( 1, MACRO // '(glued matrix size < 0)' )
ELSE IF ( MINVAL( ABS( GAMMA(1:NGAMMA) ) ) == ZERO ) THEN
        CALL HANDLER( 1, MACRO // '(glue factor = 0.0)' )
ELSE IF ( MIN( NGFORM, NGTYPE, NGSIZE, NGAMMA + 1 ) /= &
          MAX( NGFORM, NGTYPE, NGSIZE, NGAMMA + 1 ) ) THEN
        CALL HANDLER( 1, MACRO // '(number of parameters differ)' )
END IF
MAXN = MAX( MAXN, SUM( GSIZE ) )
!
! Store information in M_DATA.
!
DO I = 1, NGSIZE
   IF ( GFORM( I ) == 1 .OR. GFORM( I ) == 5  ) THEN
      M%DATA = M_DATA( GFORM( I ), GTYPE( I ), GSIZE( I ), ECOND, &
                       EDIST, ESIGN, ISEED( 1:4 ), GAMMA( I ) )  
   ELSE
      M%DATA = M_DATA( GFORM( I ), GTYPE( I ), GSIZE( I ), 1, &
                       1, 0, ISEED( 1:4 ), GAMMA( I ) )
   END IF
   ALLOCATE( M%NEXT);  M => M%NEXT; NULLIFY( M%NEXT )
   TEMP = DLARAN( ISEED )
END DO
!
NCASE = NCASE + 1
!
! Deallocate arrays.
!
DEALLOCATE( GFORM, GTYPE, GSIZE, GAMMA )
!
END SUBROUTINE DSTEMCGLUED
