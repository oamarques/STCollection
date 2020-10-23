FUNCTION GETSTRNG( STRING, LIST, N )
!
USE GSTEDEFINITIONS
!
!.. Scalar Arguments ..
CHARACTER( * ) :: STRING
INTEGER :: N
!
!.. Derived Data Type Argument ..
TYPE( DATA_FROM_RECORD ), TARGET :: LIST
!
!.. Function Result ..
CHARACTER( LEN=FIELD_LENGTH ) :: GETSTRNG( N )
!
!==============================================================================!
!                                                                              !
!  Purpose:                                                                    !
!  =======                                                                     !
!                                                                              !
!  GETSTRNG extracts strings from RECORD.                                      !
!                                                                              !
!  Arguments:                                                                  !
!  =========                                                                   !
!                                                                              !
!  STRING  (input) CHARACTER                                                   !
!          String to be printed in case of error.                              !
!                                                                              !
!  LIST    (input) DATA_FROM_RECORD                                            !
!          List of strings.                                                    !
!                                                                              !
!  N       (input) INTEGER                                                     !
!          Number of strings in LIST.                                          !
!                                                                              !
!==============================================================================!
!  
!.. Local Scalars ..
INTEGER :: I, IOERR
!
!.. Local Derived Data Type ..
TYPE( DATA_FROM_RECORD ), POINTER :: LIST_CURRENT
!
!.. External Subroutine ..
EXTERNAL HANDLER
!
!.. Executable Statements ......................................................
!
GETSTRNG( 1:N ) = ' '
LIST_CURRENT => LIST
!
DO I = 1, N
   READ ( LIST_CURRENT%FIELD, FMT='(A)', IOSTAT=IOERR ) GETSTRNG( I )
   IF ( IOERR /= 0 ) CALL HANDLER( 3, STRING )
   LIST_CURRENT => LIST_CURRENT%NEXT
END DO
!
END FUNCTION GETSTRNG
