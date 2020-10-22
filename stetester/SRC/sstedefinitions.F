MODULE SSTEDEFINITIONS
!
!==============================================================================!
!                                                                              !
! This module sets variables and defines derived data types.                   !
!                                                                              !
!==============================================================================!
!
!*******************************************
!**                                       **
!**  The following lines may be modified  **
!**                                       **
!*******************************************
!
INTEGER, PARAMETER :: PREC = KIND( 1.0E0 )
!
REAL( KIND=PREC ), PARAMETER :: HALF  = 0.5E0
REAL( KIND=PREC ), PARAMETER :: HNDRD = 1.0E2
REAL( KIND=PREC ), PARAMETER :: ONE   = 1.0E0
REAL( KIND=PREC ), PARAMETER :: TEN   = 1.0E1
REAL( KIND=PREC ), PARAMETER :: THSND = 1.0E3
REAL( KIND=PREC ), PARAMETER :: TWO   = 2.0E0
REAL( KIND=PREC ), PARAMETER :: ZERO  = 0.0E0
!
REAL( KIND=PREC ), PARAMETER :: MICRO = 1.0E-6
!
CHARACTER( LEN=19 ), DIMENSION( 1:5 ) :: FDUMP  = (/ 'sstetester_dump.T  ', &
                                                     'sstetester_dump.W  ', &
                                                     'sstetester_dump.Z  ', &
                                                     'sstetester_dump.log', &
                                                     'sstetester_dump.m  ' /)
CHARACTER( LEN=9 ), DIMENSION( 1:8 )  :: IDTEST = (/ 'SSTEQR(I)', &
                                                     'SSTEVX(A)', &
                                                     'SSTEVX(I)', &
                                                     'SSTEVX(V)', &
                                                     'SSTEDC(I)', &
                                                     'SSTEGR(A)', &
                                                     'SSTEGR(I)', &
                                                     'SSTEGR(V)' /)
!
!*****************************************
!**                                     **
!**  Do not modify the following lines  **
!**                                     **
!*****************************************
!
! Derived data type for tridiagonals.
!
TYPE M_DATA
     INTEGER :: FORM, TYPE, SIZE, COND, DIST, SIGN, SEED( 4 )
     REAL( KIND=PREC ) :: EN
END TYPE M_DATA
!
! List of tridiagonals.
!
TYPE M_LIST
     TYPE( M_DATA ) :: DATA
     TYPE( M_LIST ), POINTER :: NEXT
END TYPE M_LIST
!
! List of values read from files.
!
TYPE VALS_LIST
     REAL( KIND=PREC ), POINTER :: S( : )
     TYPE( VALS_LIST ), POINTER :: NEXT
END TYPE VALS_LIST
!
! List of tridiagonal matrices read from files.
!
TYPE MTRX_LIST
     REAL( KIND=PREC ), POINTER :: D( : ), E( : )
     TYPE( MTRX_LIST ), POINTER :: NEXT
END TYPE MTRX_LIST
!
! List of ranges, indexes.
!
TYPE I_LIST
     INTEGER :: IL, IU
     TYPE( I_LIST ), POINTER :: NEXT
END TYPE I_LIST
!
! List of ranges, values.
!
TYPE V_LIST
     REAL( KIND=PREC ) :: VL, VU
     TYPE( V_LIST ), POINTER :: NEXT
END TYPE V_LIST
!
END MODULE SSTEDEFINITIONS
