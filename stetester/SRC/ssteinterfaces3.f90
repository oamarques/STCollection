MODULE SSTEINTERFACES3
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for top level subroutines.                    !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE SSTEGETDATA
   SUBROUTINE SSTEGETDATA( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, &
                           MINTT, TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, &
                           NILIU, VLVU, NVLVU, VALS, MTRX, M )
      USE GSTEDEFINITIONS
      USE GSTEINTERFACES1
      USE SSTEDEFINITIONS
      USE SSTEINTERFACES1
      USE SSTEINTERFACES2
      LOGICAL :: CALLST( 8 ), DUMP( 8 ), SVD, TGK
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU, SETTO0
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( VALS_LIST ), POINTER :: VALS
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( M_LIST ), POINTER :: M
   END SUBROUTINE SSTEGETDATA
END INTERFACE
!-------------------------------------------------------------------------------
INTERFACE SSTERUNTESTS
   SUBROUTINE SSTERUNTESTS( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                            TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, NILIU, &
                            VLVU, NVLVU, VALS, MTRX, M )
      USE GSTEDEFINITIONS
      USE SSTEDEFINITIONS
      USE SSTEINTERFACES2, ONLY : SBDPRNRSLT, SSTEBNDGAP, SSTEMATGEN, SSTEPRNRSLT
      LOGICAL :: CALLST( 8 ), DUMP( 8 ), SVD, TGK
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU, SETTO0
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( VALS_LIST ), POINTER :: VALS
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( M_LIST ), POINTER :: M
   END SUBROUTINE SSTERUNTESTS
END INTERFACE
!-------------------------------------------------------------------------------
!
END MODULE SSTEINTERFACES3
