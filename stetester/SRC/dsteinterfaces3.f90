MODULE DSTEINTERFACES3
!
!==============================================================================!
!                                                                              !
! This module defines interfaces for top level subroutines.                    !
!                                                                              !
!==============================================================================!
!
!-------------------------------------------------------------------------------
INTERFACE DSTEGETDATA
   SUBROUTINE DSTEGETDATA( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, &
                           MINTT, TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, &
                           NILIU, VLVU, NVLVU, VALS, MTRX, M )
      USE GSTEDEFINITIONS
      USE GSTEINTERFACES1
      USE DSTEDEFINITIONS
      USE DSTEINTERFACES1
      USE DSTEINTERFACES2
      LOGICAL :: CALLST( 8 ), DUMP( 8 ), SVD, TGK
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU, SETTO0
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( VALS_LIST ), POINTER :: VALS
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( M_LIST ), POINTER :: M
   END SUBROUTINE DSTEGETDATA
END INTERFACE
!-------------------------------------------------------------------------------
INTERFACE DSTERUNTESTS
   SUBROUTINE DSTERUNTESTS( HBANDA, HBANDR, SETTO0, MAXN, NCASE, NRILIU, NRVLVU, MINTT, &
                            TOLEBZ, TOLEGR, CALLST, DUMP, SVD, TGK, ILIU, NILIU, &
                            VLVU, NVLVU, VALS, MTRX, M )
      USE GSTEDEFINITIONS
      USE DSTEDEFINITIONS
      USE DSTEINTERFACES2, ONLY : DBDPRNRSLT, DSTEBNDGAP, DSTEMATGEN, DSTEPRNRSLT
      LOGICAL :: CALLST( 8 ), DUMP( 8 ), SVD, TGK
      INTEGER :: HBANDA, HBANDR, MAXN, NCASE, NILIU, NRILIU, NRVLVU, NVLVU, SETTO0
      REAL( KIND=PREC ) :: MINTT, TOLEBZ, TOLEGR
      TYPE( VALS_LIST ), POINTER :: VALS
      TYPE( MTRX_LIST ), POINTER :: MTRX
      TYPE( I_LIST ), POINTER :: ILIU
      TYPE( V_LIST ), POINTER :: VLVU
      TYPE( M_LIST ), POINTER :: M
   END SUBROUTINE DSTERUNTESTS
END INTERFACE
!-------------------------------------------------------------------------------
!
END MODULE DSTEINTERFACES3
