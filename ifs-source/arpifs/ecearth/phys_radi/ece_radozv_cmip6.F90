SUBROUTINE ECE_RADOZV_CMIP6(KIDIA, KFDIA, KLON, KLEV,&
 & KRINT, KDLON, KSHIFT,&
 & PAPRS, PGELAM, PGEMU,&
 & POZON)

!***********************************************************************
! CAUTION: THIS ROUTINE WORKS ONLY ON A NON-ROTATED, UNSTRETCHED GRID
!***********************************************************************

!**** *RADOZV_CMIP6* - COMPUTES VARIABLE OZONE DISTRIBUTION

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *RADOZV_CMIP6* FROM *RADINT*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     Shuting Yang  EC-Earth 9/02/2010
!     (after RADOZC, J.-J. MORCRETTE  E.C.M.W.F. 95/01/25)

!-----------------------------------------------------------------------

  USE PARKIND1, ONLY: JPIM, JPRB
  USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK
  USE ECE_YOEOZOV_CMIP6, ONLY: NLON1, NLAT1, NLV1, RLONCLI, RSINC1, ROZT1, RPROC1

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)    :: KLON
  INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV
  INTEGER(KIND=JPIM), INTENT(IN)    :: KDLON
  INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN)    :: KRINT
  INTEGER(KIND=JPIM), INTENT(IN)    :: KSHIFT
  REAL(KIND=JPRB), INTENT(IN)    :: PAPRS(KLON, KLEV + 1)
  REAL(KIND=JPRB), INTENT(IN)    :: PGELAM(KLON)
  REAL(KIND=JPRB), INTENT(IN)    :: PGEMU(KLON)
  REAL(KIND=JPRB), INTENT(OUT)   :: POZON(KDLON, KLEV)
!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------

  REAL(KIND=JPRB) :: ZOZLT(KDLON, 0:NLV1 + 1), ZOZON(KDLON, KLEV + 1)
  REAL(KIND=JPRB) :: ZRRR(0:NLV1)
  REAL(KIND=JPRB) :: ZWLON(KDLON)

  INTEGER(KIND=JPIM) :: IL, INLA, JC, JK, JL, JIR
  INTEGER(KIND=JPIM) :: IINLO1(KDLON), IINLO2(KDLON)

  REAL(KIND=JPRB) :: ZPMR, ZSILAT, ZSIN, ZDLONR, ZLON, ZOZLT1, ZOZLT2
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

!*         1.     LATITUDE INDEX WITHIN OZONE CLIMATOLOGY
!                 ---------------------------------------

  IF (LHOOK) CALL DR_HOOK('RADOZV_CMIP6', 0, ZHOOK_HANDLE)

  ZSIN = PGEMU(KIDIA)
  INLA = 0
  ZSILAT = -9999._JPRB
  DO JL = NLAT1 - 1, 1, -1
    IF (ZSIN <= RSINC1(JL + 1) .AND. ZSIN >= RSINC1(JL)) THEN
      INLA = JL
    END IF
  END DO

  IF (INLA == 0) THEN
    CALL ABOR1(' Problem with lat. interpolation in radozv_cmip6!')
  END IF
  ZSILAT = (ZSIN - RSINC1(INLA))/(RSINC1(INLA + 1) - RSINC1(INLA))

!
!        1a.     LONGITUDE INDEX WITHIN OZONE CLIMATOLOGY
!    ------------------------------------------------------------------
  ZDLONR = RLONCLI(2) - RLONCLI(1)
  IINLO1(:) = 0
  IINLO2(:) = 0

  IL = KSHIFT
  DO JL = KIDIA, KFDIA, KRINT
    IL = IL + 1
    ZLON = PGELAM(JL)
    DO JIR = 1, NLON1 - 1
      IF (ZLON < RLONCLI(JIR + 1) .AND. ZLON >= RLONCLI(JIR)) THEN
        ZWLON(IL) = (ZLON - RLONCLI(JIR))/ZDLONR
        IINLO1(IL) = JIR
        IINLO2(IL) = JIR + 1
      END IF
    END DO
    IF (ZLON >= RLONCLI(NLON1)) THEN
      ZWLON(IL) = (ZLON - RLONCLI(NLON1))/ZDLONR
      IINLO1(IL) = NLON1
      IINLO2(IL) = 1
    END IF
  END DO

!EC-EARTH: The following check is incorrect. None of the indexes maybe zero!
!IF (MAXVAL(IINLO1(:)) == 0 .OR. MAXVAL(IINLO2(:)) == 0) THEN
  IF (ANY(IINLO1(KSHIFT + 1:IL) == 0) .OR. ANY(IINLO2(KSHIFT + 1:IL) == 0)) THEN
    CALL ABOR1(' Problem with lon. interpolation in radozv_cmip6!')
  END IF

!     ------------------------------------------------------------------

!*         2.     LATITUDE/LONGITUDE INTERPOLATED FIELD
!                 -------------------------------------

  IF (INLA == NLAT1 .OR. INLA == 1) THEN
    DO JC = 0, NLV1 + 1
      IL = KSHIFT
      DO JL = KIDIA, KFDIA, KRINT
        IL = IL + 1
        ZOZLT(IL, JC) = ROZT1(IINLO1(IL), INLA, JC) + ZWLON(IL)* &
          &         (ROZT1(IINLO2(IL), INLA, JC) - ROZT1(IINLO1(IL), INLA, JC))
      END DO
    END DO
  ELSE
    DO JC = 0, NLV1 + 1
      IL = KSHIFT
      DO JL = KIDIA, KFDIA, KRINT
        IL = IL + 1
        ZOZLT1 = ROZT1(IINLO1(IL), INLA, JC) + ZWLON(IL)* &
          &   (ROZT1(IINLO2(IL), INLA, JC) - ROZT1(IINLO1(IL), INLA, JC))
        ZOZLT2 = ROZT1(IINLO1(IL), INLA + 1, JC) + ZWLON(IL)* &
          &   (ROZT1(IINLO2(IL), INLA + 1, JC) - ROZT1(IINLO1(IL), INLA + 1, JC))
        ZOZLT(IL, JC) = ZOZLT1 + ZSILAT*(ZOZLT2 - ZOZLT1)
      END DO
    END DO
  END IF

!     ------------------------------------------------------------------

!*         3.     VERTICAL INTERPOLATION
!                 ----------------------

  DO JC = 0, NLV1
    ZRRR(JC) = (1.0_JPRB/(RPROC1(JC) - RPROC1(JC + 1)))
  END DO

  IL = KSHIFT
  DO JL = KIDIA, KFDIA, KRINT
    IL = IL + 1
    DO JC = 0, NLV1
      DO JK = 1, KLEV + 1
        ZPMR = PAPRS(JL, JK)
        IF (ZPMR >= RPROC1(JC) .AND. ZPMR < RPROC1(JC + 1)) &
         & ZOZON(IL, JK) = ZOZLT(IL, JC + 1) + (ZPMR - RPROC1(JC + 1))*ZRRR(JC)* &
         &             (ZOZLT(IL, JC) - ZOZLT(IL, JC + 1))
      END DO
    END DO
  END DO

  IL = KSHIFT
  DO JL = KIDIA, KFDIA, KRINT
    IL = IL + 1
    DO JK = 1, KLEV + 1
      ZPMR = PAPRS(JL, JK)
      IF (ZPMR >= RPROC1(NLV1 + 1)) ZOZON(IL, JK) = ZOZLT(IL, NLV1 + 1)
    END DO
  END DO

! INTEGRATION IN THE VERTICAL:
  IL = KSHIFT
  DO JL = KIDIA, KFDIA, KRINT
    IL = IL + 1
    DO JK = 1, KLEV
      POZON(IL, JK) = (PAPRS(JL, JK + 1) - PAPRS(JL, JK))&
                     & *(ZOZON(IL, JK) + ZOZON(IL, JK + 1))*0.5_JPRB
    END DO
  END DO

!     -----------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RADOZV_CMIP6', 1, ZHOOK_HANDLE)
END SUBROUTINE ECE_RADOZV_CMIP6
