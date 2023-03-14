SUBROUTINE ECE_SI3_GET_ICE_STATE(KSTGLO,KIDIA,KFDIA,ICE_ALBEDO,ICE_THICKNESS,SNOW_THICKNESS)

    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

    USE YOMPHYDER,ONLY: DIMENSION_TYPE
    USE YOEPHY,   ONLY: YREPHY

    USE ECEARTH,  ONLY: ECE_CPL_NEMO_WEIGHTED_ICE
    USE CPLNG

    IMPLICIT NONE

#include "surf_inq.h"

    ! Arguments
    INTEGER(KIND=JPIM),INTENT(IN)           :: KSTGLO,KIDIA,KFDIA
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_ALBEDO(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_THICKNESS(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: SNOW_THICKNESS(KIDIA:KFDIA)

    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG

    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_FRAC(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_ALB(:)

    REAL(KIND=JPRB) :: ZRCIMIN, ZRALBSEAD

    IF (LHOOK) CALL DR_HOOK('ECE_SI3_GET_ICE_STATE',0,ZHOOK_HANDLE)

    ! =========================================================================
    ! *** Pre-compute indices
    ! =========================================================================

    IL = KIDIA
    IE = KFDIA - KIDIA
    IG = KSTGLO - 1 + KIDIA

    ! =========================================================================
    ! *** Set ice/snow thickness and ice albedo from coupling fields
    ! =========================================================================

    IF(PRESENT(ICE_THICKNESS)) THEN
        ICE_THICKNESS = CPLNG_FLD(CPLNG_IDX('A_Ice_thickness'))%D(IG:IG+IE,1,1)
    ENDIF

    IF(PRESENT(SNOW_THICKNESS)) THEN
        SNOW_THICKNESS = CPLNG_FLD(CPLNG_IDX('A_Snow_thickness'))%D(IG:IG+IE,1,1)
    ENDIF

    IF(PRESENT(ICE_ALBEDO)) THEN
        IF (ECE_CPL_NEMO_WEIGHTED_ICE) THEN
            CALL SURF_INQ(YREPHY%YSURF, PRCIMIN=ZRCIMIN, PRALBSEAD=ZRALBSEAD)
            CPL_FLD_ICE_FRAC => CPLNG_FLD(CPLNG_IDX('A_Ice_frac'))%D(IG:IG+IE,1,1)
            CPL_FLD_ICE_ALB => CPLNG_FLD(CPLNG_IDX('A_Ice_albedo'))%D(IG:IG+IE,1,1)
            WHERE ( CPL_FLD_ICE_FRAC > ZRCIMIN )
                ICE_ALBEDO = MAX(ZRALBSEAD,MIN(1._JPRB,CPL_FLD_ICE_ALB/CPL_FLD_ICE_FRAC))
            ELSEWHERE
                ICE_ALBEDO = ZRALBSEAD
            ENDWHERE
        ELSE
            ICE_ALBEDO = CPLNG_FLD(CPLNG_IDX('A_Ice_albedo'))%D(IG:IG+IE,1,1)
        ENDIF
    ENDIF

    IF (LHOOK) CALL DR_HOOK('ECE_SI3_GET_ICE_STATE',1,ZHOOK_HANDLE)

END SUBROUTINE ECE_SI3_GET_ICE_STATE
