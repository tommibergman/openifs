SUBROUTINE ECE_SI3_GET_ICE_STATE(KSTGLO,KIDIA,KFDIA,ICE_ALBEDO,ICE_THICKNESS,SNOW_THICKNESS)

#ifdef WITH_CPLNG2

    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

    USE YOMPHYDER,ONLY: DIMENSION_TYPE
    USE YOEPHY,   ONLY: YREPHY

    USE ECEARTH,  ONLY: ECE_CPL_NEMO_WEIGHTED_ICE
    USE CPLNG2

    IMPLICIT NONE

#include "surf_inq.h"

    ! Arguments
    INTEGER(KIND=JPIM),INTENT(IN)           :: KSTGLO,KIDIA,KFDIA
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_ALBEDO(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_THICKNESS(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: SNOW_THICKNESS(KIDIA:KFDIA)

    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG,JL

    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_FRAC(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_ALB(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_SNW_TCK(:)

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
        ICE_THICKNESS = CPLNG2_FLD(CPLNG2_IDX('A_Ice_thickness'))%D(IG:IG+IE,1,1)
    ENDIF

    IF(PRESENT(SNOW_THICKNESS)) THEN
        SNOW_THICKNESS = CPLNG2_FLD(CPLNG2_IDX('A_Snow_thickness'))%D(IG:IG+IE,1,1)
    ENDIF

    IF(PRESENT(ICE_ALBEDO)) THEN
        IF (ECE_CPL_NEMO_WEIGHTED_ICE) THEN
            CALL SURF_INQ(YREPHY%YSURF, PRCIMIN=ZRCIMIN, PRALBSEAD=ZRALBSEAD)
            CPL_FLD_ICE_FRAC => CPLNG2_FLD(CPLNG2_IDX('A_Ice_frac'))%D(IG:IG+IE,1,1)
            CPL_FLD_ICE_ALB => CPLNG2_FLD(CPLNG2_IDX('A_Ice_albedo'))%D(IG:IG+IE,1,1)
            DO JL = 1, SIZE(CPL_FLD_ICE_FRAC)
                ! check if ice fraction is larger than 0 to avoid division by 0 
                ! otherwise set constant albedo
                IF (CPL_FLD_ICE_FRAC(JL) > ZRCIMIN .AND. CPL_FLD_ICE_FRAC(JL) > EPSILON(1._JPRB)) THEN
                    ICE_ALBEDO(JL) = MAX(ZRALBSEAD, MIN(1._JPRB, CPL_FLD_ICE_ALB(JL) / CPL_FLD_ICE_FRAC(JL)))
                ELSE
                    ICE_ALBEDO(JL) = ZRALBSEAD
                END IF
                ! ensure albedo 0 <= albedo <= 1
                ICE_ALBEDO(JL) = MIN(1.0_JPRB, MAX(0.0_JPRB, ICE_ALBEDO(JL)))
            END DO
        ELSE
            ICE_ALBEDO = CPLNG2_FLD(CPLNG2_IDX('A_Ice_albedo'))%D(IG:IG+IE,1,1)
        ENDIF
    ENDIF

    IF (LHOOK) CALL DR_HOOK('ECE_SI3_GET_ICE_STATE',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ECE_SI3_GET_ICE_STATE
