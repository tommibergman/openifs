SUBROUTINE ECE_SI3_GET_ICE_STATE(KSTGLO,KIDIA,KFDIA,ICE_ALBEDO,ICE_THICKNESS,SNOW_THICKNESS)

    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

    USE YOMPHYDER,ONLY: DIMENSION_TYPE

    USE CPLNG

    IMPLICIT NONE

    ! Arguments
    INTEGER(KIND=JPIM),INTENT(IN)           :: KSTGLO,KIDIA,KFDIA
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_ALBEDO(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_THICKNESS(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: SNOW_THICKNESS(KIDIA:KFDIA)

    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG

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
        ICE_ALBEDO = CPLNG_FLD(CPLNG_IDX('A_Ice_albedo'))%D(IG:IG+IE,1,1)
    ENDIF

    IF (LHOOK) CALL DR_HOOK('ECE_SI3_GET_ICE_STATE',1,ZHOOK_HANDLE)

END SUBROUTINE ECE_SI3_GET_ICE_STATE
