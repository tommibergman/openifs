SUBROUTINE ECE_FESIM_GET_ICE_STATE(KSTGLO,KIDIA,KFDIA,SNOW_THICKNESS,ICE_ALBEDO)

#ifdef WITH_CPLNG2

    !
    ! Called from callpar.F90 when the radiation scheme needs the ice state
    ! Ice fraction, sst, and ice temp have already been updated by awi_updclie.
    ! Here we update ice thickness, snow thickness
    !
    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMHOOK,  ONLY: LHOOK, DR_HOOK, JPHOOK

    USE YOMPHYDER,ONLY: DIMENSION_TYPE

    USE CPLNG2

    IMPLICIT NONE

    ! Arguments
    INTEGER(KIND=JPIM),INTENT(IN)             :: KSTGLO,KIDIA,KFDIA
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: SNOW_THICKNESS(KIDIA:KFDIA)
    REAL(KIND=JPRB),   INTENT(OUT),OPTIONAL :: ICE_ALBEDO(KIDIA:KFDIA)


    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG

    IF (LHOOK) CALL DR_HOOK('ECE_FESIM_GET_ICE_STATE',0,ZHOOK_HANDLE)
    ! =========================================================================
    ! *** Pre-compute indices
    ! =========================================================================

    IL = KIDIA
    IE = KFDIA - KIDIA
    IG = KSTGLO - 1 + KIDIA

    ! =========================================================================
    ! *** Ice thickness, snow thickness
    !     SST, ice fraction, ice temperature are updated in awi_updclie.F90
    ! =========================================================================

    IF(PRESENT(SNOW_THICKNESS)) SNOW_THICKNESS(:) = CPLNG2_FLD(CPLNG2_IDX('A_Snow_thickness'))%D(IG:IG+IE,1,1)
    IF(PRESENT(ICE_ALBEDO)) ICE_ALBEDO(:) = CPLNG2_FLD(CPLNG2_IDX('A_Ice_albedo'))%D(IG:IG+IE,1,1)

    IF (LHOOK) CALL DR_HOOK('ECE_FESIM_GET_ICE_STATE',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ECE_FESIM_GET_ICE_STATE
