SUBROUTINE ECE_PISCES_GET_FLUXES(KDIM,YDMODEL,PSURF,YDSURF)

#ifdef WITH_CPLNG2

    !
    ! Called from compo_apply_emissions_layer.F90 when the flux update should take place.
    !
    USE PARKIND1   ,ONLY : JPIM
    USE YOMHOOK    ,ONLY : LHOOK, DR_HOOK, JPHOOK
    USE YOMPHYDER  ,ONLY : DIMENSION_TYPE, SURF_AND_MORE_TYPE
    USE TYPE_MODEL ,ONLY : MODEL
    USE SURFACE_FIELDS_MIX ,ONLY : TSURF
    USE YOMCST     ,ONLY : RMCO2
    USE CPLNG2

    IMPLICIT NONE

    TYPE (DIMENSION_TYPE)          , INTENT(IN)    :: KDIM
    TYPE (MODEL)                   , INTENT(IN)    :: YDMODEL
    TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
    TYPE(TSURF)                    , INTENT(IN)    :: YDSURF

    ! Locals
    REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG
    INTEGER(KIND=JPIM) :: IEMIS, IEMISID

    IF (LHOOK) CALL DR_HOOK('ECE_PISCES_GET_FLUXES',0,ZHOOK_HANDLE)
    ! =========================================================================
    ! *** Pre-compute indices
    ! =========================================================================

    IL = KDIM%KIDIA
    IE = KDIM%KFDIA - KDIM%KIDIA
    IG = KDIM%KSTGLO - 1 + KDIM%KIDIA

    ! =========================================================================
    ! *** Write Ocean CO2 flux into correct emissions field.
    ! =========================================================================

    ! Loop through emissions to find CO2 ocean flux
    DO IEMIS = 1, YDMODEL%YRML_CHEM%YRCOMPO%NEMIS2D_DESC
      IF (TRIM(YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%SPECIES) == "CO2_GHG" .AND. &
        TRIM(YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%SECTOR) == "co2of") THEN
        IEMISID = YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%PARAM_INDEX
        EXIT
      ENDIF
    ENDDO

    ! convert from molC/m2/s to kgCO2/m2/s and reverse sign - since flux from PISCES is into the ocean by convention
    PSURF%PSD_VF(:,YDSURF%YSD_VF%YEMIS2D(IEMISID)%MP) = CPLNG2_FLD(CPLNG2_IDX('A_FCO2_oce'))%D(IG:IG+IE,1,1) * RMCO2 / 1000.0 * (-1)

    IF (LHOOK) CALL DR_HOOK('ECE_PISCES_GET_FLUXES',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ECE_PISCES_GET_FLUXES
