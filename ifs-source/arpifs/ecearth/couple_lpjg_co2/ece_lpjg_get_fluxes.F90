SUBROUTINE ECE_LPJG_GET_FLUXES(KDIM,YDMODEL,PSURF,YDSURF)

#ifdef WITH_CPLNG2

    !
    ! Called from compo_apply_emissions_layer.F90 when the flux update should take place.
    !
    USE PARKIND1   ,ONLY : JPIM, JPRB
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
    REAL(KIND=JPRB)    :: ZCONV, RMC

    IF (LHOOK) CALL DR_HOOK('ECE_LPJG_GET_FLUXES',0,ZHOOK_HANDLE)
    ! =========================================================================
    ! *** Pre-compute indices
    ! =========================================================================

    IL = KDIM%KIDIA
    IE = KDIM%KFDIA - KDIM%KIDIA
    IG = KDIM%KSTGLO - 1 + KDIM%KIDIA

    ! =========================================================================
    ! *** Write Vegetation CO2 fluxes into correct emissions field.
    ! =========================================================================

    ! convert from kg(C)/m2/day to kgCO2/m2/s
    RMC = 12.01115
    ZCONV = RMCO2 / RMC / ( 24 * 60 * 60 )
    !ZCONV = 0 ! for testing...

    ! Loop through emissions to find land CO2 fluxes
    DO IEMIS = 1, YDMODEL%YRML_CHEM%YRCOMPO%NEMIS2D_DESC
      IF (TRIM(YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%SPECIES) == "CO2_GHG") THEN
        IEMISID = YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%PARAM_INDEX
        IF (TRIM(YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%SECTOR) == "co2nbf") THEN
          PSURF%PSD_VF(:,YDSURF%YSD_VF%YEMIS2D(IEMISID)%MP) = CPLNG2_FLD(CPLNG2_IDX('FCO2NAT'))%D(IG:IG+IE,1,1) * ZCONV
        ELSEIF (TRIM(YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%SECTOR) == "co2apf") THEN
          PSURF%PSD_VF(:,YDSURF%YSD_VF%YEMIS2D(IEMISID)%MP) = CPLNG2_FLD(CPLNG2_IDX('FCO2ANT'))%D(IG:IG+IE,1,1) * ZCONV
        ELSEIF (TRIM(YDMODEL%YRML_CHEM%YRCOMPO%YEMIS2D_DESC(IEMIS)%SECTOR) == "co2fire") THEN
          PSURF%PSD_VF(:,YDSURF%YSD_VF%YEMIS2D(IEMISID)%MP) = CPLNG2_FLD(CPLNG2_IDX('FCO2NPP'))%D(IG:IG+IE,1,1) * ZCONV
        ENDIF
      ENDIF
    ENDDO


    IF (LHOOK) CALL DR_HOOK('ECE_LPJG_GET_FLUXES',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ECE_LPJG_GET_FLUXES
