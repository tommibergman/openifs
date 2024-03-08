SUBROUTINE ECE_NEMO_SET_OCEAN_FLUXES(YDSURF, KDIM, SURFL, PSURF, FLUX)

    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMPHY2, ONLY : YRPHY2
    USE SURFACE_FIELDS_MIX, ONLY: TSURF
    USE YOMPHYDER, ONLY: DIMENSION_TYPE, &
    &                    SURF_AND_MORE_LOCAL_TYPE, &
    &                    SURF_AND_MORE_TYPE, &
    &                    FLUX_TYPE
    USE YOMCST, ONLY: RLVTT, RLSTT, RSIGMA, RCPD
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE ECEARTH, ONLY: ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX

    USE CPLNG2

    IMPLICIT NONE

    ! Arguments
    TYPE(TSURF),                   INTENT(IN) :: YDSURF
    TYPE(DIMENSION_TYPE),          INTENT(IN) :: KDIM
    TYPE(SURF_AND_MORE_LOCAL_TYPE),INTENT(IN) :: SURFL
    TYPE(SURF_AND_MORE_TYPE),      INTENT(IN) :: PSURF
    TYPE(FLUX_TYPE),               INTENT(IN) :: FLUX

    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG
    REAL(KIND=JPRB)    :: ZAHFLTI(KDIM%KLON,2)
    REAL(KIND=JPRB)    :: ZTS2(KDIM%KLON)
    REAL(KIND=JPRB)    :: ZTS3(KDIM%KLON)
    REAL(KIND=JPRB)    :: ZU10(KDIM%KLON)
    REAL(KIND=JPRB)    :: ZCALV(KDIM%KLON)


    IF (LHOOK) CALL DR_HOOK('ECE_NEMO_SET_OCEAN_FLUXES',0,ZHOOK_HANDLE)

    ASSOCIATE(YSD_VD => YDSURF%YSD_VD, YSP_SG=>YDSURF%YSP_SG, TSPHY=>YRPHY2%TSPHY)

    ! =========================================================================
    ! *** Pre-compute indices
    ! =========================================================================

    IL = KDIM%KIDIA
    IE = KDIM%KFDIA - KDIM%KIDIA
    IG = KDIM%KSTGLO - 1 + KDIM%KIDIA

    ! =========================================================================
    ! *** Momentum fluxes (stresses)
    ! =========================================================================

    CPLNG2_FLD(CPLNG2_IDX('A_TauX_oce'))%D(IG:IG+IE,1,1) = PSURF%PUSTRTI(IL:IL+IE,1)
    CPLNG2_FLD(CPLNG2_IDX('A_TauY_oce'))%D(IG:IG+IE,1,1) = PSURF%PVSTRTI(IL:IL+IE,1)
    CPLNG2_FLD(CPLNG2_IDX('A_TauX_ice'))%D(IG:IG+IE,1,1) = PSURF%PUSTRTI(IL:IL+IE,2)
    CPLNG2_FLD(CPLNG2_IDX('A_TauY_ice'))%D(IG:IG+IE,1,1) = PSURF%PVSTRTI(IL:IL+IE,2)

    ! =========================================================================
    ! *** Radiative fluxes (solar, non-solar, dQ/dT)
    ! =========================================================================

    IF (ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX) THEN
      CPLNG2_FLD(CPLNG2_IDX('A_Qs_mix'  ))%D(IG:IG+IE,1,1) = &
      &                     SURFL%ZFRTI(IL:IL+IE,1)*SURFL%ZFRSOTI(IL:IL+IE,1) &
      &                   + SURFL%ZFRTI(IL:IL+IE,2)*SURFL%ZFRSOTI(IL:IL+IE,2)
    ELSE
      CPLNG2_FLD(CPLNG2_IDX('A_Qs_oce'  ))%D(IG:IG+IE,1,1) = &
      &                                             SURFL%ZFRSOTI(IL:IL+IE,1)
    ENDIF

    CPLNG2_FLD(CPLNG2_IDX('A_Qs_ice'  ))%D(IG:IG+IE,1,1) = &
    &                                               SURFL%ZFRSOTI(IL:IL+IE,2)

    ! Latent heat flux is computed from evaporation over water(1) and ice(2)
    ZAHFLTI(IL:IL+IE,1) = PSURF%PEVAPTI(IL:IL+IE,1) * RLVTT
    ZAHFLTI(IL:IL+IE,2) = PSURF%PEVAPTI(IL:IL+IE,2) * RLSTT

    CPLNG2_FLD(CPLNG2_IDX('A_Qns_ice'))%D(IG:IG+IE,1,1) = &
    &                                             PSURF%PAHFSTI(IL:IL+IE,2) &
    &                                           + ZAHFLTI(IL:IL+IE,2)       &
    &                                           + SURFL%ZAHFTRTI(IL:IL+IE,2)

    IF (ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX) THEN
      CPLNG2_FLD(CPLNG2_IDX('A_Qns_mix'))%D(IG:IG+IE,1,1) = &
      &   SURFL%ZFRTI(IL:IL+IE,2) * CPLNG2_FLD(CPLNG2_IDX('A_Qns_ice'))%D(IG:IG+IE,1,1) &
      & + SURFL%ZFRTI(IL:IL+IE,1) * ( PSURF%PAHFSTI(IL:IL+IE,1)    &
      &                               + ZAHFLTI(IL:IL+IE,1)        &
      &                               + SURFL%ZAHFTRTI(IL:IL+IE,1) )
    ELSE
      CPLNG2_FLD(CPLNG2_IDX('A_Qns_oce'))%D(IG:IG+IE,1,1) = &
      &                                           PSURF%PAHFSTI(IL:IL+IE,1) &
      &                                         + ZAHFLTI(IL:IL+IE,1)       &
      &                                         + SURFL%ZAHFTRTI(IL:IL+IE,1)
    ENDIF

    ! Sensitivity of non-solar heat flux (only over ice)
    ZTS2(IL:IL+IE) = PSURF%PTSKTI(IL:IL+IE,2)**2
    ZTS3(IL:IL+IE) = PSURF%PTSKTI(IL:IL+IE,2)**3
    ZU10(IL:IL+IE) = SQRT(  PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y10U%MP)**2 &
    &                     + PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y10V%MP)**2 )

    ! From NEMO core bulk formulae
    ! Pay attention to the signs from the various contributions!
    CPLNG2_FLD(CPLNG2_IDX('A_dQns_dT'))%D(IG:IG+IE,1,1) =               &
    &                   -4.00 * 0.95 * RSIGMA * ZTS3(IL:IL+IE)        &
    &                   -1.22 * RCPD * 1.63e-3 * ZU10(IL:IL+IE)       &
    &                   + RLSTT * 1.63e-3 * 11637800.                 &
    &                     * (-5897.8) * ZU10(IL:IL+IE)/ZTS2(IL:IL+IE) &
    &                     * EXP(-5897.8/PSURF%PTSKTI(IL:IL+IE,2))

    ! =========================================================================
    ! *** Mass fluxes (runoff, precipitation, evaporation)
    ! =========================================================================

    CPLNG2_FLD(CPLNG2_IDX('A_Runoff'))%D(IG:IG+IE,1,1) = &
    &                             FLUX%PFWRO1(IL:IL+IE) + FLUX%PFWROD(IL:IL+IE)

    ! remove excess snow and send it into the ocean as ice ("calving")
    ! the threshold 10000 kg/m2 is "loosely" defined at l.638
    ! in surf/module/surftstp_ctl_mod.F90
    ! transform excess snow to a mass flux

    ! TODO We have multi layer snow now PSURF%PSP_SG is 3D instead of 2D. Logic needs to be redone!

    !ZCALV(IL:IL+IE)=MAX(0._JPRB,PSURF%PSP_SG(IL:IL+IE,YSP_SG%YF%MP)-10000._JPRB)/TSPHY
    CPLNG2_FLD(CPLNG2_IDX('A_Calving'))%D(IG:IG+IE,1,1) = 0 !ZCALV(IL:IL+IE)
    ! remove the same amount from the snow tendency to keep the mass balance
    !PSURF%PSNSE1(IL:IL+IE,1:KDIM%KLEVSN) = PSURF%PSNSE1(IL:IL+IE,1:KDIM%KLEVSN)-ZCALV(IL:IL+IE,1:KDIM%KLEVSN)

    CPLNG2_FLD(CPLNG2_IDX('A_Precip_liquid'))%D(IG:IG+IE,1,1) = &
    &         FLUX%PFPLCL(IL:IL+IE,KDIM%KLEV) + FLUX%PFPLSL(IL:IL+IE,KDIM%KLEV)

    CPLNG2_FLD(CPLNG2_IDX('A_Precip_solid'))%D(IG:IG+IE,1,1) = &
    &         FLUX%PFPLCN(IL:IL+IE,KDIM%KLEV) + FLUX%PFPLSN(IL:IL+IE,KDIM%KLEV)


    CPLNG2_FLD(CPLNG2_IDX('A_Evap_total'))%D(IG:IG+IE,1,1) = &
    &                     - PSURF%PEVAPTI(IL:IL+IE,1) * SURFL%ZFRTI(IL:IL+IE,1) &
    &                     - PSURF%PEVAPTI(IL:IL+IE,2) * SURFL%ZFRTI(IL:IL+IE,2)

    CPLNG2_FLD(CPLNG2_IDX('A_Evap_ice'))%D(IG:IG+IE,1,1) = &
    &                     - PSURF%PEVAPTI(IL:IL+IE,2)

    END ASSOCIATE

    IF (LHOOK) CALL DR_HOOK('ECE_NEMO_SET_OCEAN_FLUXES',1,ZHOOK_HANDLE)

END SUBROUTINE ECE_NEMO_SET_OCEAN_FLUXES
