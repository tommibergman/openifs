SUBROUTINE ECE_FESOM_SET_OCEAN_FLUXES(YDGEOMETRY,YDSURF,KDIM,SURFL,PSURF,FLUX,PAUX,TSTEP)

#ifdef WITH_CPLNG2

    USE GEOMETRY_MOD, ONLY: GEOMETRY

    USE PARKIND1, ONLY: JPRB, JPIM
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE YOMPHYDER, ONLY: DIMENSION_TYPE, &
                         SURF_AND_MORE_LOCAL_TYPE, &
                         SURF_AND_MORE_TYPE, &
                         FLUX_TYPE, &
                         AUX_TYPE
    USE SURFACE_FIELDS_MIX, ONLY : TSURF
    USE YOMCST, ONLY: RLVTT, RLSTT, RSIGMA, RCPD
    USE CPLNG2
    USE YOETHF, ONLY: RHOH2O

    IMPLICIT NONE

    ! Arguments
    TYPE(GEOMETRY),                INTENT(IN) :: YDGEOMETRY
    TYPE(TSURF),                   INTENT(INOUT) :: YDSURF
    TYPE(DIMENSION_TYPE),          INTENT(IN) :: KDIM
    TYPE(SURF_AND_MORE_LOCAL_TYPE),INTENT(IN) :: SURFL
    TYPE(SURF_AND_MORE_TYPE),      INTENT(INOUT) :: PSURF
    TYPE(FLUX_TYPE),               INTENT(IN) :: FLUX
    TYPE(AUX_TYPE),                INTENT(IN) :: PAUX
    REAL(KIND=JPRB),               INTENT(IN) :: TSTEP

    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG,IH
    REAL(KIND=JPRB)    :: ZMASK(KDIM%KLON)
    REAL(KIND=JPIM)    :: POSMASK(KDIM%KLON)
    REAL(KIND=JPRB)    :: ZEXSNOW(KDIM%KLON)
    REAL(KIND=JPRB)    :: HEATPREC(KDIM%KLON)
    REAL(KIND=JPRB)    :: SST(KDIM%KLON)
    REAL(KIND=JPRB)    :: T2M(KDIM%KLON)
    REAL(KIND=JPRB)    :: RAIN(KDIM%KLON)
    REAL(KIND=JPRB)    :: SNOW(KDIM%KLON)
    REAL(KIND=JPRB)    :: DELTA_TEMP(KDIM%KLON)

#include "update_fields.intfb.h"

    ASSOCIATE(YSD_VD=>YDSURF%YSD_VD, YSD_VF=>YDSURF%YSD_VF, YSP_SG=>YDSURF%YSP_SG, PGELAT=>PAUX%PGELAT)

    IF (LHOOK) CALL DR_HOOK('ECE_FESOM_SET_OCEAN_FLUXES',0,ZHOOK_HANDLE)

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
    ! *** Heat fluxes (radiative and latent)
    ! =========================================================================
    SST(IL:IL+IE) = PSURF%PSD_VF(IL:IL+IE,YSD_VF%YSST%MP)
    T2M(IL:IL+IE) = PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y2T%MP)
    RAIN(IL:IL+IE) = FLUX%PFPLCL(IL:IL+IE,KDIM%KLEV) / RHOH2O + FLUX%PFPLSL(IL:IL+IE,KDIM%KLEV) / RHOH2O
    SNOW(IL:IL+IE) = FLUX%PFPLCN(IL:IL+IE,KDIM%KLEV) / RHOH2O + FLUX%PFPLSN(IL:IL+IE,KDIM%KLEV) / RHOH2O
    DO IH = IL, IL+IE
      IF (ABS(T2M(IH))- ABS(SST(IH)) < 0.1e-15) THEN
        DELTA_TEMP(IH)=0
      ELSE
        DELTA_TEMP(IH) = T2M(IH) - SST(IH)
      ENDIF
      HEATPREC = (DELTA_TEMP(IH)*RAIN(IH)*4.186_JPRB)-((DELTA_TEMP(IH))*SNOW(IH)*(4.186-2.108))*1000000
    ENDDO

    CPLNG2_FLD(CPLNG2_IDX('A_Qns_oce'))%D(IG:IG+IE,1,1) = &                        ! Heat over ocean =
       PSURF%PEVAPTI(IL:IL+IE,1) * RLVTT + &                                       ! Latent from evap +
       PSURF%PAHFSTI(IL:IL+IE,1) + &                                               ! Sensible at surface +
       SURFL%ZAHFTRTI(IL:IL+IE,1) - &                                              ! Net lw rad at surface
       (FLUX%PFPLCN(IL:IL+IE,KDIM%KLEV) + &                                        ! (Enthalpy of fusion of snow with
       FLUX%PFPLSN(IL:IL+IE,KDIM%KLEV)) * 333.5_JPRB * 1000.0_JPRB + &             ! conversion from m depth to gramm)
       HEATPREC(IL:IL+IE)                                                          ! Heatflux from precipitation temperature

    ZMASK(IL:IL+IE) = 0                                                            ! Initialze Binary sea ice mask
    WHERE(SURFL%ZFRTI(IL:IL+IE,2) .GT. 0.1e-15) ZMASK(IL:IL+IE) = 1                    ! where more than 0 sea ice -> mask = 1

    CPLNG2_FLD(CPLNG2_IDX('A_Q_ice'))%D(IG:IG+IE,1,1) = &                          ! Heat over ice =
       (PSURF%PEVAPTI(IL:IL+IE,2) * RLSTT + &                                      ! Latent from subl +
       PSURF%PAHFSTI(IL:IL+IE,2) + &                                               ! Sensible at surface +
       SURFL%ZAHFTRTI(IL:IL+IE,2) + &                                              ! Net lw rad at surface
       SURFL%ZFRSOTI(IL:IL+IE,2)) * &                                              ! Net sw rad at surface
       ZMASK(IL:IL+IE)                                                             ! Send only where there is any sea ice

    CPLNG2_FLD(CPLNG2_IDX('A_Qs_all'))%D(IG:IG+IE,1,1) = &                         ! Separatly send SW down ocean for SW-penetration
       SURFL%ZFRSOTI(IL:IL+IE,1)

    ! =========================================================================
    ! *** Mass fluxes (runoff+calving, precipitation, evaporation)
    ! =========================================================================
    ZEXSNOW(IL:IL+IE) = MAX(SUM(PSURF%PSP_SG(IL:IL+IE,1:KDIM%KLEVSN,YSP_SG%YF%MP9),DIM=2)+&  ! Snow water mass [kg/m**2]
       SUM(PSURF%PSNSE1(IL:IL+IE,1:KDIM%KLEVSN),DIM=2)*TSTEP-10000.0_JPRB,0.)                ! Snow mass tendency [kg/sm**2] - perannial snow hight [kg/m**2]

    PSURF%PSNSE1(IL:IL+IE,1) = PSURF%PSNSE1(IL:IL+IE,1)-&
                             ZEXSNOW(IL:IL+IE)/TSTEP

    CPLNG2_FLD(CPLNG2_IDX('A_Calving'))%D(IG:IG+IE,1,1) = &
       (ZEXSNOW(IL:IL+IE)/TSTEP)/1000

    CPLNG2_FLD(CPLNG2_IDX('A_Runoff'))%D(IG:IG+IE,1,1) = &
       (FLUX%PFWRO1(IL:IL+IE) + FLUX%PFWROD(IL:IL+IE))/1000

    CPLNG2_FLD(CPLNG2_IDX('A_Precip_liquid'))%D(IG:IG+IE,1,1) = &
       FLUX%PFPLCL(IL:IL+IE,KDIM%KLEV) / RHOH2O + &
       FLUX%PFPLSL(IL:IL+IE,KDIM%KLEV) / RHOH2O

    CPLNG2_FLD(CPLNG2_IDX('A_Precip_solid'))%D(IG:IG+IE,1,1) = &
       FLUX%PFPLCN(IL:IL+IE,KDIM%KLEV) / RHOH2O + &
       FLUX%PFPLSN(IL:IL+IE,KDIM%KLEV) / RHOH2O

    CPLNG2_FLD(CPLNG2_IDX('A_Evap'))%D(IG:IG+IE,1,1) = &
       PSURF%PEVAPTI(IL:IL+IE,1) / RHOH2O

    CPLNG2_FLD(CPLNG2_IDX('A_Subl'))%D(IG:IG+IE,1,1) = &
       (PSURF%PEVAPTI(IL:IL+IE,2) / RHOH2O)

    CPLNG2_FLD(CPLNG2_IDX('A_WindX'))%D(IG:IG+IE,1,1) = &
       (PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y10U%MP))

    CPLNG2_FLD(CPLNG2_IDX('A_WindY'))%D(IG:IG+IE,1,1) = &
       (PSURF%PSD_VD(IL:IL+IE,YSD_VD%Y10V%MP))

    IF (LHOOK) CALL DR_HOOK('ECE_FESOM_SET_OCEAN_FLUXES',1,ZHOOK_HANDLE)

    END ASSOCIATE

#endif

END SUBROUTINE ECE_FESOM_SET_OCEAN_FLUXES
