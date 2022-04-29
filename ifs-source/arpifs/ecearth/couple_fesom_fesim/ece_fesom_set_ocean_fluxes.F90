SUBROUTINE ECE_FESOM_SET_OCEAN_FLUXES(YDGEOMETRY,YDSURF,KDIM,SURFL,PSURF,FLUX,PAUX)

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
    USE CPLNG
    USE YOETHF, ONLY: RHOH2O
    USE YOMRIP, ONLY: YRRIP

    IMPLICIT NONE

    ! Arguments
    TYPE(GEOMETRY),                INTENT(IN) :: YDGEOMETRY
    TYPE(TSURF),                   INTENT(INOUT) :: YDSURF
    TYPE(DIMENSION_TYPE),          INTENT(IN) :: KDIM
    TYPE(SURF_AND_MORE_LOCAL_TYPE),INTENT(IN) :: SURFL
    TYPE(SURF_AND_MORE_TYPE),      INTENT(IN) :: PSURF
    TYPE(FLUX_TYPE),               INTENT(IN) :: FLUX
    TYPE(AUX_TYPE),                INTENT(IN) :: PAUX

    ! Locals
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG
    REAL(KIND=JPRB)    :: ZMASK(KDIM%KLON)
    REAL(KIND=JPIM)    :: POSMASK(KDIM%KLON)

#include "update_fields.intfb.h"

    ASSOCIATE(YSD_VD=>YDSURF%YSD_VD, TSTEP=>YRRIP%TSTEP, PGELAT=>PAUX%PGELAT)

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
    CPLNG_FLD(CPLNG_IDX('A_TauX_oce'))%D(IG:IG+IE,1,1) = PSURF%PUSTRTI(IL:IL+IE,1)
    CPLNG_FLD(CPLNG_IDX('A_TauY_oce'))%D(IG:IG+IE,1,1) = PSURF%PVSTRTI(IL:IL+IE,1)
    CPLNG_FLD(CPLNG_IDX('A_TauX_ice'))%D(IG:IG+IE,1,1) = PSURF%PUSTRTI(IL:IL+IE,2)
    CPLNG_FLD(CPLNG_IDX('A_TauY_ice'))%D(IG:IG+IE,1,1) = PSURF%PVSTRTI(IL:IL+IE,2)

    ! =========================================================================
    ! *** Heat fluxes (radiative and latent)
    ! =========================================================================
    CPLNG_FLD(CPLNG_IDX('A_Qns_oce'))%D(IG:IG+IE,1,1) = &                          ! Heat over ocean =
       PSURF%PEVAPTI(IL:IL+IE,1) * RLVTT + &                                       ! Latent from evap +
       PSURF%PAHFSTI(IL:IL+IE,1) + &                                               ! Sensible at surface +
       SURFL%ZAHFTRTI(IL:IL+IE,1)                                                  ! Net lw rad at surface

    ZMASK(IL:IL+IE) = 0                                                            ! Initialze Binary sea ice mask
    POSMASK(IL:IL+IE) = 0                                                          ! Initialze sign mask for sea ice flux
    WHERE(SURFL%ZFRTI(IL:IL+IE,2) .NE. 0.0) ZMASK = 1                              ! where more than 0 sea ice -> mask = 1
    WHERE((PSURF%PEVAPTI(IL:IL+IE,2) * RLSTT + &                                   ! where heat-to-ice flux greater 0 -> mask = 1
       PSURF%PAHFSTI(IL:IL+IE,2) + &
       SURFL%ZAHFTRTI(IL:IL+IE,2) + &
       SURFL%ZFRSOTI(IL:IL+IE,2)) * &
       ZMASK .GT. 0.0) POSMASK = 1

    CPLNG_FLD(CPLNG_IDX('A_Q_ice'))%D(IG:IG+IE,1,1) = &                            ! Heat over ice =
       (PSURF%PEVAPTI(IL:IL+IE,2) * RLSTT + &                                      ! Latent from subl +
       PSURF%PAHFSTI(IL:IL+IE,2) + &                                               ! Sensible at surface +
       SURFL%ZAHFTRTI(IL:IL+IE,2) + &                                              ! Net lw rad at surface
       SURFL%ZFRSOTI(IL:IL+IE,2)) * &                                              ! Net sw rad at surface
       ZMASK                                                                       ! Send only where there is any sea ice

    CPLNG_FLD(CPLNG_IDX('A_Qs_all'))%D(IG:IG+IE,1,1) = &                           ! Sw down ocean and ice =
       SURFL%ZFRTI(IL:IL+IE,1)*SURFL%ZFRSOTI(IL:IL+IE,1)+ &
       SURFL%ZFRTI(IL:IL+IE,2)*SURFL%ZFRSOTI(IL:IL+IE,2)

    ! =========================================================================
    ! *** Mass fluxes (runoff+calving, precipitation, evaporation)
    ! =========================================================================
    CPLNG_FLD(CPLNG_IDX('A_Runoff'))%D(IG:IG+IE,1,1) = &
       (FLUX%PFWRO1(IL:IL+IE) + FLUX%PFWROD(IL:IL+IE))

    CPLNG_FLD(CPLNG_IDX('A_Precip_liquid'))%D(IG:IG+IE,1,1) = &
       FLUX%PFPLCL(IL:IL+IE,KDIM%KLEV) / RHOH2O + &
       FLUX%PFPLSL(IL:IL+IE,KDIM%KLEV) / RHOH2O

    CPLNG_FLD(CPLNG_IDX('A_Precip_solid'))%D(IG:IG+IE,1,1) = &
       FLUX%PFPLCN(IL:IL+IE,KDIM%KLEV) / RHOH2O + &
       FLUX%PFPLSN(IL:IL+IE,KDIM%KLEV) / RHOH2O

    CPLNG_FLD(CPLNG_IDX('A_Evap'))%D(IG:IG+IE,1,1) = &
       PSURF%PEVAPTI(IL:IL+IE,1) * SURFL%ZFRTI(IL:IL+IE,1) / RHOH2O

    CPLNG_FLD(CPLNG_IDX('A_Subl'))%D(IG:IG+IE,1,1) = &
       (PSURF%PEVAPTI(IL:IL+IE,2) * SURFL%ZFRTI(IL:IL+IE,2) / RHOH2O)

    IF (LHOOK) CALL DR_HOOK('ECE_FESOM_SET_OCEAN_FLUXES',1,ZHOOK_HANDLE)

    END ASSOCIATE

END SUBROUTINE ECE_FESOM_SET_OCEAN_FLUXES
