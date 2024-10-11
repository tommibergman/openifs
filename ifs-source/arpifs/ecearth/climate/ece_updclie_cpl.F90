SUBROUTINE ECE_UPDCLIE_CPL(YDGEOMETRY, YDSURF, PTSTEP)

    USE PARKIND1, ONLY: JPRB, JPIM
    USE GEOMETRY_MOD, ONLY: GEOMETRY
    USE SURFACE_FIELDS_MIX, ONLY: TSURF
    USE YOMMCC, ONLY: YRMCC
    USE YOEPHY, ONLY: YREPHY
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK

    USE ECEARTH
    USE CPLNG

    IMPLICIT NONE

#include "surf_inq.h"

    ! Arguments
    TYPE(GEOMETRY),  INTENT(IN)    :: YDGEOMETRY
    TYPE(TSURF),     INTENT(INOUT) :: YDSURF
    REAL(KIND=JPRB), INTENT(IN)    :: PTSTEP ! Time step

    ! Locals
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    INTEGER(KIND=JPIM) :: IEND, IBL
    INTEGER(KIND=JPIM) :: JSTGLO, JROF

    REAL(KIND=JPRB),POINTER :: CPL_FLD_SST(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_FRAC(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_TEMP(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_ALB(:)
    REAL(KIND=JPRB) :: ZRTFREEZSICE, ZRCIMIN
    REAL(KIND=JPRB) :: ZTS, ZCI, ZTI


    IF (LHOOK) CALL DR_HOOK('ECE_UPDCLIE_CPL',0,ZHOOK_HANDLE)

    ASSOCIATE(NPROMA => YDGEOMETRY%YRDIM%NPROMA, &
    &         NGPTOT => YDGEOMETRY%YRGEM%NGPTOT, &
    &         LEOCWA => YREPHY%LEOCWA, &
    &         LEOCCO => YREPHY%LEOCCO, &
    &         SP_SB => YDSURF%SP_SB, &
    &         YSP_SB => YDSURF%YSP_SB, &
    &         SP_RR => YDSURF%SP_RR, &
    &         YSP_RR => YDSURF%YSP_RR, &
    &         SD_VF => YDSURF%SD_VF, &
    &         YSD_VF => YDSURF%YSD_VF)

    ! =========================================================================
    ! *** 0. Initialisation
    ! =========================================================================

    CALL SURF_INQ(YREPHY%YSURF, PRTFREEZSICE=ZRTFREEZSICE, PRCIMIN=ZRCIMIN)

    ! =========================================================================
    ! *** 1. Update coupling fields (from CPLNG coupler)
    ! =========================================================================

    CALL CPLNG_EXCHANGE(KSTAGE=ECE_CPL_STAGE_OCE_RCV)

    CPL_FLD_SST => CPLNG_FLD(CPLNG_IDX('A_SST'))%D(:,1,1)
    CPL_FLD_ICE_FRAC => CPLNG_FLD(CPLNG_IDX('A_Ice_frac'))%D(:,1,1)
    IF (ECE_CPL_NEMO_LIM) THEN
      CPL_FLD_ICE_TEMP => CPLNG_FLD(CPLNG_IDX('A_Ice_temp'))%D(:,1,1)
      CPL_FLD_ICE_ALB => CPLNG_FLD(CPLNG_IDX('A_Ice_albedo'))%D(:,1,1)
    ENDIF

    ! =========================================================================
    ! *** 2. Update IFS variables from coupling
    ! =========================================================================

    DO JSTGLO=1,NGPTOT,NPROMA

      IEND = MIN(NPROMA,NGPTOT-JSTGLO+1)
      IBL = (JSTGLO-1)/NPROMA+1

      DO JROF = 1,IEND
        ! identify land/lake gridpoints similar to the logic in
        ! src/surf/module/surfbc_ctl_mod.F90
        IF((SD_VF(JROF,YSD_VF%YLSM%MP,IBL) <= .5_JPRB) .AND. &
        &  (SD_VF(JROF,YSD_VF%YCLK%MP,IBL) <= .5_JPRB)) THEN

          ! Sea-surface temperature
          SD_VF(JROF,YSD_VF%YSST%MP,IBL) = CPL_FLD_SST(JSTGLO+JROF-1)

          ! Sea-ice fraction and temperature
          ZCI = CPL_FLD_ICE_FRAC(JSTGLO+JROF-1)
          IF (ZCI > ZRCIMIN) THEN
            SD_VF(JROF,YSD_VF%YCI%MP,IBL) = ZCI
            IF (ECE_CPL_NEMO_LIM) THEN
              SP_SB(JROF,1,YSP_SB%YTL%MP,IBL) = CPL_FLD_ICE_TEMP(JSTGLO+JROF-1)
            ENDIF
          ELSE
            SD_VF(JROF,YSD_VF%YCI%MP,IBL) = 0.
            SP_SB(JROF,1,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
            SP_SB(JROF,2,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
            SP_SB(JROF,3,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
            SP_SB(JROF,4,YSP_SB%YTL%MP,IBL) = ZRTFREEZSICE
          ENDIF

          ! Surface and soil temperature
          SP_SB(JROF,1,YSP_SB%YT%MP,IBL) = SD_VF(JROF,YSD_VF%YCI%MP,IBL)* &
          &  SP_SB(JROF,1,YSP_SB%YTL%MP,IBL)+(1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))* &
          &  SD_VF(JROF,YSD_VF%YSST%MP,IBL)
          SP_SB(JROF,2,YSP_SB%YT%MP,IBL) = SD_VF(JROF,YSD_VF%YCI%MP,IBL)* &
          &  SP_SB(JROF,2,YSP_SB%YTL%MP,IBL)+(1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))* &
          &  SD_VF(JROF,YSD_VF%YSST%MP,IBL)
          SP_SB(JROF,3,YSP_SB%YT%MP,IBL) = SD_VF(JROF,YSD_VF%YCI%MP,IBL)* &
          &  SP_SB(JROF,3,YSP_SB%YTL%MP,IBL)+(1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))* &
          &  SD_VF(JROF,YSD_VF%YSST%MP,IBL)
          SP_SB(JROF,4,YSP_SB%YT%MP,IBL) = SD_VF(JROF,YSD_VF%YCI%MP,IBL)* &
          &  SP_SB(JROF,4,YSP_SB%YTL%MP,IBL)+(1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))* &
          &  SD_VF(JROF,YSD_VF%YSST%MP,IBL)

          ! Skin temperature
          SP_RR(JROF,YSP_RR%YT%MP,IBL) = &
          & ( SD_VF(JROF,YSD_VF%YCI%MP,IBL)*SP_SB(JROF,1,YSP_SB%YTL%MP,IBL)**4 &
          & + (1.-SD_VF(JROF,YSD_VF%YCI%MP,IBL))*SD_VF(JROF,YSD_VF%YSST%MP,IBL)**4 &
          & )**.25

          ! Sea-ice albedo
          IF (ECE_CPL_NEMO_LIM) THEN
            SD_VF(JROF,YSD_VF%YALBF%MP,IBL) = CPL_FLD_ICE_ALB(JSTGLO+JROF-1)
          ENDIF

        ENDIF ! LSM <= 0.5_JPRB
      ENDDO ! JROF = 1,IEND
    ENDDO ! JSTGLO=1,NGPTOT,NPROMA

    END ASSOCIATE

    IF (LHOOK) CALL DR_HOOK('ECE_UPDCLIE_CPL',1,ZHOOK_HANDLE)

END SUBROUTINE ECE_UPDCLIE_CPL
