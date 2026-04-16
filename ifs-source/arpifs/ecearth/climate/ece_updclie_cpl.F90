SUBROUTINE ECE_UPDCLIE_CPL(YDGEOMETRY, YDSURF, YDMCC, YDDYNA, YDRIP, PTSTEP)

#ifdef WITH_CPLNG2

    USE PARKIND1, ONLY: JPRB, JPIM
    USE GEOMETRY_MOD, ONLY: GEOMETRY
    USE SURFACE_FIELDS_MIX, ONLY: TSURF
    USE YOEPHY, ONLY: YREPHY
    USE YOMMCC, ONLY : TMCC
    USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
    USE YOMDYNA, ONLY: TDYNA
    USE YOMRIP, ONLY : TRIP

    USE ECEARTH
    USE CPLNG2

    IMPLICIT NONE

#include "surf_inq.h"

    ! Arguments
    TYPE(GEOMETRY),  INTENT(IN)    :: YDGEOMETRY
    TYPE(TSURF),     INTENT(INOUT) :: YDSURF
    TYPE(TMCC) ,     INTENT(IN)    :: YDMCC
    TYPE(TDYNA),     INTENT(IN)    :: YDDYNA
    TYPE(TRIP),      INTENT(IN)    :: YDRIP
    REAL(KIND=JPRB), INTENT(IN)    :: PTSTEP ! Time step

    ! Locals
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    INTEGER(KIND=JPIM) :: IEND, IBL
    INTEGER(KIND=JPIM) :: JSTGLO, JROF

    REAL(KIND=JPRB),POINTER :: CPL_FLD_SST(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_FRAC(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_ICE_TEMP(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_OUCURR(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_OVCURR(:)
    REAL(KIND=JPRB) :: ZRTFREEZSICE, ZRCIMIN, ZPRTMELTSICE
    REAL(KIND=JPRB) :: ZTS, ZCI, ZTI

    ! LPJ-GUESS variables and structures
    REAL(KIND=JPRB) :: ZLAIL, ZLAIH, ZCVL, ZCVH, ZTVL, ZTVH, COVERSUM
    REAL(KIND=JPRB),POINTER :: CPL_FLD_VEG_LAIL(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_VEG_LAIH(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_VEG_CVL(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_VEG_CVH(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_VEG_TVL(:)
    REAL(KIND=JPRB),POINTER :: CPL_FLD_VEG_TVH(:)
    INTEGER :: dbg_i, dbg_n, dbg_n_LAIH, dbg_n_LAIL, dbg_n_CVL, dbg_n_CVH

    IF (LHOOK) CALL DR_HOOK('ECE_UPDCLIE_CPL',0,ZHOOK_HANDLE)

    ASSOCIATE(NPROMA => YDGEOMETRY%YRDIM%NPROMA, &
    &         NGPTOT => YDGEOMETRY%YRGEM%NGPTOT, &
    &         LEOCWA => YREPHY%LEOCWA, &
    &         LEOCCO => YREPHY%LEOCCO, &
    &         LNEMOLIMCUR => YDMCC%LNEMOLIMCUR, &
    &         LNEMOLIMTEMP => YDMCC%LNEMOLIMTEMP, &
    &         SP_SB => YDSURF%SP_SB, &
    &         YSP_SB => YDSURF%YSP_SB, &
    &         SP_SL => YDSURF%SP_SL, &
    &         YSP_SL => YDSURF%YSP_SL, &
    &         SP_RR => YDSURF%SP_RR, &
    &         YSP_RR => YDSURF%YSP_RR, &
    &         SD_VF => YDSURF%SD_VF, &
    &         YSD_VF => YDSURF%YSD_VF)

    ! =========================================================================
    ! *** 0. Initialisation
    ! =========================================================================

    CALL SURF_INQ(YREPHY%YSURF, PRTFREEZSICE=ZRTFREEZSICE, PRCIMIN=ZRCIMIN, &
    &             PRTMELTSICE=ZPRTMELTSICE )

    ! =========================================================================
    ! *** 1. Update coupling fields (from CPLNG coupler)
    ! =========================================================================

    CALL CPLNG2_EXCHANGE(INT(PTSTEP,KIND=JPIM),ECE_CPL_STAGE_OCE_RCV,YDDYNA,YDRIP)

    CPL_FLD_SST => CPLNG2_FLD(CPLNG2_IDX('A_SST'))%D(:,1,1)
    CPL_FLD_ICE_FRAC => CPLNG2_FLD(CPLNG2_IDX('A_Ice_frac'))%D(:,1,1)
    IF (LNEMOLIMTEMP) THEN
      CPL_FLD_ICE_TEMP => CPLNG2_FLD(CPLNG2_IDX('A_Ice_temp'))%D(:,1,1)
    ENDIF

    ! LPJ-GUESS coupling. Always receive fields when LPJ-GUESS is coupled
    CALL CPLNG2_EXCHANGE(INT(PTSTEP,KIND=JPIM),ECE_CPL_STAGE_VEG_RCV,YDDYNA,YDRIP)

    IF (ECE_CPL_LPJG) THEN
      CPL_FLD_VEG_LAIL => CPLNG2_FLD(CPLNG2_IDX('LAILVeg'))%D(:,1,1)
      CPL_FLD_VEG_LAIH => CPLNG2_FLD(CPLNG2_IDX('LAIHVeg'))%D(:,1,1)
      CPL_FLD_VEG_CVL => CPLNG2_FLD(CPLNG2_IDX('FracLVeg'))%D(:,1,1)
      CPL_FLD_VEG_CVH => CPLNG2_FLD(CPLNG2_IDX('FracHVeg'))%D(:,1,1)
      CPL_FLD_VEG_TVL => CPLNG2_FLD(CPLNG2_IDX('TypeLVeg'))%D(:,1,1)
      CPL_FLD_VEG_TVH => CPLNG2_FLD(CPLNG2_IDX('TypeHVeg'))%D(:,1,1)
    ENDIF

    IF (LNEMOLIMCUR) THEN
      CPL_FLD_OUCURR  => CPLNG2_FLD(CPLNG2_IDX('A_CurX'))%D(:,1,1)
      CPL_FLD_OVCURR  => CPLNG2_FLD(CPLNG2_IDX('A_CurY'))%D(:,1,1)
    ENDIF
    ! =========================================================================
    ! *** 2. Update IFS variables from coupling
    ! =========================================================================

    DO JSTGLO=1,NGPTOT,NPROMA

      IEND = MIN(NPROMA,NGPTOT-JSTGLO+1)
      IBL = (JSTGLO-1)/NPROMA+1

      DO JROF = 1,IEND
        ! Place ocean current velocity into IFS datastructure
        IF (LNEMOLIMCUR) THEN
          SD_VF(JROF,YSD_VF%YUCUR%MP,IBL) = CPL_FLD_OUCURR(JSTGLO+JROF-1)
          SD_VF(JROF,YSD_VF%YVCUR%MP,IBL) = CPL_FLD_OVCURR(JSTGLO+JROF-1)
        ENDIF
        ! identify land/lake gridpoints similar to the logic in
        ! src/surf/module/surfbc_ctl_mod.F90
        IF((SD_VF(JROF,YSD_VF%YLSM%MP,IBL) <= .5_JPRB) .AND. &
        &  (SD_VF(JROF,YSD_VF%YCLK%MP,IBL) <= .5_JPRB)) THEN

          ! Sea-surface temperature
          SD_VF(JROF,YSD_VF%YSST%MP,IBL) = CPL_FLD_SST(JSTGLO+JROF-1)

          ! Sea-ice fraction, temperature and albedo
          ZCI = CPL_FLD_ICE_FRAC(JSTGLO+JROF-1)
          IF (ZCI > ZRCIMIN) THEN
            SD_VF(JROF,YSD_VF%YCI%MP,IBL) = ZCI
            IF (LNEMOLIMTEMP) THEN
              IF (ECE_CPL_NEMO_WEIGHTED_ICE) THEN
                ! ensure ice fraction is above zero in JPRB precision
                ! Otherwise we could get strange numbers here and trigger unexpected behaviour
                IF (ZCI > EPSILON(1._JPRB)) THEN
                  SP_SB(JROF,1,YSP_SB%YTL%MP,IBL) = MIN(ZPRTMELTSICE,CPL_FLD_ICE_TEMP(JSTGLO+JROF-1)/ZCI)
                ELSE
                  ! If ZCI is too small, set ice temperature to melting point
                  ! Note this is ice fraction far below 0.01, so it should not have a big impact 
                  SP_SB(JROF,1,YSP_SB%YTL%MP,IBL) = ZPRTMELTSICE 
                ENDIF
              ELSE
                SP_SB(JROF,1,YSP_SB%YTL%MP,IBL) = CPL_FLD_ICE_TEMP(JSTGLO+JROF-1)
              ENDIF
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

        ELSE
          ! set SST as lake mixed layer temperature over non-ocean gridpoints
          ! this is ncessary because the open water fraction of lakes is represented
          ! by tile 1 in src/surf/module/surfbc_ctl_mod.F90
          SD_VF(JROF,YSD_VF%YSST%MP,IBL) = SP_SL(JROF,YSP_SL%YLMLT%MP,IBL)

          IF (ECE_CPL_LPJG) THEN

            ! Use LPJ-GUESS fields for land points that are NOT lakes.
            ! Lake values are set to 0 in LPJG

            ! LOW LAI
            ZLAIL=MAX(0.0_JPRB,CPL_FLD_VEG_LAIL(JSTGLO+JROF-1))
            SD_VF(JROF,YSD_VF%YLAIL%MP,IBL)=ZLAIL

            ! HIGH LAI
            ZLAIH=MAX(0.0_JPRB,CPL_FLD_VEG_LAIH(JSTGLO+JROF-1))
            SD_VF(JROF,YSD_VF%YLAIH%MP,IBL)=ZLAIH

            ! LOW COVER FRACTION
            ZCVL=MAX(0.0_JPRB,MIN(1.0_JPRB,CPL_FLD_VEG_CVL(JSTGLO+JROF-1)))

            ! HIGH COVER FRACTION
            ZCVH=MAX(0.0_JPRB,MIN(1.0_JPRB,CPL_FLD_VEG_CVH(JSTGLO+JROF-1)))

            COVERSUM = ZCVL + ZCVH
            ! Rescale if necessary to catch rounding errors leading to cover
            ! fractions > 1
            IF (COVERSUM > 1._JPRB) THEN
              ZCVL = ZCVL / COVERSUM
              ZCVH = 1._JPRB - ZCVL ! Paul M / Lars N - better rescaling
            ENDIF

            ! UPDATE COVER FRACTIONS
            SD_VF(JROF,YSD_VF%YCVL%MP,IBL)=ZCVL
            SD_VF(JROF,YSD_VF%YCVH%MP,IBL)=ZCVH

            ! LOW COVER TYPE. Restrict to be between 0 and 20
            ZTVL=CPL_FLD_VEG_TVL(JSTGLO+JROF-1)
            IF (ZTVL==8.0_JPRB) THEN
              ZTVL=0.0_JPRB ! Force deserts to have ZTVL=0
            ENDIF
            SD_VF(JROF,YSD_VF%YTVL%MP,IBL)=MIN(MAX(ZTVL,0.0_JPRB),20.0_JPRB)

            ! HIGH COVER TYPE. Restrict to be between 0 and 20
            ZTVH=CPL_FLD_VEG_TVH(JSTGLO+JROF-1)
            SD_VF(JROF,YSD_VF%YTVH%MP,IBL)=MIN(MAX(ZTVH,0.0_JPRB),20.0_JPRB)

            ! If we need to check the values of the fields, uncomment:
            ! WRITE (NULOUT, *) 'PaulM - LPJG coupling ', ZLAIL, ZLAIH, ZCVL, ZCVH, ZTVL, ZTVH

          ENDIF ! Vegetation fields

        ENDIF ! LSM <= 0.5_JPRB
      ENDDO ! JROF = 1,IEND
    ENDDO ! JSTGLO=1,NGPTOT,NPROMA

    END ASSOCIATE

    IF (LHOOK) CALL DR_HOOK('ECE_UPDCLIE_CPL',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ECE_UPDCLIE_CPL
