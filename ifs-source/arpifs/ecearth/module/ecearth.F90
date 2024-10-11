MODULE ECEARTH

    USE PARKIND1,ONLY : JPIM

    IMPLICIT NONE

    PRIVATE

    PUBLIC LECEARTH

    PUBLIC ECE_CONFIG
    PUBLIC ECE_FINALIZE

    PUBLIC ECE_CPL_STAGE_OCE_SND
    PUBLIC ECE_CPL_STAGE_OCE_RCV
    PUBLIC ECE_CPL_STAGE_CHE_SND
    PUBLIC ECE_CPL_STAGE_CHE_RCV
    PUBLIC ECE_CPL_STAGE_VEG_SND
    PUBLIC ECE_CPL_STAGE_VEG_RCV

    PUBLIC ECE_CPL_AMIP
    PUBLIC ECE_CPL_NEMO_LIM
    PUBLIC ECE_CPL_FESOM_FESIM

    PUBLIC ECE_CLIMR

    LOGICAL :: LECEARTH = .TRUE.  ! Main EC-Earth flag, always true

    LOGICAL :: ECE_CPL_AMIP = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_LIM = .FALSE.
    LOGICAL :: ECE_CPL_FESOM_FESIM = .FALSE.

    NAMELIST /NAMECECFG/ ECE_CPL_AMIP
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_LIM
    NAMELIST /NAMECECFG/ ECE_CPL_FESOM_FESIM

    LOGICAL :: ECE_CLIMR = .FALSE.

    INTEGER(KIND=JPIM),PARAMETER :: ECE_CPL_STAGE_OCE_SND = 1
    INTEGER(KIND=JPIM),PARAMETER :: ECE_CPL_STAGE_OCE_RCV = 2
    INTEGER(KIND=JPIM),PARAMETER :: ECE_CPL_STAGE_CHE_SND = 3
    INTEGER(KIND=JPIM),PARAMETER :: ECE_CPL_STAGE_CHE_RCV = 4
    INTEGER(KIND=JPIM),PARAMETER :: ECE_CPL_STAGE_VEG_SND = 5
    INTEGER(KIND=JPIM),PARAMETER :: ECE_CPL_STAGE_VEG_RCV = 6

CONTAINS

! =============================================================================
! *** ECE_CONFIG
! =============================================================================
SUBROUTINE ECE_CONFIG(YDGEOMETRY)

    USE GEOMETRY_MOD, ONLY: GEOMETRY
    USE YOMLUN, ONLY: NULOUT, NULNAM

    ! Argument
    TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY

    ! Read EC-Earth configuration namelist
    CALL POSNAM(NULNAM,'NAMECECFG')
    READ(NULNAM,NAMECECFG)
    WRITE(NULOUT, nml=NAMECECFG)

    ! Sanity check of coupling configuration
    IF (COUNT((/ECE_CPL_AMIP, ECE_CPL_NEMO_LIM, ECE_CPL_FESOM_FESIM/)) /= 1) THEN
        WRITE(NULOUT, *) &
        'Exactly one of ECE_CPL_AMIP, ECE_CPL_NEMO_LIM, ECE_CPL_FESOM_FESIM must be &
        &set to true in namelist NAMECECFG!'
        CALL ABOR1('ECE_CONFIG: ABOR1 CALLED (invalid EC-Earth configuration)')
    ENDIF

    ! Configure CPLNG
    IF (ECE_CPL_AMIP .OR. ECE_CPL_NEMO_LIM .OR. ECE_CPL_FESOM_FESIM) THEN
        CALL ECE_CONFIG_COUPLING(YDGEOMETRY)
    ENDIF

    ! Configure reading of CLIMR (ICMCL) files
    ECE_CLIMR = .TRUE.

END SUBROUTINE ECE_CONFIG

! =============================================================================
! *** ECE_CONFIG_COUPLING
! =============================================================================
SUBROUTINE ECE_CONFIG_COUPLING(YDGEOMETRY)

    USE GEOMETRY_MOD, ONLY: GEOMETRY
    USE YOMLUN, ONLY: NULOUT
    USE YOMMCC, ONLY: YRMCC

    USE CPLNG
    USE MOD_OASIS

    ! Argument
    TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY

    ASSOCIATE(LNEMOLIMALB => YRMCC%LNEMOLIMALB, &
    &         LNEMOLIMTEMP => YRMCC%LNEMOLIMTEMP, &
    &         LNEMOLIMTHK => YRMCC%LNEMOLIMTHK, &
    &         LNEMOCOUP => YRMCC%LNEMOCOUP, &
    &         LNEMO1WAY => YRMCC%LNEMO1WAY, &
    &         LNEMOFLUXNC => YRMCC%LNEMOFLUXNC, &
    &         LMCCDYNSEAICE => YRMCC%LMCCDYNSEAICE, &
    &         LNEMOLIMGET => YRMCC%LNEMOLIMGET, &
    &         LNEMOLIMPUT => YRMCC%LNEMOLIMPUT)

    ! -------------------------------------------------------------------------
    ! * (1) Configure coupling fields
    ! -------------------------------------------------------------------------

    ! 1.0 SST and sea-ice fraction: Needed in all setups
    CALL CPLNG_ADD_FLD('A_SST',     CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN,ECE_CPL_STAGE_OCE_RCV)
    CALL CPLNG_ADD_FLD('A_Ice_frac',CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN,ECE_CPL_STAGE_OCE_RCV)

    ! 1.1a IFS-NEMO coupling
    IF (ECE_CPL_NEMO_LIM) THEN

        LNEMOLIMALB  = .TRUE.
        LNEMOLIMTEMP = .TRUE.
        LNEMOLIMTHK  = .TRUE.

        WRITE(UNIT=NULOUT,FMT='(" Resetting some YOMCC variables to configure EC-Earth coupling with NEMO")')
        WRITE(UNIT=NULOUT,FMT='(&
        &  " LNEMOCOUP = ",L2, &
        &  " LNEMO1WAY = ",L2," LNEMOFLUXNC = ",L2," LMCCDYNSEAICE = ",L2, &
        &  " LNEMOLIMGET = ",L2," LNEMOLIMPUT = ",L2, &
        &  " LNEMOLIMALB = ",L2," LNEMOLIMTEMP = ",L2," LNEMOLIMTHK = ",L2)') &
        &  LNEMOCOUP, LNEMO1WAY, LNEMOFLUXNC, LMCCDYNSEAICE, &
        &  LNEMOLIMGET, LNEMOLIMPUT,LNEMOLIMALB, LNEMOLIMTEMP,LNEMOLIMTHK

        ! Fields sent from atmosphere to ocean
        CALL CPLNG_ADD_FLD('A_TauX_oce',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_TauY_oce',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_TauX_ice',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_TauY_ice',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Qs_mix',        CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Qs_ice',        CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Qns_mix',       CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Qns_ice',       CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Precip_liquid', CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Precip_solid',  CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Evap_total',    CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Evap_ice',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_dQns_dT',       CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)

        ! Fields sent from atmosphere to ocean via runoff-mapper
        CALL CPLNG_ADD_FLD('A_Runoff',        CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)

        ! Fields received by the atmosphere from the ocean
        CALL CPLNG_ADD_FLD('A_Ice_temp',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)
        CALL CPLNG_ADD_FLD('A_Ice_albedo',    CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)
        CALL CPLNG_ADD_FLD('A_Ice_thickness', CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)
        CALL CPLNG_ADD_FLD('A_Snow_thickness',CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)

    ENDIF ! ECE_CPL_NEMO_LIM

    
    ! 1.1b IFS-FESOM coupling
    IF (ECE_CPL_FESOM_FESIM) THEN

        LNEMOLIMALB  = .TRUE.
        LNEMOLIMTEMP = .TRUE.
        LNEMOLIMTHK  = .TRUE.

          WRITE(UNIT=NULOUT,FMT='(" Resetting some YOMCC variables to configure OIFS-FESOM2 coupling")')
          WRITE(UNIT=NULOUT,FMT='(&
          &  " LNEMOCOUP = ",L2, &
          &  " LNEMO1WAY = ",L2," LNEMOFLUXNC = ",L2," LMCCDYNSEAICE = ",L2, &
          &  " LNEMOLIMGET = ",L2," LNEMOLIMPUT = ",L2, &
          &  " LNEMOLIMALB = ",L2," LNEMOLIMTEMP = ",L2," LNEMOLIMTHK = ",L2)')&
          &  LNEMOCOUP, LNEMO1WAY, LNEMOFLUXNC, LMCCDYNSEAICE, &
          &  LNEMOLIMGET, LNEMOLIMPUT,LNEMOLIMALB, LNEMOLIMTEMP,LNEMOLIMTHK

        ! Fields sent from atmosphere to ocean
        CALL CPLNG_ADD_FLD('A_TauX_oce',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_TauY_oce',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_TauX_ice',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_TauY_ice',      CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Q_ice',         CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Qns_oce',       CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Qs_all',        CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Precip_liquid', CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Precip_solid',  CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Evap',          CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
        CALL CPLNG_ADD_FLD('A_Subl',          CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)

        ! Fields sent from atmosphere to ocean via runoff-mapper
        CALL CPLNG_ADD_FLD('A_Runoff',        CPLNG_FLD_TYPE_GRIDPOINT,OASIS_OUT,ECE_CPL_STAGE_OCE_SND)
                                                                                                      
        ! Fields received by the atmosphere from the ocean                                            
        CALL CPLNG_ADD_FLD('A_Ice_temp'      ,CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)
        CALL CPLNG_ADD_FLD('A_Ice_albedo'    ,CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)
        CALL CPLNG_ADD_FLD('A_Snow_thickness',CPLNG_FLD_TYPE_GRIDPOINT,OASIS_IN, ECE_CPL_STAGE_OCE_RCV)

    ENDIF ! ECE_CPL_FESOM_FESIM


    ! -------------------------------------------------------------------------
    ! * (2) Complete CPLNG configuration
    ! -------------------------------------------------------------------------
    CALL CPLNG_ADD_FLD_COMPLETED(YDGEOMETRY)

    END ASSOCIATE

END SUBROUTINE ECE_CONFIG_COUPLING

! =============================================================================
! *** ECE_FINALIZE
! =============================================================================
SUBROUTINE ECE_FINALIZE

    USE CPLNG

    ! -------------------------------------------------------------------------
    ! * (1) Shut down CPLNG
    ! -------------------------------------------------------------------------
    CALL CPLNG_FINALIZE

END SUBROUTINE ECE_FINALIZE

END MODULE ECEARTH
