MODULE SURFECE
    ! Some logic switches for EC-Earth need to be available in surf.
    ! As the main module ECEARTH is part of the later compiled arpifs,
    ! which depends on surf, we create this new module.

    ! Authors:
    ! 2023-11-22: Jan Streffing (AWI): Created separate module from ECEARTH

    USE PARKIND1,ONLY : JPIM

    IMPLICIT NONE

    PRIVATE

    PUBLIC LECEARTH

    PUBLIC ECE_CPL_AMIP
    PUBLIC ECE_CPL_NEMO_LIM
    PUBLIC ECE_CPL_FESOM_FESIM

    PUBLIC ECE_CLIMR

    PUBLIC ECE_CPL_NEMO_WEIGHTED_ICE
    PUBLIC ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX
    PUBLIC ECE_CPL_LPJG
    PUBLIC ECE_CPL_ISMM

    LOGICAL :: LECEARTH = .TRUE.  ! Main EC-Earth flag, always true

    LOGICAL :: ECE_CPL_AMIP = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_LIM = .FALSE.
    LOGICAL :: ECE_CPL_FESOM_FESIM = .FALSE.

    LOGICAL :: ECE_CPL_NEMO_WEIGHTED_ICE = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX = .FALSE.
    LOGICAL :: ECE_CPL_LPJG = .FALSE.
    LOGICAL :: ECE_CPL_ISMM = .FALSE.

    NAMELIST /NAMECECFG/ ECE_CPL_AMIP
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_LIM
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_WEIGHTED_ICE
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX
    NAMELIST /NAMECECFG/ ECE_CPL_LPJG
    NAMELIST /NAMECECFG/ ECE_CPL_FESOM_FESIM
    NAMELIST /NAMECECFG/ ECE_CPL_ISMM

    LOGICAL :: ECE_CLIMR = .FALSE.

    PUBLIC SURFECE_CONFIG

CONTAINS

! =============================================================================
! *** ECE_CONFIG
! =============================================================================
SUBROUTINE SURFECE_CONFIG()

    USE YOMLUN_IFSAUX, ONLY: NULOUT, NULNAM

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

    ! Configure reading of CLIMR (ICMCL) files
    ECE_CLIMR = .TRUE.

END SUBROUTINE SURFECE_CONFIG

END MODULE SURFECE
