MODULE SURFECE
    ! Some logic switches for EC-Earth need to be available in surf.
    ! As the main module ECEARTH is part of the later compiled arpifs,
    ! which depends on surf, we create this new module.

    ! Authors:
    ! 2023-11-22: Jan Streffing (AWI): Created separate module from ECEARTH
    ! 2026-02-28: Jorge Bernales (DMI): Add storage/API for LANDICE features

    USE PARKIND1,ONLY : JPIM, JPRB

    IMPLICIT NONE

    PRIVATE

    PUBLIC LECEARTH

    PUBLIC ECE_CPL_AMIP
    PUBLIC ECE_CPL_NEMO_LIM
    PUBLIC ECE_CPL_NEMO_PISCES
    PUBLIC ECE_CPL_FESOM_FESIM
    PUBLIC ECE_CPL_FESOM_RECOM
    PUBLIC ECE_CPL_LPJG
    PUBLIC ECE_CPL_LPJG_CO2
    PUBLIC ECE_CPL_ISMM

    PUBLIC ECE_CLIMR

    PUBLIC ECE_CPL_NEMO_WEIGHTED_ICE
    PUBLIC ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX

    PUBLIC NISMRGNS
    PUBLIC ISMRGNPFX
    PUBLIC ECE_LANDICE
    PUBLIC ECE_LANDICE_THRESH
    PUBLIC ECE_LANDICE_READY
    PUBLIC ECE_LANDICE_ACTIVE_BLOCK

    ! Ice sheet albedo evolution parameters (shared by srfsn_asn and srfsn_lwimp)
    PUBLIC ECE_LANDICE_ALB_MIN, ECE_LANDICE_ALB_REFROZ
    PUBLIC ECE_LANDICE_ALB_FIRN, ECE_LANDICE_ALB_FRESH
    PUBLIC ECE_LANDICE_TAU_DRY, ECE_LANDICE_TAU_WET

    LOGICAL :: LECEARTH = .TRUE.  ! Main EC-Earth flag, always true

    LOGICAL :: ECE_CPL_AMIP = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_LIM = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_PISCES = .FALSE.
    LOGICAL :: ECE_CPL_FESOM_FESIM = .FALSE.
    LOGICAL :: ECE_CPL_FESOM_RECOM = .FALSE.
    LOGICAL :: ECE_CPL_LPJG = .FALSE.
    LOGICAL :: ECE_CPL_LPJG_CO2 = .FALSE.
    LOGICAL :: ECE_CPL_ISMM = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_WEIGHTED_ICE = .FALSE.
    LOGICAL :: ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX = .FALSE.

    INTEGER(KIND=JPIM) :: NISMRGNS
    CHARACTER(len=10), DIMENSION(:), ALLOCATABLE :: ISMRGNPFX
    LOGICAL :: ECE_LANDICE = .FALSE.
    REAL(KIND=JPRB) :: ECE_LANDICE_THRESH = 0.5_JPRB  ! Threshold for ice sheet mask [0-1]; cells above treated as glacier
    LOGICAL :: ECE_LANDICE_READY = .FALSE.
    INTEGER(KIND=JPIM) :: ECE_LANDICE_ACTIVE_BLOCK = -1_JPIM

    ! Ice sheet albedo evolution parameters
    REAL(KIND=JPRB), PARAMETER :: ECE_LANDICE_ALB_MIN    = 0.60_JPRB  ! Minimum albedo over ice sheets
    REAL(KIND=JPRB), PARAMETER :: ECE_LANDICE_ALB_REFROZ = 0.65_JPRB  ! Albedo for refreezing snow
    REAL(KIND=JPRB), PARAMETER :: ECE_LANDICE_ALB_FIRN   = 0.75_JPRB  ! Albedo for firn
    REAL(KIND=JPRB), PARAMETER :: ECE_LANDICE_ALB_FRESH  = 0.85_JPRB  ! Albedo for fresh snow
    REAL(KIND=JPRB), PARAMETER :: ECE_LANDICE_TAU_DRY    = 1._JPRB/30._JPRB  ! Dry decay rate [1/day]
    REAL(KIND=JPRB), PARAMETER :: ECE_LANDICE_TAU_WET    = 1._JPRB           ! Wet decay rate [1/day]
!$OMP THREADPRIVATE(ECE_LANDICE_ACTIVE_BLOCK)
    REAL(KIND=JPRB), ALLOCATABLE :: ECE_LANDICE_MASK_B(:,:) ! In blocked form (NPROMA,NGPBLKS)

    NAMELIST /NAMECECFG/ ECE_CPL_AMIP
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_LIM
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_PISCES
    NAMELIST /NAMECECFG/ ECE_CPL_FESOM_FESIM
    NAMELIST /NAMECECFG/ ECE_CPL_FESOM_RECOM
    NAMELIST /NAMECECFG/ ECE_CPL_LPJG
    NAMELIST /NAMECECFG/ ECE_CPL_LPJG_CO2
    NAMELIST /NAMECECFG/ ECE_CPL_ISMM
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_WEIGHTED_ICE
    NAMELIST /NAMECECFG/ ECE_CPL_NEMO_CONSERVATIVE_HEATFLUX
    NAMELIST /NAMECECFG/ NISMRGNS, ISMRGNPFX
    NAMELIST /NAMECECFG/ ECE_LANDICE
    NAMELIST /NAMECECFG/ ECE_LANDICE_THRESH

    LOGICAL :: ECE_CLIMR = .FALSE.

    PUBLIC SURFECE_CONFIG
    PUBLIC SURFECE_SET_LANDICE_LOCAL
    PUBLIC SURFECE_SET_ACTIVE_BLOCK
    PUBLIC SURFECE_GET_LANDICE
    PUBLIC SURFECE_GET_LANDICE_BLOCK
    PUBLIC SURFECE_RESET_LANDICE

CONTAINS

! =============================================================================
! *** ECE_CONFIG
! =============================================================================
SUBROUTINE SURFECE_CONFIG()

    USE YOMLUN_IFSAUX, ONLY: NULOUT, NULNAM

    !Allocate ISM regions
    ALLOCATE(ISMRGNPFX(10))

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

! =============================================================================
! *** SURFECE_SET_LANDICE_LOCAL
! =============================================================================
SUBROUTINE SURFECE_SET_LANDICE_LOCAL(PLANDICE_MASK_B)
    REAL(KIND=JPRB), INTENT(IN) :: PLANDICE_MASK_B(:,:)

    IF (.NOT. ECE_LANDICE) THEN
        CALL SURFECE_RESET_LANDICE()
        RETURN
    ENDIF

    ! Avoid stale wrong-sized arrays across restarts/re-init
    IF (ALLOCATED(ECE_LANDICE_MASK_B)) THEN
        IF (SIZE(ECE_LANDICE_MASK_B,1) /= SIZE(PLANDICE_MASK_B,1) .OR. &
            SIZE(ECE_LANDICE_MASK_B,2) /= SIZE(PLANDICE_MASK_B,2)) THEN
            DEALLOCATE(ECE_LANDICE_MASK_B)
        ENDIF
    ENDIF

    ! Allocate exact incoming shape
    IF (.NOT. ALLOCATED(ECE_LANDICE_MASK_B)) THEN
        ALLOCATE(ECE_LANDICE_MASK_B(SIZE(PLANDICE_MASK_B,1), SIZE(PLANDICE_MASK_B,2)))
    ENDIF

    ! Copy in data (to decouple caller lifetime) and publish that data is ready
    ECE_LANDICE_MASK_B(:,:) = PLANDICE_MASK_B(:,:)
    ECE_LANDICE_READY = .TRUE.
    ECE_LANDICE_ACTIVE_BLOCK = -1_JPIM
END SUBROUTINE SURFECE_SET_LANDICE_LOCAL

! =============================================================================
! *** SURFECE_SET_ACTIVE_BLOCK
! =============================================================================
SUBROUTINE SURFECE_SET_ACTIVE_BLOCK(KBL)
    INTEGER(KIND=JPIM), INTENT(IN) :: KBL
    ECE_LANDICE_ACTIVE_BLOCK = KBL
END SUBROUTINE SURFECE_SET_ACTIVE_BLOCK

! =============================================================================
! *** SURFECE_GET_LANDICE
! =============================================================================
SUBROUTINE SURFECE_GET_LANDICE(PLANDICE)
    REAL(KIND=JPRB), INTENT(OUT) :: PLANDICE(:)

    CALL SURFECE_GET_LANDICE_BLOCK(ECE_LANDICE_ACTIVE_BLOCK, PLANDICE)
END SUBROUTINE SURFECE_GET_LANDICE

! =============================================================================
! *** SURFECE_GET_LANDICE_BLOCK
! =============================================================================
SUBROUTINE SURFECE_GET_LANDICE_BLOCK(KBL, PLANDICE)
    INTEGER(KIND=JPIM), INTENT(IN) :: KBL
    REAL(KIND=JPRB),    INTENT(OUT) :: PLANDICE(:)

    ! Default so we ensure deterministic output
    PLANDICE(:) = 0.0_JPRB

    ! Safety
    IF (.NOT. ECE_LANDICE) RETURN
    IF (.NOT. ECE_LANDICE_READY) RETURN
    IF (.NOT. ALLOCATED(ECE_LANDICE_MASK_B)) RETURN
    IF (KBL < 1 .OR. KBL > SIZE(ECE_LANDICE_MASK_B,2)) RETURN

    ! Fill from stored block; limit to avoid bounds issues
    PLANDICE(1:MIN(SIZE(PLANDICE),SIZE(ECE_LANDICE_MASK_B,1))) = &
        ECE_LANDICE_MASK_B(1:MIN(SIZE(PLANDICE),SIZE(ECE_LANDICE_MASK_B,1)),KBL)
END SUBROUTINE SURFECE_GET_LANDICE_BLOCK

! =============================================================================
! *** SURFECE_RESET_LANDICE
! =============================================================================
SUBROUTINE SURFECE_RESET_LANDICE()
    ! Cleanup/reinit hook for job phases/restarts/tests
    IF (ALLOCATED(ECE_LANDICE_MASK_B)) DEALLOCATE(ECE_LANDICE_MASK_B)
    ECE_LANDICE_READY = .FALSE.
    ECE_LANDICE_ACTIVE_BLOCK = -1_JPIM
END SUBROUTINE SURFECE_RESET_LANDICE

END MODULE SURFECE
