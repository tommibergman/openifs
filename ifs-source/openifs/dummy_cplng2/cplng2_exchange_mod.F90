MODULE CPLNG2_EXCHANGE_MOD

    IMPLICIT NONE

    PRIVATE

    PUBLIC CPLNG2_EXCHANGE

CONTAINS

SUBROUTINE CPLNG2_EXCHANGE(TSTEP,KSTAGE,YDDYNA,YDRIP)

    USE PARKIND1, ONLY: JPIM, JPRB
    USE YOMCT2,   ONLY: NSTAR2,NSTOP2
    USE YOERAD,   ONLY: YRERAD
    USE YOMRIP,   ONLY: TRIP
    USE YOMDYNA,  ONLY: TDYNA
    USE CPLNG2_DATA_MOD
    
    IMPLICIT NONE
    
    ! Argument
    INTEGER(KIND=JPIM), INTENT(IN) :: TSTEP
    INTEGER(KIND=JPIM), INTENT(IN) :: KSTAGE
    TYPE(TRIP),         INTENT(IN) :: YDRIP
    TYPE(TDYNA),        INTENT(IN) :: YDDYNA
    ! Locals
    INTEGER(KIND=JPIM) :: II
    INTEGER(KIND=JPIM) :: ILVL,ICAT
    INTEGER(KIND=JPIM) :: KINFO
    INTEGER(KIND=JPIM) :: ITIME_IN_SECONDS
    CHARACTER(LEN=3)   :: CERRSTR

    ASSOCIATE(RSTATI=>YDRIP%RSTATI, LPERPET=>YRERAD%LPERPET, LTWOTL=>YDDYNA%LTWOTL)

    ! Early return if KSTAGE is set to zero (which means ignore this field)
    IF (KSTAGE==0) RETURN

    ! If LPERPET is true, this is a perpetual run and RSTATI (needed later)
    ! doesn't represent the time since the model started. Can't handle this.
    ! Note that LPERPET is true for Aqua planet (LAQUA)
    IF (LPERPET) THEN
        CALL ABOR1("CPLNG2_EXCHANGE: Coupling doesn't work for perpetual runs.")
    ENDIF

    IF (LTWOTL) THEN
        ! Return if time is beyond last step
        IF (RSTATI>TSTEP*NSTOP2) RETURN
        ! Compute time at start of current step
        ITIME_IN_SECONDS = NINT(RSTATI - NSTAR2*TSTEP - 0.5_JPRB*TSTEP,JPIM)
    ELSE
        CALL ABOR1("CPLNG2_EXCHANGE: Can't handle LTWOTL==.FALSE. yet.")
    ENDIF

    END ASSOCIATE

END SUBROUTINE CPLNG2_EXCHANGE

END MODULE CPLNG2_EXCHANGE_MOD
