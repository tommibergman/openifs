MODULE CPLNG2_INIT_MOD

    IMPLICIT NONE

    PRIVATE

    PUBLIC CPLNG2_INIT

CONTAINS

    SUBROUTINE CPLNG2_INIT
        USE PARKIND1, ONLY: JPIM
        USE MPL_MODULE, ONLY: LMPLUSERCOMM, MPLUSERCOMM
        USE MOD_OASIS
        USE CPLNG2_DATA_MOD
        ! parameters
        CHARACTER(LEN=*), PARAMETER :: CMODEL_NAME = "OpenIFS"
        ! locals
        INTEGER(KIND=JPIM) :: comp_id
        INTEGER(KIND=JPIM) :: error
        CHARACTER(LEN=3) :: error_str

        ! Initialise OASIS coupling
        CALL OASIS_INIT_COMP(comp_id, CMODEL_NAME, error)
        IF (error /= OASIS_OK) THEN
            WRITE (error_str, '(I3)') error
            CALL ABOR1("CPLNG2_INIT: Error in OASIS_INIT_COMM: "//error_str)
        END IF

        ! Get internal communicator from OASIS
        CALL OASIS_GET_LOCALCOMM(MPLUSERCOMM, error)
        IF (error /= OASIS_OK) THEN
            WRITE (error_str, '(I3)') error
            CALL ABOR1("CPLNG2_INIT: Error in OASIS_GET_LOCALCOMM: "//error_str)
        END IF

        ! Let IFS know we have set an internal communicator
        LMPLUSERCOMM = .TRUE.
    END SUBROUTINE CPLNG2_INIT

END MODULE CPLNG2_INIT_MOD
