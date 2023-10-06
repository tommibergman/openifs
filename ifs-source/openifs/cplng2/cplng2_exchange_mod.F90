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
    USE MOD_OASIS
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

    ! Compute the current time in seconds since the start of the current leg.
    ! By convention, what we want to send to OASIS is the time at the beginning
    ! of the current time step intervall.
    !
    ! NOTE: Since the OASIS interface requires an integer of
    !       SELECTED_INT_KIND(9) for it's time/date argument in OASIS_Put/Get,
    !       the length of a leg is limited to 10**9 seconds (about 31 years)!
    !
    ! To figure out what the current time of the model is, we rely on RSTATI
    ! from module YOMRIP, which is updated early in the time step by UPDTIM.
    ! However, RSTATI does not contain the time at the beginning of the time
    ! step intervall. Instead, it is at time step n:
    !     IF (LTWOTL): RSTATI == (n+1/2)*dt
    !     ELSE       : ???
    ! and, hence, it is corrected accordingly.
    !   We also need to correct for previous restart legs since RSTATI is not
    ! reset during a restart. This is done by substracting NSTAR2*TSTEP, where
    ! NSTAR2 from YOMCT2 has the first time step number of the leg.
    !
    ! Note that NSTEP, which could also be used to compute the current time, is
    ! updated much later (too late) in the time step loop in CNT4!
    !
    ! Moreover, a test is implemented to avoid calling the coupler after the end
    ! of the leg, which is defined by TSTEP*NSTOP2. IFS is doing one additional
    ! time step, probably due to the two-level time stepping scheme.
    !
    IF (LTWOTL) THEN
        ! Return if time is beyond last step
        IF (RSTATI>TSTEP*NSTOP2) RETURN
        ! Compute time at start of current step
        ITIME_IN_SECONDS = NINT(RSTATI - NSTAR2*TSTEP - 0.5_JPRB*TSTEP,JPIM)
    ELSE
        CALL ABOR1("CPLNG2_EXCHANGE: Can't handle LTWOTL==.FALSE. yet.")
    ENDIF

    DO II=1,SIZE(CPLNG2_FLD)

        ! Do nothing if this field is not exchanged in the current stage
        IF (CPLNG2_FLD(II)%STAGE/=KSTAGE) CYCLE

        ! Decide whether the coupling field is to be sent (put) or
        ! received (get)
        SELECT CASE (CPLNG2_FLD(II)%INOUT)

        CASE (OASIS_Out) ! Call OASIS_PUT for couple fields that are sent

            DO ICAT=1,CPLNG2_FLD(II)%NUM_CAT
                DO ILVL=1,CPLNG2_FLD(II)%NUM_LVL
                    !WRITE(*,*) "CPLNG2_EXCHANGE_MOD: OASIS_PUT: II, NAME ",II,CPLNG2_FLD(II)%NAME, &
                    !& CPLNG2_FLD(II)%ID(ILVL,ICAT)
                    CALL OASIS_PUT(CPLNG2_FLD(II)%ID(ILVL,ICAT),  &
                    &              ITIME_IN_SECONDS,             &
                    &              CPLNG2_FLD(II)%D(:,ILVL,ICAT), &
                    &              KINFO)

                    SELECT CASE (KINFO)

                    CASE (OASIS_Sent,      &
                    &     OASIS_LocTrans,  &
                    &     OASIS_ToRest,    &
                    &     OASIS_Output,    &
                    &     OASIS_SentOut,   &
                    &     OASIS_ToRestOut, &
                    &     OASIS_Waitgroup, &
                    &     OASIS_Ok         )

                        CONTINUE

                    CASE DEFAULT

                        WRITE (CERRSTR,'(I3)') KINFO
                        CALL ABOR1("CPLNG2_EXCHANGE: Error in OASIS_PUT: "//CERRSTR)

                    END SELECT
                ENDDO
            ENDDO

        CASE (OASIS_In) ! Call OASIS_GET for couple fields that are received

            DO ICAT=1,CPLNG2_FLD(II)%NUM_CAT
                DO ILVL=1,CPLNG2_FLD(II)%NUM_LVL
                    !WRITE(*,*) "CPLNG2_EXCHANGE_MOD: OASIS_GET NAME, time: ",II,CPLNG2_FLD(II)%NAME, &
                    !            & CPLNG2_FLD(II)%ID(ILVL,ICAT), ITIME_IN_SECONDS
                    CALL OASIS_GET(CPLNG2_FLD(II)%ID(ILVL,ICAT),  &
                    &              ITIME_IN_SECONDS,             &
                    &              CPLNG2_FLD(II)%D(:,ILVL,ICAT), &
                    &              KINFO)

                    SELECT CASE (KINFO)

                    CASE (OASIS_Recvd,       &
                    &     OASIS_FromRest,    &
                    &     OASIS_Input,       &
                    &     OASIS_RecvOut,     &
                    &     OASIS_FromRestOut, &
                    &     OASIS_Ok           )

                        CONTINUE

                    CASE DEFAULT

                        WRITE (CERRSTR,'(I3)') KINFO
                        CALL ABOR1("CPLNG2_EXCHANGE: Error in OASIS_GET: "//CERRSTR)

                    END SELECT
                ENDDO
            ENDDO

        CASE DEFAULT ! Everthing else is an error

            WRITE (CERRSTR,'(I3,A,I5,A)') II,' (INOUT=',CPLNG2_FLD(II)%INOUT,')'
            CALL ABOR1("CPLNG2_EXCHANGE: Error in definition of field no. "//CERRSTR)

        END SELECT

    ENDDO
    
    END ASSOCIATE

END SUBROUTINE CPLNG2_EXCHANGE

END MODULE CPLNG2_EXCHANGE_MOD
