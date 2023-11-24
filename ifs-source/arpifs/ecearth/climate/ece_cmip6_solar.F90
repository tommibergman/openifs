SUBROUTINE ECE_CMIP6_SOLAR(YEAR, MONTH, YDMODEL)

!**** *ECE_CMIP6_SOLAR*

!     PURPOSE.
!     --------

!     Updates the climatological solar forcing data for EC-Earth

!**   INTERFACE.
!     ----------

!     CALL ECE_CMIP6_SOLAR    from *UPDRGAS*

!        Explicit arguments:
!        -------------------
!        YEAR:    Year of current call
!        MONTH:   Month of current call
!        YDMODEL: Datastructure containing solar incliniation angle RSOLINC

!     METHOD.
!     -------
!     The function reads solar forcing data from a file on the first call, 
!     broadcasts it to other processes if applicable, and then provides the 
!     solar forcing value based on the specified year and month. 
!     The solar forcing data is obtained either from a fixed value for the 
!     period before 1850 or from monthly values thereafter. 

!     EXTERNALS.
!     ----------

!     ECE_CMIP6_SOLAR

!     AUTHORS.
!     --------
!     K. Wyser 2023-12

!     MODIFICATIONS.
!     --------------
!     J. Streffing 2024-01 Adapt to OpenIFS 48r1
!     ------------------------------------------------------------------

  USE PARKIND1,   ONLY: JPIM, JPRB
  USE YOMLUN,     ONLY: RESERVE_LUN, FREE_LUN
  USE YOMMP0,     ONLY: MYPROC, NPROC
  USE MPL_MODULE, ONLY: MPL_BROADCAST
  USE TYPE_MODEL, ONLY: MODEL
  USE ECE_CMIP6,  ONLY: CMIP6DATADIR, NCMIPFIXYR

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN) :: YEAR, MONTH
  TYPE(MODEL)       , INTENT(INOUT) :: YDMODEL

  LOGICAL, SAVE :: LLFIRSTCALL = .TRUE.

  REAL(KIND=JPRB), SAVE :: ZSOLMONTH(5400)

  INTEGER(KIND=JPIM), PARAMETER  :: MPITAG = 12345
  INTEGER(KIND=JPIM) :: I, IUNIT, IYR
  REAL(KIND=JPRB)    :: SKIP
  CHARACTER(LEN=132) :: SKIP_LINE

  ASSOCIATE (RSOLINC => YDMODEL%YRML_PHY_RAD%YRERDI%RSOLINC)

    IF (LLFIRSTCALL) THEN
      IF (MYPROC == 1) THEN
        IUNIT = RESERVE_LUN()
        OPEN (UNIT=IUNIT, FILE=TRIM(CMIP6DATADIR)//'/solarforcing_ref_monthly.txt', &
          &  STATUS='OLD', FORM='FORMATTED')
        !Skip the first three lines
        READ (UNIT=IUNIT, FMT=*) SKIP_LINE
        READ (UNIT=IUNIT, FMT=*) SKIP_LINE
        READ (UNIT=IUNIT, FMT=*) SKIP_LINE
        DO I = 1, SIZE(ZSOLMONTH)
          READ (UNIT=IUNIT, FMT=*) SKIP, SKIP, ZSOLMONTH(I)
        END DO
        CLOSE (UNIT=IUNIT)
        CALL FREE_LUN(IUNIT)
      END IF

      IF (NPROC > 1) THEN
        CALL MPL_BROADCAST(ZSOLMONTH, MPITAG, 1, CDSTRING='ECE_CMIP6_SOLAR:')
      END IF

      LLFIRSTCALL = .FALSE.
    END IF

    IF (YEAR < 1850 .OR. NCMIPFIXYR == 1850) THEN ! Use the mean of the period 1850-01-01 to 1873-01-28
      RSOLINC = 1360.7470703125 ! Based on solarforcing_picontrol_fx_3.2.nc from http://solarisheppa.geomar.de/cmip6
    ELSE ! Use monthly values
      IF (NCMIPFIXYR > 0) THEN
        IYR = NCMIPFIXYR
      ELSE IF (YEAR > 2299) THEN
        IYR = 2299 - 10 + Modulo(YEAR - 2300, 11) ! Repeat last solar cycle data after 2299
      ELSE
        IYR = YEAR
      END IF
      I = (IYR - 1850)*12 + MONTH
      RSOLINC = ZSOLMONTH(I)
    END IF

  END ASSOCIATE

END SUBROUTINE

