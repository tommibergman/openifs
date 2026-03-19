SUBROUTINE ECE_CMIP_SOLAR(YEAR, MONTH, YDERDI)

!**** *ECE_CMIP_SOLAR*

!     PURPOSE.
!     --------

!     Updates the climatological solar forcing data for EC-Earth

!**   INTERFACE.
!     ----------

!     CALL ECE_CMIP_SOLAR    from *UPDRGAS*

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

!     ECE_CMIP_SOLAR

!     AUTHORS.
!     --------
!     K. Wyser 2023-12

!     MODIFICATIONS.
!     --------------
!     J. Streffing 2024-01 Adapt to OpenIFS 48r1
!     ------------------------------------------------------------------

  USE PARKIND1,   ONLY: JPIM, JPRB
  USE YOMLUN,     ONLY: RESERVE_LUN, FREE_LUN, NULOUT
  USE YOMMP0,     ONLY: MYPROC, NPROC
  USE MPL_MODULE, ONLY: MPL_BROADCAST
  USE YOERDI,     ONLY: TERDI
  USE ECE_CMIP,   ONLY: CMIP6DATADIR, CMIP7DATADIR, NCMIPFIXYR
  USE ECE_CMIP,   ONLY: LCMIP6, LCMIP7
  USE NETCDF

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN) :: YEAR, MONTH
  TYPE(TERDI),     INTENT(INOUT) :: YDERDI

  LOGICAL, SAVE :: LLFIRSTCALL = .TRUE.

  REAL(KIND=JPRB), SAVE :: ZSOLMONTH(5400)

  INTEGER(KIND=JPIM), PARAMETER  :: MPITAG = 12345
  INTEGER(KIND=JPIM) :: I, IUNIT, IYR, IYREF, ISIZE
  REAL(KIND=JPRB)    :: SKIP
  CHARACTER(LEN=132) :: SKIP_LINE

  REAL(KIND=JPRB), ALLOCATABLE, SAVE :: CMIP7_TSI(:)
  CHARACTER(LEN=255) :: FILENAME
  INTEGER            :: TSI_VARID

  ASSOCIATE (RSOLINC => YDERDI%RSOLINC)

    IF (LCMIP6) THEN
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
          CALL MPL_BROADCAST(ZSOLMONTH, MPITAG, 1, CDSTRING='ECE_CMIP_SOLAR:')
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

    ELSEIF (LCMIP7) THEN

      IYR = YEAR
      IF (NCMIPFIXYR > 0) IYR = NCMIPFIXYR
      IF (IYR > 2299) IYR = 2299 - 10 + MODULO(IYR - 2300, 11) ! Repeat last solar cycle data after 2299

      ! CMIP7 forcing dataset
      ! - Note historical and scenario files overlap over 2022-2023
      ! 
      FILENAME = TRIM(CMIP7DATADIR)//'/solar/multiple_input4MIPs_solar_'
      IF (NCMIPFIXYR == 1850 .OR. IYR < 1850) THEN
        FILENAME = TRIM(FILENAME)//'CMIP_SOLARIS-HEPPA-CMIP-4-6_gn.nc'
        ISIZE = 1
        IYREF = 0
      ELSE IF (IYR >= 1850 .AND. IYR <= 2022)
        FILENAME = TRIM(FILENAME)//'CMIP_SOLARIS-HEPPA-CMIP-4-6_gn_185001-202312.nc'
        ISIZE = 2088
        IYREF = 1850
      ELSE IF (IYR >= 2023)
        FILENAME = TRIM(FILENAME)//'ScenarioMIP_SOLARIS-HEPPA-ScenarioMIP-4-6_gn_202201-229912.nc'
        ISIZE = 3336
        IYREF = 2022
      END IF

      IF (LLFIRSTCALL) THEN
        ALLOCATE (CMIP7_TSI(ISIZE))
        IF (MYPROC == 1) THEN
          WRITE (NULOUT, '(A)') 'ECE_CMIP_SOLAR read file:', TRIM(FILENAME)
          CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
          CALL CHECKGHG(NF90_INQ_VARID(IUNIT, 'tsi', TSI_VARID))
          CALL CHECKGHG(NF90_GET_VAR(IUNIT, TSI_VARID, CMIP7_TSI))
          CALL CHECKGHG(NF90_CLOSE(IUNIT))
        END IF
        IF (NPROC > 1) THEN
          CALL MPL_BROADCAST(CMIP7_TSI, MPITAG, 1, CDSTRING='ECE_CMIP_SOLAR:')
        END IF
        LLFIRSTCALL = .FALSE.
      END IF

      ! piControl uses annual mean, otherwise monthly
      IF (NCMIPFIXYR == 1850 .OR. IYR < 1850) THEN
        I = 1
      ELSE
        I = (IYR - IYREF)*12 + MONTH
      END IF
      RSOLINC = CMIP7_TSI(I)
    END IF

  END ASSOCIATE

  RETURN

CONTAINS

  SUBROUTINE CHECKGHG(STATUS)
    INTEGER, INTENT(IN) :: STATUS

    IF (STATUS /= NF90_NOERR) THEN
      CALL ABOR1('ECE_CMIP_SOLAR: '//TRIM(NF90_STRERROR(STATUS)))
    END IF
    RETURN
  END SUBROUTINE CHECKGHG

END SUBROUTINE ECE_CMIP_SOLAR

