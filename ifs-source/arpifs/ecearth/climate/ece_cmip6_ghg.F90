SUBROUTINE ECE_CMIP6_GHG(IYR, IMN)

!**** *ECE_CMIP6_GHG*

!     PURPOSE.
!     --------

!     Updates the climatological greenhouse gas forcing data for EC-Earth

!**   INTERFACE.
!     ----------

!     CALL ECE_CMIP6_GHG    from *UPDRGAS*

!        Explicit arguments:
!        -------------------
!        IYR:     Year of current call
!        IMN:     Month of current call
!        YDMODEL: Datastructure containing solar incliniation angle RSOLINC

!     METHOD.
!     -------
!     The function reads greehouse gas forcing data from a file on the first call,
!     broadcasts it to other processes if applicable, and then provides the
!     greenhouse gas forcing values based on the specified year and month.

!     EXTERNALS.
!     ----------

!     ECE_CMIP6_GHG

!     AUTHORS.
!     --------
!     K. Wyser 2023-12

!     MODIFICATIONS.
!     --------------
!     J. Streffing 2024-01 Adapt to OpenIFS 48r1
!     ------------------------------------------------------------------

  USE PARKIND1,      ONLY: JPIM, JPRB
  USE YOMLUN, ONLY: NULOUT
  USE MPL_MODULE,    ONLY: MPL_BROADCAST
  USE YOERDI, ONLY: YRERDI
  USE ECE_CMIP6,     ONLY: CMIP6DATADIR, NCMIPFIXYR, NCMIPFIXYR_CH4, &
                       & SSPNAME, LA4xCO2, L1PCTCO2, LGHGMONTHLY

  USE NETCDF

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN) :: IYR, IMN
  INTEGER(KIND=JPIM) :: IYR1, IYR2, IMN0
  INTEGER(KIND=JPIM), SAVE :: IYR2OLD, IMNOLD

  REAL(KIND=JPRB) :: ZCO2RMWG, ZCH4RMWG, ZN2ORMWG, ZNO2RMWG, ZC11RMWG, ZC12RMWG

  INTEGER(KIND=JPIM) :: JGAS
  REAL(KIND=JPRB) :: ZFIXNO2, ZCONC(6), ZGRADC
  REAL(KIND=JPRB), SAVE :: ZFCONC(2, 5) = 0._JPRB

! Increase rate of co2 per year in 1pctCO2 experiment
  REAL(KIND=JPRB), PARAMETER :: RCO2INC = 0.01_JPRB
! co2 level of co2 the Abrupt4xCO2 experiment
  REAL(KIND=JPRB), PARAMETER :: R4xCO2 = 4.0_JPRB

  LOGICAL, SAVE :: FIRST_CALL = .TRUE.

  ASSOCIATE (RCARDI => YRERDI%RCARDI, RCFC11 => YRERDI%RCFC11, &
   & RCFC12 => YRERDI%RCFC12, RCH4 => YRERDI%RCH4, RN2O => YRERDI%RN2O, &
   & RNO2 => YRERDI%RNO2)

! set constants

!ZGASRMWG = ZGASMWG / ZAIRMWG
    ZCO2RMWG = 1.5191923_JPRB
    ZCH4RMWG = 0.5537798_JPRB
    ZN2ORMWG = 1.5192613_JPRB
    ZNO2RMWG = 1.5880566_JPRB
    ZC11RMWG = 4.7417535_JPRB
    ZC12RMWG = 4.1737660_JPRB

    ZFIXNO2 = 500._JPRB

! find year before or after
    IF (IMN <= 6) THEN
      IYR1 = IYR - 1
      IYR2 = IYR
      IMN0 = IMN + 6
    ELSEIF (IMN > 6) THEN
      IYR1 = IYR
      IYR2 = IYR + 1
      IMN0 = IMN - 6
    END IF

    IF (IMN .NE. IMNOLD) THEN
      IF (LGHGMONTHLY) THEN
        CALL READCMIP6GHGDATA(IYR, IMN)
      ELSE IF (IYR2 .NE. IYR2OLD) THEN
        IF (FIRST_CALL) THEN
          CALL READCMIP6GHGDATA(IYR1)
        END IF
        CALL READCMIP6GHGDATA(IYR2)
        IYR2OLD = IYR2
      END IF
      FIRST_CALL = .FALSE.

      IF (LGHGMONTHLY) THEN
        ZCONC(1:5) = ZFCONC(2, 1:5)
      ELSE
        DO JGAS = 1, 5
          ZGRADC = (ZFCONC(2, JGAS) - ZFCONC(1, JGAS))/12._JPRB
          ZCONC(JGAS) = ZFCONC(1, JGAS) + IMN0*ZGRADC
        END DO
      END IF
      ZCONC(6) = ZFIXNO2

      WRITE (NULOUT, *) 'ECE_CMIP6_GHG:'
      IF (NCMIPFIXYR <= 0) THEN
        WRITE (NULOUT, FMT='('' IYR ='',I4,'' IMN ='',I4,'' IMN0 ='',I4 &
          & ,'' IYR1='',I4,'' IYR2='',I4)') &
          & IYR, IMN, IMN0, IYR1, IYR2
      ELSE
        WRITE (NULOUT, FMT='('' NCMIPFIXYR ='',I4,'' IMN ='',I4,'' IMN0 ='',I4 &
          & ,'' LA4xCO2='',L4,'' L1pctCO2='',L4)') &
          & NCMIPFIXYR, IMN, IMN0, LA4xCO2, L1pctCO2
      END IF
      DO JGAS = 1, 6
        WRITE (NULOUT, *) 'JGAS=', JGAS, ' ZCONC(JGAS)=', ZCONC(JGAS)
      END DO

      RCARDI = ZCONC(1)*1.E-06_JPRB*ZCO2RMWG
      RCH4 = ZCONC(2)*1.E-09_JPRB*ZCH4RMWG
      RN2O = ZCONC(3)*1.E-09_JPRB*ZN2ORMWG
      RCFC11 = ZCONC(4)*1.E-12_JPRB*ZC11RMWG
      RCFC12 = ZCONC(5)*1.E-12_JPRB*ZC12RMWG
      RNO2 = ZCONC(6)*1.E-13_JPRB*ZNO2RMWG

      IMNOLD = IMN
    END IF

  END ASSOCIATE

  RETURN

CONTAINS

  SUBROUTINE READCMIP6GHGDATA(IYEAR, IMONTH)

    USE YOMMP0, ONLY: MYPROC, NPROC
    USE MPL_MODULE, ONLY: MPL_BROADCAST

    INTEGER, PARAMETER :: RPRC = 1, ITAG = 12345

    INTEGER(KIND=JPIM), INTENT(IN)  :: IYEAR
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: IMONTH

    REAL(KIND=JPRB), ALLOCATABLE :: ZFCONC_BUF(:)

    CHARACTER(LEN=255) :: FILENAME, FILEID, FILEID_CH4
    CHARACTER(LEN=10)  :: NC_NAME
    INTEGER(KIND=JPIM) :: I, IFIRSTYR, JYEAR
    INTEGER            :: IUNIT, IDX, IDX_CH4, NY_VARID

    CHARACTER(LEN=13)  :: TIMEPERIOD
    INTEGER            :: NMON, IMON

    CHARACTER(LEN=*), PARAMETER :: NY_NAME = 'time'
    REAL(KIND=JPRB) :: ZZ_YEARS

    CHARACTER(LEN=*), PARAMETER :: CO2_FILE = 'mole-fraction-of-carbon-dioxide-in-air'
    CHARACTER(LEN=*), PARAMETER :: CO2_NAME = 'mole_fraction_of_carbon_dioxide_in_air'
    INTEGER(KIND=JPIM) :: CO2_VARID
    REAL(KIND=JPRB) :: ZZCO2

    CHARACTER(LEN=*), PARAMETER :: CH4_FILE = 'mole-fraction-of-methane-in-air'
    CHARACTER(LEN=*), PARAMETER :: CH4_NAME = 'mole_fraction_of_methane_in_air'
    INTEGER(KIND=JPIM) :: CH4_VARID
    REAL(KIND=JPRB) :: ZZCH4

    CHARACTER(LEN=*), PARAMETER :: N2O_FILE = 'mole-fraction-of-nitrous-oxide-in-air'
    CHARACTER(LEN=*), PARAMETER :: N2O_NAME = 'mole_fraction_of_nitrous_oxide_in_air'
    INTEGER(KIND=JPIM) :: N2O_VARID
    REAL(KIND=JPRB) :: ZZN2O

    CHARACTER(LEN=*), PARAMETER :: CFC11_FILE = 'mole-fraction-of-cfc11eq-in-air'
    CHARACTER(LEN=*), PARAMETER :: CFC11_NAME = 'mole_fraction_of_cfc11eq_in_air'
    INTEGER(KIND=JPIM) :: CFC11_VARID
    REAL(KIND=JPRB) :: ZZCFC11

    CHARACTER(LEN=*), PARAMETER :: CFC12_FILE = 'mole-fraction-of-cfc12-in-air'
    CHARACTER(LEN=*), PARAMETER :: CFC12_NAME = 'mole_fraction_of_cfc12_in_air'
    INTEGER(KIND=JPIM) :: CFC12_VARID
    REAL(KIND=JPRB) :: ZZCFC12

    IF (LGHGMONTHLY) THEN
      WRITE (NULOUT, *) 'ECE_CMIP6_GHG: read monthly CMIP6 GHG concentrations for ', IYEAR, IMONTH
    ELSE
      WRITE (NULOUT, *) 'ECE_CMIP6_GHG: read yearly CMIP6 GHG concentrations for ', IYEAR
    ENDIF

    ALLOCATE (ZFCONC_BUF(SIZE(ZFCONC)))

    IF (MYPROC == RPRC) THEN

      ! All GHG except CH4
      IF (NCMIPFIXYR .GT. 0) THEN
        JYEAR = NCMIPFIXYR
      ELSE
        JYEAR = IYEAR
      END IF

      IF (JYEAR < 2015) THEN
        ! yearly or monthly
        IF (LGHGMONTHLY) THEN
          TIMEPERIOD = '000001-201412'
        ELSE
          TIMEPERIOD = '0000-2014'
        END IF

        FILEID = '_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        ! The first year in the data set is 0
        IFIRSTYR = 0
      ELSE
        ! yearly or monthly
        IF (LGHGMONTHLY) THEN
          TIMEPERIOD = '201501-250012'
        ELSE
          TIMEPERIOD = '2015-2500'
        END IF

        SELECT CASE (TRIM(SSPNAME))
        CASE ("SSP1-1.9")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp119-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP1-2.6")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP2-4.5")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP3-7.0")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP3-LowNTCF")
          FILEID = '_input4MIPs_GHGConcentrations_AerChemMIP_UoM-AIM-ssp370-lowNTCF-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP4-3.4")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-GCAM4-ssp434-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP4-6.0")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-GCAM4-ssp460-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP5-3.4-OS")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp534-over-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP5-8.5")
          FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE DEFAULT
          ! default to SSP3-7.0 but only to finish historical runs (last step of historical)
          IF (.NOT. FIRST_CALL) THEN
            FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
          ELSE
            CALL ABOR1('ECE_CMIP6_GHG : unknown scenario')
          END IF
        END SELECT

        ! The first year in the data set is 2015
        IFIRSTYR = 2015
      END IF

      IF (LGHGMONTHLY) THEN
        NMON = 12
        IMON = IMONTH
      ELSE
        NMON = 1
        IMON = 1
      END IF
      IF (NCMIPFIXYR .GT. 0) THEN
        IDX = (NCMIPFIXYR - IFIRSTYR)*NMON + IMON
      ELSE
        IDX = (MIN(IYEAR, 2500) - IFIRSTYR)*NMON + IMON
      END IF

      ! Specific case of CH4
      IF (NCMIPFIXYR_CH4 .GT. 0) THEN
        JYEAR = NCMIPFIXYR_CH4
      ELSE
        JYEAR = IYEAR
      END IF
      IF (JYEAR < 2015) THEN
        FILEID_CH4 = '_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        ! The first year in the data set is 0
        IFIRSTYR = 0
      ELSE
        SELECT CASE (TRIM(SSPNAME))
        CASE ("SSP1-1.9")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp119-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP1-2.6")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP2-4.5")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP3-7.0")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP3-LowNTCF")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_AerChemMIP_UoM-AIM-ssp370-lowNTCF-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP4-3.4")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-GCAM4-ssp434-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP4-6.0")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-GCAM4-ssp460-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP5-3.4-OS")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp534-over-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP5-8.5")
          FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE DEFAULT
          ! default to SSP3-7.0 but only to finish historical runs (last step of historical)
          IF (.NOT. FIRST_CALL) THEN
            FILEID_CH4 = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
          ELSE
            CALL ABOR1('ECE_CMIP6_GHG : unknown scenario')
          END IF
        END SELECT
        ! The first year in the data set is 2015
        IFIRSTYR = 2015
      END IF

      IF (NCMIPFIXYR_CH4 .GT. 0) THEN
        IDX_CH4 = (NCMIPFIXYR_CH4 - IFIRSTYR)*NMON + IMON
      ELSE
        IDX_CH4 = (MIN(IYEAR, 2500) - IFIRSTYR)*NMON + IMON
      END IF

      ! Read in CO2 data
      WRITE (FILENAME, '(A)') TRIM(CMIP6DATADIR)//'/'//TRIM(CO2_FILE)//TRIM(FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))

      ! Check timestamp only from CO2 file
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, NY_NAME, NY_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, NY_VARID, ZZ_YEARS, start=(/IDX/)))
      WRITE (NULOUT, *) 'ZZ_YEARS=', ZZ_YEARS

      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CO2_NAME, CO2_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CO2_VARID, ZZCO2, start=(/1, IDX/)))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCO2=', ZZCO2

      ! Read in CH4 data
      WRITE (FILENAME, '(A)') TRIM(CMIP6DATADIR)//'/'//TRIM(CH4_FILE)//TRIM(FILEID_CH4)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CH4_NAME, CH4_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CH4_VARID, ZZCH4, start=(/1, IDX_CH4/)))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCH4=', ZZCH4

      ! Read in N2O data
      WRITE (FILENAME, '(A)') TRIM(CMIP6DATADIR)//'/'//TRIM(N2O_FILE)//TRIM(FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, N2O_NAME, N2O_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, N2O_VARID, ZZN2O, start=(/1, IDX/)))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZN2O=', ZZN2O

      ! Read in CFC11 data
      WRITE (FILENAME, '(A)') TRIM(CMIP6DATADIR)//'/'//TRIM(CFC11_FILE)//TRIM(FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CFC11_NAME, CFC11_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CFC11_VARID, ZZCFC11, start=(/1, IDX/)))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCFC11=', ZZCFC11

      ! Read in CFC12 data
      WRITE (FILENAME, '(A)') TRIM(CMIP6DATADIR)//'/'//TRIM(CFC12_FILE)//TRIM(FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CFC12_NAME, CFC12_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CFC12_VARID, ZZCFC12, start=(/1, IDX/)))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCFC12=', ZZCFC12

      ! replace IYR1 with IYR2 and fill IYR2 with new values
      ZFCONC(1, :) = ZFCONC(2, :)

      ZFCONC(2, 2) = ZZCH4
      ZFCONC(2, 3) = ZZN2O
      ZFCONC(2, 4) = ZZCFC11
      ZFCONC(2, 5) = ZZCFC12

      IF (LA4xCO2) THEN
        ! Abrupt4xCO2 - The CO2 abruptly quadrupled of the NCMIPFIXYR level and held constant
        ZFCONC(2, 1) = ZZCO2*R4xCO2
        WRITE (NULOUT, '(''LA4xCO2='',L1,'' CO2 concentration prescribed as '',F8.2)') &
          &   LA4xCO2, ZFCONC(2, 1)

      ELSE IF (L1PCTCO2) THEN
        ! 1pctCO2 - the CO2 increase at the rate of 1%/year. In contrast to the CMIP5 protocol,
        ! in CMIP6 the increase is continued to the end of the simulation,
        ! even after CO2 has quadrupled.

        I = IYEAR - NCMIPFIXYR
        ZFCONC(2, 1) = ZZCO2*EXP(I*LOG(1._JPRB + RCO2INC))

        WRITE (NULOUT, '(''L1PCTCO2='',L1,'' CO2 concentration prescribed as '',F8.2)') &
          &   L1PCTCO2, ZFCONC(2, 1)

      ELSE
        ! default, use the CO2 value read from the file
        ZFCONC(2, 1) = ZZCO2
      END IF

      ZFCONC_BUF = RESHAPE(ZFCONC, SHAPE(ZFCONC_BUF))
    END IF
    ! distribute to all PROCS
    CALL MPL_BROADCAST(ZFCONC_BUF, KTAG=ITAG + 3, KROOT=RPRC, CDSTRING='CMIP6GHG-1: ')
    ZFCONC = RESHAPE(ZFCONC_BUF, SHAPE(ZFCONC))
    DEALLOCATE (ZFCONC_BUF)

    RETURN
  END SUBROUTINE READCMIP6GHGDATA

  SUBROUTINE CHECKGHG(STATUS)
    INTEGER, INTENT(IN) :: STATUS

    IF (STATUS /= NF90_NOERR) THEN
      CALL ABOR1('ECE_CMIP6_GHG: '//TRIM(NF90_STRERROR(STATUS)))
    END IF
    RETURN
  END SUBROUTINE CHECKGHG

END SUBROUTINE ECE_CMIP6_GHG
