SUBROUTINE ECE_CMIP_GHG(IYR, IMN, YDERDI)

!**** *ECE_CMIP_GHG*

!     PURPOSE.
!     --------

!     Updates the climatological greenhouse gas forcing data for EC-Earth

!**   INTERFACE.
!     ----------

!     CALL ECE_CMIP_GHG    from *UPDRGAS*

!        Explicit arguments:
!        -------------------
!        IYR:    Year of current call
!        IMN:    Month of current call
!        YDERDI: Type with mass mixing ratios of trace gases

!     METHOD.
!     -------
!     The function reads greehouse gas forcing data from a file on the first call,
!     broadcasts it to other processes if applicable, and then provides the
!     greenhouse gas forcing values based on the specified year and month.

!     EXTERNALS.
!     ----------

!     ECE_CMIP_GHG

!     AUTHORS.
!     --------
!     K. Wyser 2023-12

!     MODIFICATIONS.
!     --------------
!     J. Streffing/P. Le Sager 2024 - Adapt to OpenIFS 48r1; more robust
!     ------------------------------------------------------------------

  USE PARKIND1,      ONLY: JPIM, JPRB
  USE YOMLUN_IFSAUX, ONLY: NULOUT
  USE MPL_MODULE,    ONLY: MPL_BROADCAST
  USE YOERDI,        ONLY: TERDI
  USE ECE_CMIP,      ONLY: CMIP6DATADIR, CMIP7DATADIR, NCMIPFIXYR, NCMIPFIXYR_CH4, &
                       & SCENARIONAME, LANxCO2, RNxCO2, L1PCTCO2, LGHGMONTHLY, LCMIP6, LCMIP7
  USE YOEMETH,       ONLY : RQLIM
  USE NETCDF

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN) :: IYR, IMN
  TYPE(TERDI),     INTENT(INOUT) :: YDERDI ! Output gas concentrations

  INTEGER(KIND=JPIM) :: IYR1, IYR2, IMN0
  INTEGER(KIND=JPIM), SAVE :: IYR2OLD = 0_JPIM, IMNOLD = 0_JPIM, IYROLD = 0_JPIM
  LOGICAL :: LRDGHG

  REAL(KIND=JPRB) :: ZCO2RMWG, ZCH4RMWG, ZN2ORMWG, ZNO2RMWG, ZC11RMWG, ZC12RMWG

  INTEGER(KIND=JPIM) :: JGAS
  REAL(KIND=JPRB) :: ZFIXNO2, ZCONC(6), ZGRADC
  REAL(KIND=JPRB), SAVE :: ZFCONC(2, 5) = 0._JPRB

! Increase rate of co2 per year in 1pctCO2 experiment
  REAL(KIND=JPRB), PARAMETER :: RCO2INC = 0.01_JPRB

  ! Reference value of the global mean volume mixing ratio of CH4 at the surface,
  ! used in the parameterization of the water vapour production from methane oxidation in the stratosphere.
  ! In the CMIP6 version of EC-Earth3 a value of 1.78 ppmv was used
  ! corresponding to the reference value in the Cariolle climatology (IPCC, 2002), set in module SU_GHGCLIM.
  ! In the current model, a default value of 1.72 ppmv from IPCC (1990) is set, 
  ! but this value is "unlikely to be used in any configuration" (see module SUECRAD).
  ! We therefore use a value from the time series from the CMIP7 GHG forcing data set;
  ! and pick a recent year that is more or less consistent with the updated value of RQLIM (see below):
  ! 1.869522 ppmv for the year 2018: updated constant RQLIM operational
  ! 1.877265 ppmv for 2019
  ! 1.889479 ppmv for 2020
  ! 1.906181 ppmv for 2021: last year CMIP7 historical simulation 
  REAL(KIND=JPRB),PARAMETER :: RCH4_CL = 1.906181_JPRB 

  LOGICAL, SAVE :: FIRST_CALL = .TRUE.

  ! Generic GHG File Info for a given year to read and CMIP version
  TYPE GHGFILEINFO
     CHARACTER(LEN=255) :: DATADIR
     CHARACTER(LEN=255) :: FILEID   ! filename without leading GHG name
     INTEGER(KIND=JPIM) :: IFIRSTYR ! 1st year in the file
     INTEGER(KIND=JPIM) :: ILASTYR  ! last year in the file
  END TYPE GHGFILEINFO

  ASSOCIATE (RCARDI => YDERDI%RCARDI, RCFC11 => YDERDI%RCFC11, &
   & RCFC12 => YDERDI%RCFC12, RCH4 => YDERDI%RCH4, RN2O => YDERDI%RN2O, &
   & RNO2 => YDERDI%RNO2)

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

    ! do we need to update GHG conc ?
    LRDGHG = .FALSE.
    IF (LGHGMONTHLY) THEN
      LRDGHG = (IMN .NE. IMNOLD).OR.(IYR .NE. IYROLD)
    ELSE
      LRDGHG = (IMN .NE. IMNOLD).OR.(IYR2 .NE. IYR2OLD) ! Test on IYR2 in case we jump a full year between calls
    END IF

    IF (LRDGHG) THEN
      WRITE (NULOUT, *) 'ECE_CMIP_GHG:'
      IF (NCMIPFIXYR <= 0) THEN
        WRITE (NULOUT, FMT='('' IYR ='',I4,'' IMN ='',I4,'' IMN0 ='',I4 &
             & ,'' IYR1='',I4,'' IYR2='',I4,'' IYR2OLD='',I4,'' IMNOLD='',I4)') &
             & IYR, IMN, IMN0, IYR1, IYR2, IYR2OLD, IMNOLD
      ELSE
        WRITE (NULOUT, FMT='('' NCMIPFIXYR ='',I4,'' IMN ='',I4,'' IMN0 ='',I4 &
             & ,'' L1pctCO2='',L4,'' LANxCO2='',L4,'' RNxCO2='',F4.2)') &
             & NCMIPFIXYR, IMN, IMN0, L1pctCO2, LANxCO2, RNxCO2
      END IF

      IF (LGHGMONTHLY) THEN
        CALL READCMIPGHGDATA(IYR, IMN)
      ELSE IF ((IYR2 .NE. IYR2OLD).OR.FIRST_CALL) THEN ! FIRST_CALL is needed for corner case: starting an experiment in the first 6 months of year 0000
        IF (FIRST_CALL.OR.(IYR2OLD.NE.IYR1)) THEN      ! Test on IYR1 in case we jump more than one month between calls
          CALL READCMIPGHGDATA(IYR1)
        END IF
        CALL READCMIPGHGDATA(IYR2) ! IYR2OLD pushed into IYR1
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

      DO JGAS = 1, 6
        WRITE (NULOUT, *) 'JGAS=', JGAS, ' ZCONC(JGAS)=', ZCONC(JGAS)
      END DO

      RCARDI = ZCONC(1)*1.E-06_JPRB*ZCO2RMWG
      RCH4 = ZCONC(2)*1.E-09_JPRB*ZCH4RMWG
      RN2O = ZCONC(3)*1.E-09_JPRB*ZN2ORMWG
      RCFC11 = ZCONC(4)*1.E-12_JPRB*ZC11RMWG
      RCFC12 = ZCONC(5)*1.E-12_JPRB*ZC12RMWG
      RNO2 = ZCONC(6)*1.E-13_JPRB*ZNO2RMWG

      ! Modify RQLIM so to scale the production of water vapour from CH4 oxidation
      !   RQLIM = RQLIM * [1.0 + (2.0/7.7)*Δ[CH4]_gm(t)].
      ! The constant RQLIM was increased from 4.25E-6 to 4.81E-6 (corresponding to an increase from 6.8 to 7.7 ppmv)
      ! in march 2017; this update became operational in cycle 45r1 in 2018.
      ! Hence, the value of RCH4_CL has also been updated (see above).
      !
      ! Ideally this scaling is done based on the interpolated annual global mean mixing ratio at the surface
      ! as in the CMIP6 version of EC-Earth3
      ! Currently, it uses ZCONC(2), even if LGHGMONTHLY is TRUE
      ! This still needs to be changed.
      RQLIM = 4.81E-6_JPRB * (1.0_JPRB + (2.0_JPRB/7.7_JPRB)* (ZCONC(2) * 1.E-3_JPRB - RCH4_CL ))
      WRITE(NULOUT,*) 'Updating RQLIM=4.81e-6 * (1.+2./7.7 * Delt[CH4] )=',RQLIM

      IMNOLD = IMN
      IYROLD = IYR
    END IF

  END ASSOCIATE

  RETURN

CONTAINS

  SUBROUTINE READCMIPGHGDATA(IYEAR, IMONTH)

    USE YOMMP0, ONLY: MYPROC, NPROC
    USE MPL_MODULE, ONLY: MPL_BROADCAST

    INTEGER, PARAMETER :: RPRC = 1, ITAG = 12345

    INTEGER(KIND=JPIM), INTENT(IN)  :: IYEAR
    INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL  :: IMONTH

    REAL(KIND=JPRB), ALLOCATABLE :: ZFCONC_BUF(:)
    INTEGER, ALLOCATABLE :: START(:), STARTCH4(:)

    CHARACTER(LEN=255) :: FILENAME
    INTEGER(KIND=JPIM) :: I, JYEAR, NDIM
    INTEGER            :: IUNIT, IDX, IDX_CH4, NY_VARID

    INTEGER            :: NMON, IMON

    CHARACTER(LEN=*), PARAMETER :: NY_NAME = 'time'
    REAL(KIND=JPRB) :: ZZ_YEARS

    CHARACTER(LEN=45) :: CO2_FILE, CH4_FILE, N2O_FILE, CFC11_FILE, CFC12_FILE
    CHARACTER(LEN=45) :: CO2_NAME, CH4_NAME, N2O_NAME, CFC11_NAME, CFC12_NAME
    INTEGER(KIND=JPIM):: CO2_VARID, CH4_VARID, N2O_VARID, CFC11_VARID, CFC12_VARID
    REAL(KIND=JPRB)   :: ZZCO2, ZZCH4, ZZN2O, ZZCFC11, ZZCFC12
    TYPE(GHGFILEINFO) :: ZTGHG, ZTGHGCH4

    IF (LCMIP7) THEN
      CO2_FILE = 'co2'
      CO2_NAME = 'co2'
      CH4_FILE = 'ch4'
      CH4_NAME = 'ch4'
      N2O_FILE = 'n2o'
      N2O_NAME = 'n2o'
      CFC11_FILE = 'cfc11eq'
      CFC11_NAME = 'cfc11eq'
      CFC12_FILE = 'cfc12'
      CFC12_NAME = 'cfc12'
      NDIM=1
    ELSEIF (LCMIP6) THEN
      CO2_FILE = 'mole-fraction-of-carbon-dioxide-in-air'
      CO2_NAME = 'mole_fraction_of_carbon_dioxide_in_air'
      CH4_FILE = 'mole-fraction-of-methane-in-air'
      CH4_NAME = 'mole_fraction_of_methane_in_air'
      N2O_FILE = 'mole-fraction-of-nitrous-oxide-in-air'
      N2O_NAME = 'mole_fraction_of_nitrous_oxide_in_air'
      CFC11_FILE = 'mole-fraction-of-cfc11eq-in-air'
      CFC11_NAME = 'mole_fraction_of_cfc11eq_in_air'
      CFC12_FILE = 'mole-fraction-of-cfc12-in-air'
      CFC12_NAME = 'mole_fraction_of_cfc12_in_air'
      NDIM=2
    ENDIF
    ALLOCATE(START(NDIM))
    ALLOCATE(STARTCH4(NDIM))
    START(1)=1                  ! Matters only for CMIP6, for which its the sector index,
    STARTCH4(1)=1               ! and ignored (overwritten hereafter) for CMIP7

    IF (LGHGMONTHLY) THEN
      WRITE (NULOUT, *) 'ECE_CMIP_GHG: read monthly CMIP GHG concentrations for ', IYEAR, IMONTH
    ELSE
      WRITE (NULOUT, *) 'ECE_CMIP_GHG: read yearly CMIP GHG concentrations for ', IYEAR
    ENDIF

    ALLOCATE (ZFCONC_BUF(SIZE(ZFCONC)))

    IF (MYPROC == RPRC) THEN

      ! -- All GHG except CH4
      IF (NCMIPFIXYR .GT. 0) THEN
        JYEAR = NCMIPFIXYR
      ELSE
        JYEAR = IYEAR
      END IF
      
      ! -- GHG data only exists until 2100 
      ! -- Repeat 2100 if necessary
      IF (LCMIP7) THEN
          JYEAR = MIN(JYEAR,2100)
          WRITE(NULOUT, *) "ECE_CMIP_GHG: Set JYEAR=MIN(JYEAR,2100)=",JYEAR
      END IF

      ZTGHG = TGHG(JYEAR)

      IF (LGHGMONTHLY) THEN
        NMON = 12
        IMON = IMONTH
      ELSE
        NMON = 1
        IMON = 1
      END IF
      IF (NCMIPFIXYR .GT. 0) THEN
        IDX = (NCMIPFIXYR - ZTGHG%IFIRSTYR)*NMON + IMON
      ELSE
        IDX = (MIN(IYEAR, ZTGHG%ILASTYR) - ZTGHG%IFIRSTYR)*NMON + IMON
      END IF
      START(NDIM) = IDX

      ! -- Specific case of CH4
      IF (NCMIPFIXYR_CH4 .GT. 0) THEN
        JYEAR = NCMIPFIXYR_CH4
      ELSE
        JYEAR = IYEAR
      END IF
      
      ! -- CH4 data ends in 2100
      ! -- Repeat 2100 if necessary
      IF (LCMIP7) THEN
          JYEAR = MIN(JYEAR,2100)
          WRITE(NULOUT,*) "ECE_CMIP_GHG: Set JYEAR=MIN(JYEAR,2100)=",JYEAR 
      END IF

      ZTGHGCH4 = TGHG(JYEAR)

      IF (NCMIPFIXYR_CH4 .GT. 0) THEN
        IDX_CH4 = (NCMIPFIXYR_CH4 - ZTGHGCH4%IFIRSTYR)*NMON + IMON
      ELSE
        IDX_CH4 = (MIN(IYEAR, ZTGHGCH4%ILASTYR) - ZTGHGCH4%IFIRSTYR)*NMON + IMON
      END IF
      STARTCH4(NDIM) = IDX_CH4

      ! Read in CO2 data
      WRITE (FILENAME, '(A)') TRIM(ZTGHG%DATADIR)//'/'//TRIM(CO2_FILE)//TRIM(ZTGHG%FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))

      ! Check timestamp only from CO2 file
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, NY_NAME, NY_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, NY_VARID, ZZ_YEARS, start=(/IDX/)))
      WRITE (NULOUT, *) 'ZZ_YEARS=', ZZ_YEARS

      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CO2_NAME, CO2_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CO2_VARID, ZZCO2, start=START))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCO2=', ZZCO2

      ! Read in CH4 data
      WRITE (FILENAME, '(A)') TRIM(ZTGHG%DATADIR)//'/'//TRIM(CH4_FILE)//TRIM(ZTGHGCH4%FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CH4_NAME, CH4_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CH4_VARID, ZZCH4, start=STARTCH4))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCH4=', ZZCH4

      ! Read in N2O data
      WRITE (FILENAME, '(A)') TRIM(ZTGHG%DATADIR)//'/'//TRIM(N2O_FILE)//TRIM(ZTGHG%FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, N2O_NAME, N2O_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, N2O_VARID, ZZN2O, start=START))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZN2O=', ZZN2O

      ! Read in CFC11 data
      WRITE (FILENAME, '(A)') TRIM(ZTGHG%DATADIR)//'/'//TRIM(CFC11_FILE)//TRIM(ZTGHG%FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CFC11_NAME, CFC11_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CFC11_VARID, ZZCFC11, start=START))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCFC11=', ZZCFC11

      ! Read in CFC12 data
      WRITE (FILENAME, '(A)') TRIM(ZTGHG%DATADIR)//'/'//TRIM(CFC12_FILE)//TRIM(ZTGHG%FILEID)//'.nc'
      WRITE (NULOUT, *) TRIM(FILENAME)

      CALL CHECKGHG(NF90_OPEN(FILENAME, NF90_NOWRITE, IUNIT))
      CALL CHECKGHG(NF90_INQ_VARID(IUNIT, CFC12_NAME, CFC12_VARID))
      CALL CHECKGHG(NF90_GET_VAR(IUNIT, CFC12_VARID, ZZCFC12, start=START))
      CALL CHECKGHG(NF90_CLOSE(IUNIT))
      WRITE (NULOUT, *) 'ZZCFC12=', ZZCFC12

      ! replace IYR1 with IYR2 and fill IYR2 with new values
      ZFCONC(1, :) = ZFCONC(2, :)

      ZFCONC(2, 2) = ZZCH4
      ZFCONC(2, 3) = ZZN2O
      ZFCONC(2, 4) = ZZCFC11
      ZFCONC(2, 5) = ZZCFC12

      IF (LANxCO2) THEN
        ! AbruptNxCO2 - The CO2 abruptly multiplied by RNxCO2 of the NCMIPFIXYR level and held constant
        ZFCONC(2, 1) = ZZCO2*RNxCO2
        WRITE (NULOUT, '(''LANxCO2='',L1,'' RNxCO2='',F4.2,'' CO2 concentration prescribed as '',F8.2)') &
          &   LANxCO2, RNxCO2, ZFCONC(2, 1)

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
    CALL MPL_BROADCAST(ZFCONC_BUF, KTAG=ITAG + 3, KROOT=RPRC, CDSTRING='CMIPGHG-1: ')
    ZFCONC = RESHAPE(ZFCONC_BUF, SHAPE(ZFCONC))
    DEALLOCATE (ZFCONC_BUF)

    RETURN
  END SUBROUTINE READCMIPGHGDATA

  SUBROUTINE CHECKGHG(STATUS)
    INTEGER, INTENT(IN) :: STATUS

    IF (STATUS /= NF90_NOERR) THEN
      CALL ABOR1('ECE_CMIP_GHG: '//TRIM(NF90_STRERROR(STATUS)))
    END IF
    RETURN
  END SUBROUTINE CHECKGHG


  FUNCTION TGHG(IY)
    !
    ! Build generic (partial) GHG filename according to CMIP version and requested date
    !
    INTEGER(KIND=JPIM), INTENT(IN) :: IY
    TYPE(GHGFILEINFO) :: TGHG

    CHARACTER(LEN=13) :: TIMEPERIOD
    CHARACTER(LEN=15) :: MIP
    CHARACTER(LEN=15) :: SOURCE_ID

    IF (LCMIP7) THEN

      IF (IY < 2023) THEN
        TGHG%DATADIR = TRIM(CMIP7DATADIR)//'/ghg'
        MIP = 'CMIP'
        SOURCE_ID = 'CR-CMIP-1-0-0'
        SELECT CASE (IY)
        CASE (1:999)
          TGHG%IFIRSTYR = 1        ! The first year in the data set
          TGHG%ILASTYR = 999
        CASE (1000:1749)
          TGHG%IFIRSTYR = 1000
          TGHG%ILASTYR = 1749
        CASE (1750:2022)
          TGHG%IFIRSTYR = 1750
          TGHG%ILASTYR = 2022
        END SELECT
        ! yearly or monthly
        !IF (LGHGMONTHLY) THEN
        !  WRITE(TIMEPERIOD, '(I4.4,I2.2,"-",I4.4,I2.2)') TGHG%IFIRSTYR,1,TGHG%ILASTYR,12
        !ELSE
        !  WRITE(TIMEPERIOD, '(I4.4,"-",I4.4)') TGHG%IFIRSTYR,TGHG%ILASTYR
        !END IF
        !TGHG%FILEID = '_input4MIPs_GHGConcentrations_CMIP_CR-CMIP-1-0-0_gm_'//TRIM(TIMEPERIOD)
      
      ELSE IF (IY >= 2023 .AND. IY <= 2100) THEN 
        TGHG%DATADIR = TRIM(CMIP7DATADIR)//'/ghg/ScenarioMIP/'
        MIP = 'ScenarioMIP'
        SOURCE_ID = 'CR-'//TRIM(SCENARIONAME)//'-1-1-0'
        TGHG%IFIRSTYR = 2022
        TGHG%ILASTYR = 2100
        
      ELSE IF (IY > 2100 .AND. IY <= 2300) THEN 
        ! extended scenarios
        TGHG%DATADIR = TRIM(CMIP7DATADIR)//'/ghg/ScenarioMIP/'
        MIP = 'ScenarioMIP'
        SOURCE_ID = 'CR-'//TRIM(SCENARIONAME)//'-ext-1-1-0'
        SELECT CASE (IY)
          CASE (2101:2200)
            TGHG%IFIRSTYR = 2101
            TGHG%ILASTYR = 2200
          CASE (2201:2300)
            TGHG%IFIRSTYR = 2201
            TGHG%ILASTYR = 2300
          CASE (2301:2400)
            TGHG%IFIRSTYR = 2301
            TGHG%ILASTYR = 2400
          CASE (2401:2500)
            TGHG%IFIRSTYR = 2401
            TGHG%ILASTYR = 2500
        END SELECT 
      ENDIF

      ! yearly or monthly
      IF (LGHGMONTHLY) THEN
        WRITE(TIMEPERIOD, '(I4.4,I2.2,"-",I4.4,I2.2)') TGHG%IFIRSTYR,1,TGHG%ILASTYR,12
      ELSE
        WRITE(TIMEPERIOD, '(I4.4,"-",I4.4)') TGHG%IFIRSTYR,TGHG%ILASTYR
      END IF
      TGHG%FILEID = '_input4MIPs_GHGConcentrations_'//TRIM(MIP)//'_'//TRIM(SOURCE_ID)//'_gm_'//TRIM(TIMEPERIOD)

    ELSEIF (LCMIP6) THEN

      TGHG%DATADIR = CMIP6DATADIR

      IF (IY < 2015) THEN
        ! yearly or monthly
        IF (LGHGMONTHLY) THEN
          TIMEPERIOD = '000001-201412'
        ELSE
          TIMEPERIOD = '0000-2014'
        END IF

        TGHG%FILEID = '_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        ! The first year in the data set is 0
        TGHG%IFIRSTYR = 0
        TGHG%ILASTYR = 2014
      ELSE
        ! yearly or monthly
        IF (LGHGMONTHLY) THEN
          TIMEPERIOD = '201501-250012'
        ELSE
          TIMEPERIOD = '2015-2500'
        END IF

        SELECT CASE (TRIM(SCENARIONAME))
        CASE ("SSP1-1.9")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp119-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP1-2.6")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP2-4.5")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP3-7.0")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP3-LowNTCF")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_AerChemMIP_UoM-AIM-ssp370-lowNTCF-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP4-3.4")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-GCAM4-ssp434-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP4-6.0")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-GCAM4-ssp460-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP5-3.4-OS")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp534-over-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE ("SSP5-8.5")
          TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
        CASE DEFAULT
          !! THIS SHOULD NOT BE NEEDED ANYMORE, SINCE WE ARE NOW USING THE LAST YEAR IN
          !!    MIN(IYEAR, ZTGHG%ILASTYR)
          !! TO FIND IDX IN THE CALLING ROUTINE. I.E. WE REPEAT 2014 FOR 2015 WHEN RUNNING HISTORICAL.
          !! COMMENT AS OBSOLETE FOR NOW, REMOVE LATER:

          !OBSOLETE   ! default to SSP3-7.0 but only to finish historical runs (last step of historical)
          !OBSOLETE   IF (.NOT. FIRST_CALL) THEN
          !OBSOLETE     TGHG%FILEID = '_input4MIPs_GHGConcentrations_ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr1-GMNHSH_'//TRIM(TIMEPERIOD)
          !OBSOLETE   ELSE
          CALL ABOR1('ECE_CMIP_GHG : unknown CMIP6 scenario')
          !OBSOLETE   END IF
        END SELECT

        ! The first year in the data set is 2015
        TGHG%IFIRSTYR = 2015
        TGHG%ILASTYR = 2500
      END IF
    ENDIF

  END FUNCTION TGHG

END SUBROUTINE ECE_CMIP_GHG
