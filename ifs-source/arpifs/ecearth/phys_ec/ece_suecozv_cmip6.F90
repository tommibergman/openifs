SUBROUTINE ECE_SUECOZV_CMIP6(KINDAT, KMINUT)

!**** *ECE_SUECOZV_CMIP6* - GETS VARIATIONAL DISTRIBUTION OF OZONE
!                 (FOR CMIP6 RUNS)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *ECE_SUECOZV_CMIP6* FROM *UPDTIM*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
!     ==== OUTPUTS ===
! ROZT1   :                : AMOUNT OF OZONE (KG/KG)

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     Shuting Yang  EC-Earth 9/02/2010
!     (after SUECOZC, J.-J. MORCRETTE  E.C.M.W.F. 95/01/25)

!     MODIFICATIONS.
!     --------------
!     S. Stefanescu - 14/10/2010 Update input fields when year changes;
!                     set IYR to 2099 after 2099 (instead of 2100),
!                     corresponding to the last ozone scenario file
!     S. Stefanescu - 19/10/2010 Correct weights for time interpolation
!     Michiel van Weele - 16 aug 2016 CMIP6 historical ozone forcings
!     Twan van Noije - May 2022, adapted from EC-Earth3
!     Klaus Wyser - June 2023, minor tweaks for integration in ECE4

!---------------------------------------------------------------------

  USE PARKIND1, ONLY: JPIM, JPIB, JPRB, JPRD
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMLUN, ONLY: NULOUT
  USE YOMCST, ONLY: RPI, RDAY
  USE YOERDI, ONLY: YRERDI
  USE ECE_YOEOZOV_CMIP6, ONLY: NLON1, NLAT1, NLV1, RLATCLI, RLONCLI, RLONCLI_DEG, &
       & RSINC1, ROZT1, RPROC1, NMONTH1
  USE YOERAD, ONLY: YRERAD

  USE YOMMP0, ONLY: MYPROC
  USE MPL_MODULE, ONLY: MPL_BROADCAST

  USE NETCDF

  USE ECE_CMIP6, ONLY: CMIP6DATADIR, SSPNAME, NCMIPFIXYR

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(IN)    :: KINDAT
  INTEGER(KIND=JPIB), INTENT(IN)    :: KMINUT
!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------

!*       0.2   LOCAL ARRAYS.
!              -------------

  REAL(KIND=JPRB), SAVE :: ZOZCL(NLON1, NLAT1, NLV1, 0:13)

  REAL(KIND=JPRB) :: ZOZCL_BUF(NLON1*NLAT1*NLV1*14)
  INTEGER, PARAMETER :: RPRC = 1, ITAG = 34567

  INTEGER(KIND=JPIM) :: IDY, IM, IM1, IM2, IMN, JK, JL, JI, JM
  INTEGER(KIND=JPIM) :: IH0, IJ0, IM0, IA0, IDD, ISS, IHR, IMIN, ISC, IYR, ILMOIS(12)
  INTEGER(KIND=JPIM) :: IREADY, IREADM

  INTEGER(KIND=JPIM) :: I, IUNIT, IDIR, IFIL, IRT
  LOGICAL            :: IS_OPEN
  REAL(KIND=JPRB)    :: SKIP
  CHARACTER(LEN=132) :: SKIP_LINE

  CHARACTER(LEN=200) ::  CLFN
  CHARACTER(LEN=100) ::  CO3DATAFIL

  REAL(KIND=JPRB) :: ZTIMI, ZXTIME, ZINCH, ZINCL, ZINCROZ
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  INTEGER(KIND=JPIM), SAVE :: IYROLD = -999, IM1OLD = -999

#include "abor1.intfb.h"
#include "updcalsec.intfb.h"

#include "fcttim.func.h"
!     ------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('ECE_SUECOZV_CMIP6', 0, ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*         1.     TIME INDEX WITHIN OZONE CLIMATOLOGY
!                 -----------------------------------

!ECEARTH: Time interpolation of climatology for any type of calendar. (See also suecaec.F90!)

  IA0 = NCCAA(KINDAT)
  IM0 = NMM(KINDAT)
  IJ0 = NDD(KINDAT)
  IH0 = 0
  IDD = KMINUT/1440_JPIM ! The number of days since KINDAT
  ISS = Modulo(KMINUT, 1440_JPIM)*60_JPIM ! The number of seconds since the start of the current day
  CALL UPDCALSEC(IH0, IJ0, IM0, IA0, IDD, ISS, IHR, IMIN, ISC, IDY, IMN, IYR, ILMOIS, -1)

  ZXTIME = IDY - 1 + (60*(60*IHR + IMIN) + ISC)/86400.0_JPRB - 0.5*ILMOIS(IMN) ! Number of days relative to center of the current month
  IF (ZXTIME < 0.0_JPRB) THEN
    IM1 = Modulo(IMN - 2, 12) + 1
    IM2 = IMN
    ZXTIME = ZXTIME + 0.5*(ILMOIS(IM1) + ILMOIS(IM2)) ! Adjust relative time to the center of the previous month
  ELSE
    IM1 = IMN
    IM2 = Modulo(IMN, 12) + 1
  END IF
  ZTIMI = ZXTIME/(0.5*(ILMOIS(IM1) + ILMOIS(IM2))) ! Compute interpolation weight between the centers of the two months

! Yearly files include dec data from previous year (index 0), and jan data from next year (index
!  13). Use those if we are not in a perpetual forcing.
  IF (NCMIPFIXYR <= 0) THEN
    IF (IMN == 1 .AND. IM1 == 12) THEN
      IM1 = 0
    ELSEIF (IMN == 12 .AND. IM2 == 1) THEN
      IM2 = 13
    END IF
  END IF

!*         1.     TIME INDEX WITHIN OZONE CLIMATOLOGY
!                 -----------------------

! SET TIME INTERVAL

! use perpetual year if NCMIPFIXYR is set
  IF (NCMIPFIXYR > 0) IYR = NCMIPFIXYR

! limit IYR to available dataset
  IYR = MIN(2099, MAX(1850, IYR))

  IF (IYR /= IYROLD) THEN

    ! READ OZONE FORCING FILE
    ! Original CMIP6 data have been processed into yearly files containing NMONTH=14 months
    ! e.g. cmip6-data/o3_histo/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_1850.nc
    ! Store in ZOZCL(JI,JL,JK,JM) with CMIP6 ozone forcing data dimensions (144,96,66,14)
    IF (NCMIPFIXYR == 1850) THEN
      WRITE (CO3DATAFIL, '(''o3_pi/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_clim_'',I4.4,''.nc'')') IYR
    ELSEIF (IYR <= 2014) THEN
      WRITE (CO3DATAFIL, '(''o3_histo/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_'',I4.4,''.nc'')') IYR
    ELSE ! IYR < 2100
      SELECT CASE (TRIM(SSPNAME))
      CASE ("SSP1-1.9", "SSP1-2.6", "SSP1-2.6-Ext")
        ! For SSP1-1.9, use ozone from SSP1-2.6
        WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp126_'',I4.4,''.nc'')') IYR
      CASE ("SSP2-4.5")
        WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp245_'',I4.4,''.nc'')') IYR
      CASE ("SSP3-7.0")
        WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp370_'',I4.4,''.nc'')') IYR
      CASE ("SSP3-LowNTCF")
        CALL ABOR1('ECE_SUECOZV_CMIP6: This scenario is to be used for AerChemMIP only, which does not use this forcing.')
      CASE ("SSP4-3.4")
        WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp434_'',I4.4,''.nc'')') IYR
      CASE ("SSP5-3.4-OS")
        WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp534os_'',I4.4,''.nc'')') IYR
      CASE ("SSP4-6.0", "SSP5-3.4-OS-Ext")
        CALL ABOR1('ECE_SUECOZV_CMIP6: No data provided for this Tier 2 scenario')
      CASE ("SSP5-8.5", "SSP5-8.5-Ext")
        WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp585_'',I4.4,''.nc'')') IYR
      CASE DEFAULT
        ! default to SSP3-7.0 but only for the first time step of 2015
        IF (IYROLD == 2014) THEN
          WRITE (CO3DATAFIL, '(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_ssp370_'',I4.4,''.nc'')') IYR
          IYROLD = -1
        ELSE
          CALL ABOR1('ECE_SUECOZV_CMIP6: Unknown scenario')
        END IF
      END SELECT
    END IF

    IF (IYROLD /= -1) IYROLD = IYR

    IF (MYPROC == RPRC) THEN
      IDIR = LEN_TRIM(CMIP6DATADIR)
      IFIL = LEN_TRIM(CO3DATAFIL)
      CLFN = CMIP6DATADIR(1:IDIR)//'/'//CO3DATAFIL(1:IFIL)
      WRITE (NULOUT, '("ECE_SUECOZV_CMIP6: READ IN CMIP6 OZONE DATA FROM FILE ",A)') TRIM(CLFN)
      CALL READ_NC_FILE_OZONE_CMIP6(CLFN, IRT)

      ZOZCL_BUF = RESHAPE(ZOZCL, SHAPE(ZOZCL_BUF))
    END IF

    CALL MPL_BROADCAST(IRT, KTAG=ITAG, KROOT=RPRC, CDSTRING='SUECOZV-0: ')
    CALL CHECK(IRT)

    ! distribute to all PROCs
    CALL MPL_BROADCAST(ZOZCL_BUF, KTAG=ITAG, KROOT=RPRC, CDSTRING='SUECOZV-1: ')
    ZOZCL = RESHAPE(ZOZCL_BUF, SHAPE(ZOZCL))

    DO JI = 1, NLON1
      RLONCLI(JI) = RLONCLI_DEG(JI)*RPI/180._JPRB
    END DO

    IM1OLD = -999
  END IF

!*         2.0    TIME INTERPOLATED FIELD
!                 -----------------------

  IF (IM1 /= IM1OLD) THEN
    IM1OLD = IM1
!*( Field is also transformed in kg/kg! )
    DO JK = 1, NLV1
      DO JL = 1, NLAT1
        DO JI = 1, NLON1
          ROZT1(JI, JL, JK) = YRERDI%RO3*(ZOZCL(JI, JL, JK, IM1)&
            & + ZTIMI*(ZOZCL(JI, JL, JK, IM2) - ZOZCL(JI, JL, JK, IM1)))
        END DO
      END DO
    END DO

    DO JL = 1, NLAT1
      DO JI = 1, NLON1
        ROZT1(JI, JL, 0) = 0.0_JPRB
        ROZT1(JI, JL, NLV1 + 1) = ROZT1(JI, JL, NLV1)
      END DO
    END DO

!     VECTOR OF LATITUDES FOR OZONE CLIMATOLOGY:

    DO JL = 1, NLAT1
      RSINC1(JL) = SIN(RLATCLI(JL)*RPI/180.0_JPRB)
    END DO
  END IF

  IF (LHOOK) CALL DR_HOOK('ECE_SUECOZV_CMIP6', 1, ZHOOK_HANDLE)

CONTAINS

  SUBROUTINE READ_NC_FILE_OZONE_CMIP6(CC, IRC)

    CHARACTER(LEN=200), INTENT(IN)  :: CC
    INTEGER(KIND=JPIM), INTENT(OUT) :: IRC

    ! use precision of vmro3 in the file
    REAL(KIND=4) :: OZO_CMIP6(NLON1, NLAT1, NLV1, 1)

    INTEGER(KIND=JPIM)               :: INCUNIT, OZO_VARID
    INTEGER(KIND=JPIM)               :: IMONTH1, ILV
    INTEGER(KIND=JPIM), DIMENSION(4) :: ISTART, ISIZE
    CHARACTER(LEN=*), PARAMETER      :: OZO_NAME = 'vmro3'

    ! OPEN NETCDF FILE
    IRC = NF90_OPEN(CC, NF90_NOWRITE, INCUNIT)
    IF (IRC /= NF90_NOERR) RETURN

    ! READ 3-D FIELD PER MONTH (14 MONTHS in total)
    IRC = NF90_INQ_VARID(INCUNIT, OZO_NAME, OZO_VARID)
    IF (IRC /= NF90_NOERR) RETURN

    ISIZE = (/NLON1, NLAT1, NLV1, 1/)

    DO IMONTH1 = 1, NMONTH1
      ISTART = (/1, 1, 1, IMONTH1/)

      IRC = NF90_GET_VAR(INCUNIT, OZO_VARID, OZO_CMIP6, start=ISTART, count=ISIZE)
      IF (IRC /= NF90_NOERR) RETURN

      ! store monthly ozone field in ZOZCL and convert from mole/mole to ppm
      ! reverse CMIP6 ozone inputdata (from surf to TOA) to (from TOA to surf)
      DO ILV = 1, NLV1
        ZOZCL(:, :, NLV1 - ILV + 1, IMONTH1 - 1) = OZO_CMIP6(:, :, ILV, 1)*1.E+06_JPRB
      END DO

    END DO

    ! CLOSE NETCDF FILE
    IRC = NF90_CLOSE(INCUNIT)

  END SUBROUTINE READ_NC_FILE_OZONE_CMIP6

  SUBROUTINE CHECK(STATUS)

    INTEGER, INTENT(IN) :: STATUS

    IF (STATUS /= NF90_NOERR) THEN
      CALL ABOR1('READ_NC_FILE_OZONE_CMIP6: '//TRIM(NF90_STRERROR(STATUS)))
    END IF

  END SUBROUTINE CHECK

END SUBROUTINE ECE_SUECOZV_CMIP6
