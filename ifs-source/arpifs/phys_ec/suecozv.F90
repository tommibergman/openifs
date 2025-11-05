! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
!pgi$r opt=0 
SUBROUTINE SUECOZV(YDECMIP,KINDAT)

!**** *SUECOZV* - GETS VARIATIONAL DISTRIBUTION OF OZONE
!                 (FOR CMIP RUNS)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SUECOZV* FROM *UPDTIM*

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
!     Hans Hersbach ECMWF 09/09/2010 Update input fields at year changes
!     P. Bechtold         14/05/2012 replace 86400 by RDAY
!     R. Senan/C. Roberts 26/01/2017 Support for CMIP6 forcings
!     O. Marsden          30/01/2018 Split the update out from the setup, new UPDECOZV routine
!     R. Senan            15/12/2021 CMIP6: Support for single precision
!     J. Kjellsson        25/09/2025 Support for CMIP7 forcings 
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD, JPIB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RPI, RDAY
USE YOECMIP  , ONLY : TECMIP,NLON1_CMIP5, NLAT1_CMIP5, NLV1_CMIP5, NMONTH1, &
 &                    NLON1_CMIP6, NLAT1_CMIP6, NLV1_CMIP6, & 
 &                    NLON1_CMIP7, NLAT1_CMIP7, NLV1_CMIP7 
! &                    CO3DATADIR, CO3DATAFIL, NRCP, NO3CMIP,NCMIPFIXYR
USE ECE_CMIP,  ONLY : NCMIPFIXYR, SCENARIONAME 

IMPLICIT NONE

!     -----------------------------------------------------------------

!*       0.1   ARGUMENTS.
!              ----------
TYPE(TECMIP)      ,INTENT(INOUT) :: YDECMIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDAT
!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------

INTEGER(KIND=JPIM) :: JK, JL, JI, JM, IYR

INTEGER(KIND=JPIM) :: I, IUNIT, IDIR, IFIL
LOGICAL            :: LLIS_OPEN
CHARACTER(LEN=132) :: CLSKIP_LINE

CHARACTER (LEN = 300) ::  CLFN       ! full file name (dir + file)
CHARACTER (LEN = 80)  ::  ZO3DATAFIL ! temporary string for file name


CHARACTER (LEN = 10)  ::  CO3SCEN  ! scenario name 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL(KIND=JPRB), ALLOCATABLE :: ZOZO_DATA(:,:,:,:) ! temporary array for ozone 

INTEGER(KIND=JPIM) :: NLON1, NLAT1 ,NLV1
INTEGER(KIND=JPIM) :: IYEAR1, IYEAR2 

LOGICAL            :: LFIRSTYEAR  ! if first year of the ozone file 
LOGICAL            :: LLASTYEAR   ! if last year of the ozone file 

INTEGER(KIND=JPIM),SAVE :: IYROLD=-999
#include "abor1.intfb.h"

#include "fcttim.func.h"      
#include "netcdf.inc"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUECOZV',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IYR = NCCAA(KINDAT)

!ECEARTH: Time interpolation of climatology for any type of calendar. (See also suecaec.F90!)


!*         1.     TIME INDEX WITHIN OZONE CLIMATOLOGY
!                 -----------------------


! SET TIME INTERVAL
IF (YDECMIP%NO3CMIP == 7) THEN ! CMIP7
  
  IF (IYR < 1850) THEN
      WRITE(NULOUT,*) "SUECOZV: For year < 1850 we set year = 1850" 
      IYR=1850
  ENDIF

  IF (IYR > 2022) THEN
      WRITE(NULOUT,*) "SUECOZV: For year > 2022 we set year = 2022" 
      IYR=2022
  ENDIF

  NLON1 = NLON1_CMIP7
  NLAT1 = NLAT1_CMIP7
  NLV1  = NLV1_CMIP7

ELSE IF (YDECMIP%NO3CMIP == 6) THEN ! CMIP6 
  IF(IYR < 1850) THEN
    IYR=1850
  ENDIF

  NLON1=NLON1_CMIP6
  NLAT1=NLAT1_CMIP6
  NLV1=NLV1_CMIP6
  
ELSEIF (YDECMIP%NO3CMIP == 5) THEN ! CMIP5
  IF(IYR < 1850) THEN 
   IYR=1850
  ELSEIF(IYR >= 2100) THEN
    IYR=2100
  ENDIF

  NLON1=NLON1_CMIP5
  NLAT1=NLAT1_CMIP5
  NLV1=NLV1_CMIP5

ELSE
  WRITE(NULOUT,*)"NO3CMIP:",YDECMIP%NO3CMIP
  CALL ABOR1('SUECOZV: Value of NO3CMIP not supported')
ENDIF

IF (.NOT. ALLOCATED(YDECMIP%RSINC1))    ALLOCATE(YDECMIP%RSINC1(NLAT1))
IF (.NOT. ALLOCATED(YDECMIP%ROZT1))     ALLOCATE(YDECMIP%ROZT1(NLON1,NLAT1,0:NLV1+1))
IF (.NOT. ALLOCATED(YDECMIP%RPROC1))    ALLOCATE(YDECMIP%RPROC1(0:NLV1+1))
IF (.NOT. ALLOCATED(YDECMIP%RLATCLI))   ALLOCATE(YDECMIP%RLATCLI(NLAT1))
IF (.NOT. ALLOCATED(YDECMIP%RLONCLI))   ALLOCATE(YDECMIP%RLONCLI(NLON1))

IF (.NOT. ALLOCATED(YDECMIP%ZOZCL))     ALLOCATE(YDECMIP%ZOZCL(NLON1,NLAT1,NLV1,0:NMONTH1-1))
  
!IF(IYR/=IYROLD) THEN
!  IYROLD=IYR


! OPEN OZONE FORCING FILE
IF (YDECMIP%NO3CMIP == 7) THEN ! read CMIP7 ozone data
  
  ! Note: Using CMIPFIXYR from NAMECECMIP rather than YDECMIP so that
  ! we dont have to set the variable twice in EC-Earth
  IF (NCMIPFIXYR > 0) THEN
     IYR = NCMIPFIXYR
     WRITE(NULOUT,*) "SUECOZV: NCMIPFIXYR to IYR ",NCMIPFIXYR,IYR
  ENDIF

  WRITE(NULOUT,*)"SUECOZV: CO3DATADIR = ",YDECMIP%CO3DATADIR 

  ! For a given year and scenario, find the right file to read
  ! and its start and end year
  ! Returns file name and whether the year is the first or last of the file
  ! Set CMIP7_SCEN=SCENARIONAME to read a scenario
  CALL FIND_NC_FILE_OZONE_CMIP7( IYR, IYEAR1, IYEAR2, &
                               & YDECMIP%CO3DATAFIL, & ! file name
                               & LFIRSTYEAR, LLASTYEAR, & ! first and last year of file 
                               & CCMIP7_SCEN=SCENARIONAME) 

  IDIR=LEN_TRIM(YDECMIP%CO3DATADIR)
  IFIL=LEN_TRIM(YDECMIP%CO3DATAFIL)
  CLFN=YDECMIP%CO3DATADIR(1:IDIR)//'/'//YDECMIP%CO3DATAFIL(1:IFIL)
  
  WRITE(NULOUT,'("SUECOZV: READ IN CMIP7 OZONE DATA FROM FILE ",A)')CLFN
  WRITE(NULOUT,'("SUECOZV: READ IN CMIP7 OZONE DATA DIMENSIONS ",4I4)') NLON1,NLAT1,NLV1,NMONTH1
  
  ! Read data from netCDF file. 
  ! We will read the last month of previous year and first month of next year, 
  ! e.g. for 2000 we read Dec 1999 - Jan 2001 as index 0:13. 
  ! But if it is the first or last year of the file, this is not possible so we just 
  ! repeat the first or last month. 
  CALL READ_NC_FILE_OZONE_CMIP7(CLFN, NLON1, NLAT1, NLV1, IYR, IYEAR1, IYEAR2, NMONTH1, YDECMIP%ZOZCL)
  
  ! If we are reading the first year of the file, we also need the last month
  ! of the previous file
  IF (LFIRSTYEAR) THEN
      WRITE(NULOUT,*) "SUECOZV: Reading previous year ",IYR-1
      ALLOCATE(ZOZO_DATA(NLON1, NLAT1, NLV1, 0:NMONTH1-1))
      CALL FIND_NC_FILE_OZONE_CMIP7(IYR-1, IYEAR1, IYEAR2, ZO3DATAFIL, LFIRSTYEAR, LLASTYEAR)
      IFIL=LEN_TRIM(ZO3DATAFIL) 
      CLFN=YDECMIP%CO3DATADIR(1:IDIR)//'/'//ZO3DATAFIL(1:IFIL)  
      CALL READ_NC_FILE_OZONE_CMIP7(CLFN, NLON1, NLAT1, NLV1, IYR-1, IYEAR1, IYEAR2, NMONTH1, ZOZO_DATA)
      ! Put time index 12 (Dec) of previous year at index 0 (Dec)
      YDECMIP%ZOZCL(:,:,:,12) = ZOZO_DATA(:,:,:,0)
      DEALLOCATE(ZOZO_DATA)  
      
  ! If we are reading the last year of the file, we also need the first month 
  ! of the next file
  ELSE IF (LLASTYEAR) THEN
      ! if year = 2022, then there is no more data
      ! We will simply repeat 2022 again. 
      ! When scenarios become available we can read that instead
      IF (IYR == 2022) THEN
          WRITE(NULOUT,*) "SUECOZV: Last year of data. Not reading next year." 
      ELSE    
          WRITE(NULOUT,*) "SUECOZV: Reading next year ",IYR+1 
          ALLOCATE(ZOZO_DATA(NLON1, NLAT1, NLV1, 0:NMONTH1-1))
          CALL FIND_NC_FILE_OZONE_CMIP7(IYR+1, IYEAR1, IYEAR2, ZO3DATAFIL, LFIRSTYEAR, LLASTYEAR)
          IFIL=LEN_TRIM(ZO3DATAFIL)
          CLFN=YDECMIP%CO3DATADIR(1:IDIR)//'/'//ZO3DATAFIL(1:IFIL)
          CALL READ_NC_FILE_OZONE_CMIP7(CLFN, NLON1, NLAT1, NLV1, IYR+1, IYEAR1, IYEAR2, NMONTH1, ZOZO_DATA)
          ! Put tme index 1 (Jan) of next year at index 13 (Jan) 
          YDECMIP%ZOZCL(:,:,:,1) = ZOZO_DATA(:,:,:,13)
          DEALLOCATE(ZOZO_DATA)
      ENDIF
  ENDIF
  
  ! end if CMIP7

ELSE IF (YDECMIP%NO3CMIP == 6) THEN ! READ CMIP6 OZONE DATA
  ! Joakim: Use CMIPFIXYR (from NAMECECMIP, rather than YDECMIP)
  IF (NCMIPFIXYR>0) THEN 
    IYR=NCMIPFIXYR ! If using perpetual CMIP forcing
    WRITE(NULOUT,*) "SUECOZV: NCMIPFIXYR to IYR ",NCMIPFIXYR,IYR
  ENDIF
  
  WRITE(NULOUT,*)"CO3DATADIR:",YDECMIP%CO3DATADIR

  IF (IYR < 1850 .OR. IYR >= 2100) THEN
    ! CMIP6 scenario forcing data not yet available for 2015+ (Feb 2017)
    WRITE(NULOUT,'("SUECOZV: UNABLE TO OPEN CMIP6 OZONE FORCING FILE ",A)')CLFN
    CALL ABOR1('SUECOZV: UNABLE TO OPEN CMIP6 OZONE FORCING FILE')
  ENDIF

  IF (IYR == 1850) THEN

    WRITE(YDECMIP%CO3DATAFIL,'(''o3_pi/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_clim_'',I4.4,''.nc'')') IYR

  ELSE IF (IYR <= 2014) THEN
    
    ! historical 
    WRITE(YDECMIP%CO3DATAFIL,'(''o3_histo/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_'',I4.4,''.nc'')') IYR
  
  ELSE IF (IYR >= 2015) THEN
    
    ! scenario
    SELECT CASE (TRIM(SCENARIONAME))
    CASE("SSP1-2.6")
       CO3SCEN="ssp126"
    CASE("SSP2-4.5")
       CO3SCEN="ssp245"
    CASE("SSP3-7.0")
       CO3SCEN="ssp370"
    CASE("SSP5-8.5")
       CO3SCEN="ssp585" 
    CASE DEFAULT
       CALL ABOR1("SUECOZV: SCENARIONAME FOR CMIP6 MUST BE SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5")
    END SELECT 
    WRITE(YDECMIP%CO3DATAFIL,'(''o3_scenarios/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_'',A,''_'',I4.4,''.nc'')') TRIM(CO3SCEN),IYR

  ENDIF  

  IDIR=LEN_TRIM(YDECMIP%CO3DATADIR)
  IFIL=LEN_TRIM(YDECMIP%CO3DATAFIL)
  CLFN=YDECMIP%CO3DATADIR(1:IDIR)//'/'//YDECMIP%CO3DATAFIL(1:IFIL)
  WRITE(NULOUT,'("SUECOZV: READ IN CMIP6 OZONE DATA FROM FILE ",A)')CLFN
  WRITE(NULOUT,'("SUECOZV: READ IN CMIP6 OZONE DATA DIMENSIONS ",4I4)') NLON1,NLAT1,NLV1,NMONTH1
  CALL READ_NC_FILE_OZONE_CMIP6(CLFN,NLON1,NLAT1,NLV1,NMONTH1,YDECMIP%ZOZCL)  


  CALL CMIP6_OZONE_COORD
 
ELSE IF (YDECMIP%NO3CMIP == 5) THEN ! READ CMIP5 OZONE DATA

  IF (YDECMIP%NCMIPFIXYR>0) IYR=YDECMIP%NCMIPFIXYR ! If using perpetual CMIP forcing

  !Find a free unit number
  DO I=10,1000
    INQUIRE(UNIT=I,OPENED=LLIS_OPEN)
    IF (.NOT.LLIS_OPEN) THEN
      IUNIT = I
      EXIT
    ENDIF
  ENDDO

  IF (IYR < 2010) THEN  
    WRITE(YDECMIP%CO3DATAFIL,'(''Ozone_CMIP5_ACC_SPARC_'',I4.4,''.dat'')') IYR
  ENDIF
  
  IF (IYR >= 2010) THEN
    IF (YDECMIP%NRCP == 1) THEN
      WRITE(YDECMIP%CO3DATAFIL,'(''Ozone_CMIP5_ACC_SPARC_'',I4.4,''_RCP2.6.dat'')') IYR
    ELSEIF (YDECMIP%NRCP == 2) THEN  
      WRITE(YDECMIP%CO3DATAFIL,'(''Ozone_CMIP5_ACC_SPARC_'',I4.4,''_RCP4.5.dat'')') IYR
    ELSEIF (YDECMIP%NRCP == 3) THEN
      WRITE(NULOUT,'("SUECOZV: CMIP5 RCP 6.0 OZONE DATA IS NOT AVAILABLE ",A)')CLFN
    ELSEIF (YDECMIP%NRCP == 4) THEN
      WRITE(YDECMIP%CO3DATAFIL,'(''Ozone_CMIP5_ACC_SPARC_'',I4.4,''_RCP8.5.dat'')') IYR    
    ENDIF
  ENDIF

  IDIR=LEN_TRIM(YDECMIP%CO3DATADIR)
  IFIL=LEN_TRIM(YDECMIP%CO3DATAFIL)
  CLFN=YDECMIP%CO3DATADIR(1:IDIR)//YDECMIP%CO3DATAFIL(1:IFIL)
  OPEN(UNIT=IUNIT,FILE=CLFN,STATUS='OLD',FORM='FORMATTED',ERR=999)
  GOTO 1000
  999 CONTINUE
  WRITE(NULOUT,'("SUECOZV: UNABLE TO OPEN FILE ",A)')CLFN
  CALL ABOR1('SUECOZV: UNABLE TO OPEN OZONE FORCING FILE')
  1000 CONTINUE
  WRITE(NULOUT,'("SUECOZV: READ IN CMIP5 OZONE DATA FROM FILE ",A)')CLFN

  ! Read in 14 months data
  DO JM=0, 13
    READ(UNIT=IUNIT,FMT=*) CLSKIP_LINE
    DO JK=1,NLV1
      READ(UNIT=IUNIT,FMT=*) CLSKIP_LINE
      DO JL=1,NLAT1
        READ(UNIT=IUNIT,FMT=*) (YDECMIP%ZOZCL(JI,JL,JK,JM),JI=1,NLON1)
      ENDDO
    ENDDO
  ENDDO
  CLOSE(UNIT=IUNIT)

  CALL CMIP5_OZONE_COORD
  !WRITE(NULOUT,'("CMIP5 OZONE: ZOZCL(JI=1,JL=19,JK=2,JM=12)=",F10.5)') YDECMIP%ZOZCL(1,19,2,12)

ELSE
  WRITE(NULOUT,*)"NO3CMIP:",YDECMIP%NO3CMIP
  CALL ABOR1('SUECOZV: Value of NO3CMIP not recognized.')
ENDIF

! store the year we've read in to ZOZCL in the TECMIP instance (used in UPDECOZV)
YDECMIP%NCURRYR = IYR


DO JI=1,NLON1
  YDECMIP%RLONCLI(JI)=YDECMIP%RLONCLI(JI)*RPI/180._JPRB
ENDDO
!ENDIF  !! IF IYR /= IYROLD


!     VECTOR OF LATITUDES FOR OZONE CLIMATOLOGY:

DO JL=1,NLAT1
  YDECMIP%RSINC1(JL)=SIN(YDECMIP%RLATCLI(JL)*RPI/180.0_JPRB) 
ENDDO

IF (LHOOK) CALL DR_HOOK('SUECOZV',1,ZHOOK_HANDLE)

!-----------------------------------------------------------------------------

CONTAINS


SUBROUTINE CMIP5_OZONE_COORD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SUECOZV:CMIP5_OZONE_COORD',0,ZHOOK_HANDLE)

YDECMIP%RPROC1 = (/&
&    0.0_JPRB, &
&   100._JPRB,  150._JPRB,    200._JPRB,   300._JPRB,   500._JPRB, &
&   700._JPRB, 1000._JPRB,   1500._JPRB,  2000._JPRB,  3000._JPRB, &
&  5000._JPRB, 7000._JPRB,   8000._JPRB, 10000._JPRB, 15000._JPRB, &
& 20000._JPRB, 25000._JPRB, 30000._JPRB, 40000._JPRB, 50000._JPRB, &
& 60000._JPRB, 70000._JPRB, 85000._JPRB,100000._JPRB, &
&110000._JPRB &
& /)

YDECMIP%RLATCLI= (/&
& -90.0, -85.0, -80.0, -75.0, -70.0, -65.0, -60.0, -55.0, &
& -50.0, -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, &
& -10.0,  -5.0,   0.0,   5.0,  10.0,  15.0,  20.0,  25.0, &
&  30.0,  35.0,  40.0,  45.0,  50.0,  55.0,  60.0,  65.0, &
&  70.0,  75.0,  80.0,  85.0,  90.0 /)

YDECMIP%RLONCLI=(/&
&   0.0,   5.0,  10.0,  15.0,  20.0,  25.0,  30.0,  35.0, &
&  40.0,  45.0,  50.0,  55.0,  60.0,  65.0,  70.0,  75.0, &
&  80.0,  85.0,  90.0,  95.0, 100.0, 105.0, 110.0, 115.0, &
& 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, &
& 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, &
& 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, &
& 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, &
& 280.0, 285.0, 290.0, 295.0, 300.0, 305.0, 310.0, 315.0, &
& 320.0, 325.0, 330.0, 335.0, 340.0, 345.0, 350.0, 355.0 /)

IF (LHOOK) CALL DR_HOOK('SUECOZV:CMIP5_OZONE_COORD',1,ZHOOK_HANDLE)
END SUBROUTINE CMIP5_OZONE_COORD


SUBROUTINE CMIP6_OZONE_COORD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SUECOZV:CMIP6_OZONE_COORD',0,ZHOOK_HANDLE)

YDECMIP%RPROC1 = (/&
&      0.0_JPRB, &
&     0.01_JPRB,    0.03_JPRB,   0.05_JPRB,   0.08_JPRB,      0.1_JPRB, &
&     0.15_JPRB,     0.2_JPRB,    0.3_JPRB,    0.4_JPRB,      0.5_JPRB, &
&      0.7_JPRB,      1._JPRB,    1.5_JPRB,     2._JPRB,       3._JPRB, &
&       4._JPRB,      5._JPRB,     7._JPRB,    10._JPRB,      15._JPRB, &
&      20._JPRB,     30._JPRB,    40._JPRB,    50._JPRB,      70._JPRB, &
&     100._JPRB,    150._JPRB,   200._JPRB,   300._JPRB,     400._JPRB, &
&     500._JPRB,    700._JPRB,  1000._JPRB,  1500._JPRB,    2000._JPRB, &
&    2500._JPRB,   3000._JPRB,  3500._JPRB,  4000._JPRB,    5000._JPRB, &
&    6000._JPRB,   7000._JPRB,  8000._JPRB,  9000._JPRB,   10000._JPRB, &
&   11500._JPRB,  13000._JPRB, 15000._JPRB, 17000._JPRB,   20000._JPRB, &
&   25000._JPRB,  28500._JPRB, 30000._JPRB, 35000._JPRB,   40000._JPRB, &
&   45000._JPRB,  50000._JPRB, 60000._JPRB, 65000._JPRB,   70000._JPRB, &
&   75000._JPRB,  78000._JPRB, 80000._JPRB, 85000._JPRB,   92500._JPRB, &
&  100000._JPRB, 110000._JPRB /)

YDECMIP%RLATCLI= (/&
&    -90.00000, -88.10526,   -86.21053, -84.31579, -82.42105, &
&    -80.52631, -78.63158,   -76.73684, -74.84210, -72.94736, &
&    -71.05264, -69.15790,   -67.26316, -65.36842, -63.47368, &
&    -61.57895, -59.68421,   -57.78947, -55.89474, -54.00000, & 
&    -52.10526, -50.21053,   -48.31579, -46.42105, -44.52632, &
&    -42.63158, -40.73684,   -38.84211, -36.94737, -35.05263, &
&    -33.15789, -31.26316,   -29.36842, -27.47368, -25.57895, &
&    -23.68421, -21.78947,   -19.89474, -18.00000, -16.10526, &
&    -14.21053, -12.31579,   -10.42105, -8.526316, -6.631579, &
&    -4.736842, -2.842105,  -0.9473684, 0.9473684,  2.842105, &
&     4.736842,  6.631579,    8.526316,  10.42105,  12.31579, &
&     14.21053,  16.10526,    18.00000,  19.89474,  21.78947, &
&     23.68421,  25.57895,    27.47368,  29.36842,  31.26316, &
&     33.15789,  35.05263,    36.94737,  38.84211,  40.73684, &
&     42.63158,  44.52632,    46.42105,  48.31579,  50.21053, &
&     52.10526,  54.00000,    55.89474,  57.78947,  59.68421, &
&     61.57895,  63.47368,    65.36842,  67.26316,  69.15790, &
&     71.05264,  72.94736,    74.84210,  76.73684,  78.63158, &
&     80.52631,  82.42105,    84.31579,  86.21053,  88.10526, &
&     90.00000 /)

YDECMIP%RLONCLI=(/&
&       0.0,   2.5,   5.0,   7.5,  10.0,  12.5,  15.0,  17.5, &
&      20.0,  22.5,  25.0,  27.5,  30.0,  32.5,  35.0,  37.5, &
&      40.0,  42.5,  45.0,  47.5,  50.0,  52.5,  55.0,  57.5, &
&      60.0,  62.5,  65.0,  67.5,  70.0,  72.5,  75.0,  77.5, &
&      80.0,  82.5,  85.0,  87.5,  90.0,  92.5,  95.0,  97.5, &
&     100.0, 102.5, 105.0, 107.5, 110.0, 112.5, 115.0, 117.5, &
&     120.0, 122.5, 125.0, 127.5, 130.0, 132.5, 135.0, 137.5, &
&     140.0, 142.5, 145.0, 147.5, 150.0, 152.5, 155.0, 157.5, &
&     160.0, 162.5, 165.0, 167.5, 170.0, 172.5, 175.0, 177.5, &
&     180.0, 182.5, 185.0, 187.5, 190.0, 192.5, 195.0, 197.5, &
&     200.0, 202.5, 205.0, 207.5, 210.0, 212.5, 215.0, 217.5, &
&     220.0, 222.5, 225.0, 227.5, 230.0, 232.5, 235.0, 237.5, &
&     240.0, 242.5, 245.0, 247.5, 250.0, 252.5, 255.0, 257.5, &
&     260.0, 262.5, 265.0, 267.5, 270.0, 272.5, 275.0, 277.5, &
&     280.0, 282.5, 285.0, 287.5, 290.0, 292.5, 295.0, 297.5, &
&     300.0, 302.5, 305.0, 307.5, 310.0, 312.5, 315.0, 317.5, &
&     320.0, 322.5, 325.0, 327.5, 330.0, 332.5, 335.0, 337.5, &
&     340.0, 342.5, 345.0, 347.5, 350.0, 352.5, 355.0, 357.5 /)

IF (LHOOK) CALL DR_HOOK('SUECOZV:CMIP6_OZONE_COORD',1,ZHOOK_HANDLE)
END SUBROUTINE CMIP6_OZONE_COORD

!SUBROUTINE CMIP7_OZONE_COORD(CFILE, NLON1NC, NLAT1NC, NLV1NC)
! 
! Purpose: 
!   Populate the arrays RLONCLI, RLATCLI, RPROC1 in YDECMIP 
!   which contain the lon, lat and pressure levels for ozone
!
! Joakim Kjellsson, SMHI, 25/09/2025
!

!    CHARACTER(LEN=150), INTENT(IN)   :: CFILE
!    INTEGER(KIND=JPIM), INTENT(IN)   :: NLON1NC, NLAT1NC, NLV1NC 
!    INTEGER(KIND=JPIM)               :: INCUNIT, ILONVARID, ILATVARID, IPREVARID
!    INTEGER(KIND=JPIM)               :: ISTART
!    REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
    
!    IF (LHOOK) CALL DR_HOOK('SUECOZV:CMIP6_OZONE_COORD',0,ZHOOK_HANDLE)

    ! Open file and check for variables
!    CALL CHECK( NF_OPEN(CC,NF_NOWRITE,INCUNIT) ) 

!    CALL CHECK( NF_INQ_VARID(INCUNIT, "lon", ILONVARID) )
!    CALL CHECK( NF_INQ_VARID(INCUNIT, "lat", ILATVARID) )
!    CALL CHECK( NF_INQ_VARID(INCUNIT, "plev", IPREVARID) )  
    
    ! read lon, lat, plev
!    ISTART = 1
!    IF(JPRB==JPRD) THEN ! if double precision 
!        CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, ILONVARID, ISTART, NLON1NC, YDECMIP%RLONCLI) )
!        CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, ILONVARID, ISTART, NLAT1NC, YDECMIP%RLATCLI) )
!        CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, ILONVARID, ISTART, NLV1NC,  YDECMIP%RPROC1) )
!    ELSE ! if single precision
!        CALL CHECK( NF_GET_VARA_REAL(INCUNIT, ILONVARID, ISTART, NLON1NC, YDECMIP%RLONCLI) )
!        CALL CHECK( NF_GET_VARA_REAL(INCUNIT, ILONVARID, ISTART, NLAT1NC, YDECMIP%RLATCLI) )
!        CALL CHECK( NF_GET_VARA_REAL(INCUNIT, ILONVARID, ISTART, NLV1NC,  YDECMIP%RPROC1) )
!    ENDIF

    ! done
!    CALL CHECK( NF_CLOSE(INCUNIT) )

!IF (LHOOK) CALL DR_HOOK('SUECOZV:CMIP6_OZONE_COORD',1,ZHOOK_HANDLE)

!END SUBROUTINE CMIP7_OZONE_COORD

SUBROUTINE READ_NC_FILE_OZONE_CMIP6(CC,NLON1NC,NLAT1NC,NLV1NC,NMONTH1NC,ZOZCL)
  ! Based on EC-Earth code from Michiel van Weele.
  ! Read CMIP6 ozone forcing from NetCDF datafiles (NMONTH1=14 months)
  ! e.g. histo/vmro3_input4MIPs_ozone_CMIP6_UReading-CCMI_1850.nc
  ! store in ZOZCL(JI,JL,JK,JM) with CMIP6 ozone forcing data dimensions (144,96,66,14)
  
  CHARACTER(LEN=150),INTENT(IN)    :: CC
  INTEGER(KIND=JPIM),INTENT(IN)    :: NLON1NC, NLAT1NC ,NLV1NC, NMONTH1NC
  REAL(KIND=JPRB), INTENT(INOUT)   :: ZOZCL(NLON1NC, NLAT1NC ,NLV1NC, 0:NMONTH1NC-1)
  REAL(KIND=JPRB), ALLOCATABLE     :: OZO_CMIP6(:,:,:)
  INTEGER(KIND=JPIM)               :: INCUNIT, OZO_VARID
  INTEGER(KIND=JPIM)               :: IMONTH1, ILEV1
  INTEGER(KIND=JPIM), DIMENSION(4) :: ISTART,ISIZE
  CHARACTER(LEN=*),PARAMETER       :: OZO_NAME='vmro3'
  
  ! OPEN NETCDF FILE
  CALL CHECK( NF_OPEN(CC,NF_NOWRITE,INCUNIT) )
  ALLOCATE(OZO_CMIP6(NLON1NC,NLAT1NC,NLV1NC))
  ! READ 3-D FIELD PER MONTH (14 MONTHS in total)
  CALL CHECK( NF_INQ_VARID(INCUNIT, OZO_NAME, OZO_VARID) )
  ISIZE  = (/NLON1NC,NLAT1NC,NLV1NC,1/)
  DO IMONTH1=1,NMONTH1NC
     ISTART = (/1,1,1,IMONTH1/)
     IF(JPRB==JPRD)THEN
        CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, OZO_VARID, ISTART, ISIZE, OZO_CMIP6) )
     ELSE
        CALL CHECK( NF_GET_VARA_REAL(INCUNIT, OZO_VARID, ISTART, ISIZE, OZO_CMIP6) )
     ENDIF              
     ! Reverse vertical levels, convert from mole/mole to ppm and store in ZOZCL
     DO ILEV1= 1,NLV1NC
        ZOZCL(:,:,ILEV1,IMONTH1-1) = OZO_CMIP6(:,:,NLV1NC-ILEV1+1)*1.E+06_JPRB
     ENDDO
  ENDDO  
  ! CLOSE NETCDF FILE

  CALL CHECK( NF_CLOSE(INCUNIT) )

  DEALLOCATE(OZO_CMIP6)

END SUBROUTINE READ_NC_FILE_OZONE_CMIP6


SUBROUTINE FIND_NC_FILE_OZONE_CMIP7(IYR1, IYEAR1F, IYEAR2F, CC, LFIRSTYR, LLASTYR, CCMIP7_SCEN)

    INTEGER(KIND=JPIM),           INTENT(IN)   :: IYR1              ! year
    CHARACTER(LEN=*), OPTIONAL,   INTENT(IN)   :: CCMIP7_SCEN       ! scenario
    INTEGER(KIND=JPIM),           INTENT(OUT)  :: IYEAR1F, IYEAR2F  ! first and last year of file
    CHARACTER(LEN=80),            INTENT(OUT)  :: CC                ! file name
    LOGICAL,                      INTENT(OUT)  :: LFIRSTYR, LLASTYR ! is it first or last year of file?

    ! Determine which file to read
    ! CMIP7 files (so far) cover periods
    ! 182901-184912
    ! 185001-189912
    ! 190001-194912 
    ! 195001-199912 
    ! 200001-202212 
    ! For piControl (or any fixed year) we just repeat
    SELECT CASE ( IYR1 )
        CASE ( 1829:1849 ) ! 1829 <= IYR <= 1849 
            IYEAR1F = 1829
            IYEAR2F = 1849
        CASE ( 1850:1899 ) ! 1850 <= IYR <= 1899
            IYEAR1F = 1850
            IYEAR2F = 1899
        CASE ( 1900:1949 ) ! 1900 <= IYR <= 1949
            IYEAR1F = 1900
            IYEAR2F = 1949
        CASE ( 1950:1999 ) ! 1950 <= IYR <= 1999
            IYEAR1F = 1950
            IYEAR2F = 1999
        CASE ( 2000:2022 ) ! 2000 <= IYR <= 2022 
            IYEAR1F = 2000
            IYEAR2F = 2022
        CASE ( 2023:2100 ) ! 2023 <= IYR <= 2100
            IYEAR1F = 2023
            IYEAR2F = 2100
            IF ( .NOT. PRESENT(CCMIP7_SCEN) ) THEN
                WRITE(NULOUT,*) "SUECOZV: 2023 <= YEAR <= 2100 but no CMIP7 scenario specified "
                CALL ABOR1("SUECOZV: Can not find ozone file ") 
            END IF 
        CASE DEFAULT       ! else: not included in historical forcing
            ! todo: add scenarios as they become available later
            WRITE(NULOUT,*) "SUECOZV: Can not find ozone data for year ",IYR1
            WRITE(NULOUT,*) "SUECOZV: CMIP7 ozone only works for years 1829-2022 "  
            CALL ABOR1("SUECOZV: No CMIP7 ozone data found ")
    END SELECT
    
    ! Check if IYR1 is the first or last year of file
    ! If so, we will later read the last month of previous year
    ! or first month of next year
    LFIRSTYR = .FALSE.
    LLASTYR  = .FALSE.
    
    IF ( IYR1 == 1849 .OR. &
       & IYR1 == 1899 .OR. &
       & IYR1 == 1949 .OR. &
       & IYR1 == 1999 .OR. &
       & IYR1 == 2022 ) THEN

        LLASTYR = .TRUE.
    
    ELSE IF ( IYR1 == 1850 .OR. &
            & IYR1 == 1900 .OR. & 
            & IYR1 == 1950 .OR. & 
            & IYR1 == 2000 ) THEN

        LFIRSTYR = .TRUE.

    END IF
    
    WRITE(NULOUT,*) "SUECOZV: IYR1, LFIRSTYR, LLASTYR = ",IYR1,LFIRSTYR,LLASTYR 

    ! set file name for historical ozone 
    ! if year1 = 1850 and year2 = 1899 we need to write 185001 and 189912 
    WRITE(CC,'(''ozone/vmro3_input4MIPs_ozone_CMIP_FZJ-CMIP-ozone-1-0_gn_'',I6.6,''-'',I6.6,''.nc'')') &
            & IYEAR1F*100+1, IYEAR2F*100+12
    WRITE(NULOUT,*) "SUECOZV: CC = ",CC

END SUBROUTINE FIND_NC_FILE_OZONE_CMIP7


SUBROUTINE READ_NC_FILE_OZONE_CMIP7(CC,NLON1NC,NLAT1NC,NLV1NC,IYR1NC,IYEAR1NC,IYEAR2NC,NMONTH1NC,ZOZCL)
    !
    ! Based on READ_NC_FILE_OZONE_CMIP6 above 
    !
    ! Purpose: 
    !    Read CMIP7 ozone forcing from NetCDF datafiles. 
    !    50 years for most files 
    !    Dimensions should be nlon=144, nlat=96, nlev=66
    ! 
    ! Method: 
    !    Read ozone data for the entire year at the beginning of each year.
    !    Include Dec from year-1 and Jan from year+1 (total 14 months) to allow
    !    time interpolations to cover < 15 Jan and > 15 Dec.  
    !    Store in ZOZCL(JI,JL,JK,JM) with CMIP7 ozone forcing data dimensions (144,96,66,14)
    !
    ! Note: 
    !    CMIP7 ozone has the same name and dimensions as CMIP6 so we could have re-used CMIP6 routines
    !    But it felt better to re-do it in case CMIP7 later upgrades to higher spatial or temporal resolution
    !
    CHARACTER(LEN=*),INTENT(IN)      :: CC ! file to read
    INTEGER(KIND=JPIM),INTENT(IN)    :: NLON1NC, NLAT1NC ,NLV1NC ! nlon, nlat, nlev
    INTEGER(KIND=JPIM),INTENT(IN)    :: IYR1NC, NMONTH1NC ! year to read and number of months (usually 14 months)
    INTEGER(KIND=JPIM),INTENT(IN)    :: IYEAR1NC, IYEAR2NC ! first and last year of file
    REAL(KIND=JPRB), INTENT(INOUT)   :: ZOZCL(NLON1NC, NLAT1NC ,NLV1NC, 0:NMONTH1NC+1) ! ozone field for 2 extra months 
    REAL(KIND=JPRD), ALLOCATABLE     :: ZLON(:), ZLAT(:), ZLV(:) 
    REAL(KIND=JPRB), ALLOCATABLE     :: OZO_CMIP7(:,:,:) ! 3d field of ozone to read from file 
    INTEGER(KIND=JPIM)               :: INCUNIT, ILONVARID, ILATVARID, IPREVARID, IOZOVARID ! ids for netcdf reading
    INTEGER(KIND=JPIM)               :: IMONTH1, ILEV1 ! month and level indices
    INTEGER(KIND=JPIM), DIMENSION(4) :: ISTART,ISIZE ! start index and size of netcdf arrays
    CHARACTER(LEN=*),PARAMETER       :: COZONAME='vmro3' ! name of ozone variable in netcdf files

    ! open netcdf file and read coordinates
    CALL CHECK( NF_OPEN(CC,NF_NOWRITE,INCUNIT) )
    
    CALL CHECK( NF_INQ_VARID(INCUNIT, "lon", ILONVARID) )
    CALL CHECK( NF_INQ_VARID(INCUNIT, "lat", ILATVARID) )
    CALL CHECK( NF_INQ_VARID(INCUNIT, "plev", IPREVARID) )
    
    ALLOCATE( ZLON(NLON1NC), ZLAT(NLAT1NC), ZLV(NLV1NC) )

    ! read lon, lat, plev
    ! These are DOUBLE in the data, so always read double precision
    CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, ILONVARID, (/1/), (/NLON1NC/), ZLON) )
    CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, ILATVARID, (/1/), (/NLAT1NC/), ZLAT) )
    CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, IPREVARID, (/1/), (/NLV1NC/),  ZLV) )
    
    !WRITE(NULOUT, *) "SUECOZV: ZLON = ",ZLON
    !WRITE(NULOUT, *) "SUECOZV: ZLAT = ",ZLAT
    !WRITE(NULOUT, *) "SUECOZV: ZLV  = ",ZLV

    !
    ! read ozone for this year as well as Dec of year-1 and Jan of year+1
    !

    ALLOCATE(OZO_CMIP7(NLON1NC,NLAT1NC,NLV1NC)) ! 3d field to read

    ! indices to read from file 
    ISIZE  = (/ NLON1NC, NLAT1NC, NLV1NC, 1 /)
  
    ! READ 3-D FIELD 
    CALL CHECK( NF_INQ_VARID(INCUNIT, COZONAME, IOZOVARID) )
    
    DO IMONTH1 = 1,NMONTH1NC ! loop from 1 to 14
        
        ! If IYEAR1NC is 1850 and IYR is 1860, 
        ! then (IYR1NC - IYEAR1NC)*12 + 1 = 10*12 = 120 
        ! which would be Dec 1859. Then we read 1 month at a time
        ISTART = (/       1,       1,      1, (IYR1NC - IYEAR1NC)*12+IMONTH1-1 /)
        
        ! If it is the first year and first month of file, the previous index does not exist
        ! Read first index anyway
        IF (IYR1NC == IYEAR1NC .AND. IMONTH1 == 1) THEN
            ISTART(4) = 1
        ! If it is the last year and last month of file, the next index does not exist
        ! Read the last again
        ELSE IF (IYR1NC == IYEAR2NC .AND. IMONTH1 == NMONTH1NC) THEN
            ISTART(4) = (IYR1NC - IYEAR1NC)*12+IMONTH1-2
        END IF
            
        WRITE(NULOUT,*) "SUECOZV: ISTART ", ISTART
        WRITE(NULOUT,*) "SUECOZV: ISIZE ", ISIZE 
        
        IF(JPRB==JPRD)THEN ! if double precision
            ! Joakim: I am skeptical about this read. 
            ! The vmro3 field in the netCDF file is single precision but here we read into 
            ! double precision. I suppose it just pads with garbage decimals at the end. 
            ! A better solution would be to always read SP and then cast to DP if needed. 
            WRITE(NULOUT,*) "SUECOZV: Read 3D ozone in DP "
            CALL CHECK( NF_GET_VARA_DOUBLE(INCUNIT, IOZOVARID, ISTART, ISIZE, OZO_CMIP7(:,:,:)) )
        ELSE ! single precision
            WRITE(NULOUT,*) "SUECOZV: Read 3D ozone in SP "
            CALL CHECK( NF_GET_VARA_REAL  (INCUNIT, IOZOVARID, ISTART, ISIZE, OZO_CMIP7(:,:,:)) )
        ENDIF
        
        ! Reverse vertical levels, convert from mole/mole to ppm and store in ZOZCL
        ! Put first month on index 0 of ZOZCL so that second month (Jan of YEAR) is index 1
        DO ILEV1= 1,NLV1NC
            ZOZCL(:,:,ILEV1,IMONTH1-1) = OZO_CMIP7(:,:,NLV1NC-ILEV1+1)*1.E+06_JPRB
        ENDDO
          
    ENDDO
    
    ! Reverse order of plev as well 
    DO ILEV1= 1,NLV1NC
        YDECMIP%RPROC1(ILEV1)      = ZLV(NLV1NC-ILEV1+1)
    ENDDO
    
    ! plev is size NLV1NC from 0.01 Pa to 100000 Pa
    ! RPROC1 is 0:NLV1NC+1 and should start at 0 Pa 
    ! and end at 110000 Pa
    YDECMIP%RPROC1(0)        = 0.0_JPRB
    YDECMIP%RPROC1(NLV1NC+1) = 110000.0_JPRB 
    
    ! lon and lat
    ! Joakim: ZLON, ZLAT and ZLV are double precision in the netCDF file
    ! If we run single precision, I think these will be cast to single precision
    ! here, but not sure if this is the best way to do it. 
    YDECMIP%RLONCLI(:) = ZLON(:)
    YDECMIP%RLATCLI(:) = ZLAT(:) 
    WRITE(NULOUT, *) "SUECOZV: RLONCLI = ",YDECMIP%RLONCLI(:)
    WRITE(NULOUT, *) "SUECOZV: RLATCLI = ",YDECMIP%RLATCLI(:)
    WRITE(NULOUT, *) "SUECOZV: RPROC1  = ",YDECMIP%RPROC1(:) 

    ! CLOSE NETCDF FILE
    CALL CHECK( NF_CLOSE(INCUNIT) )

    DEALLOCATE(ZLON, ZLAT, ZLV, OZO_CMIP7)

END SUBROUTINE READ_NC_FILE_OZONE_CMIP7

SUBROUTINE CHECK(STATUS)
  INTEGER(KIND=JPIM), INTENT (IN) :: STATUS
  IF(STATUS /= NF_NOERR) THEN
     CALL ABOR1('SUECOZV:  '//TRIM(NF_STRERROR(STATUS)))
  ENDIF
END SUBROUTINE CHECK

END SUBROUTINE SUECOZV


