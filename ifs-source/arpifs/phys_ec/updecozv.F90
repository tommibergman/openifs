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
SUBROUTINE UPDECOZV(YDERDI,YDECMIP,KINDAT,KMINUT)

!**** *UPDECOZV* - GETS VARIATIONAL DISTRIBUTION OF OZONE
!                 (FOR CMIP RUNS)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *UPDECOZV* FROM *UPDTIM*

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

!        Separated out from SUECOZV : SUEC does the initial setup, this is called at appropriate time steps from updtim

!     AUTHOR.
!     -------
!     O. Marsden  - January 2018
 
!     MODIFICATIONS.
!     --------------
!     J. Kjellsson 26/09/2025 - Support for CMIP7 forcings 
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD, JPIB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RPI, RDAY
USE YOERDI   , ONLY : TERDI
USE YOECMIP  , ONLY : TECMIP,NLON1_CMIP5, NLAT1_CMIP5, NLV1_CMIP5, NMONTH1, &
 &                    NLON1_CMIP6, NLAT1_CMIP6, NLV1_CMIP6, &
 &                    NLON1_CMIP7, NLAT1_CMIP7, NLV1_CMIP7 
USE ECE_CMIP,  ONLY : NCMIPFIXYR

IMPLICIT NONE

!*       0.1   ARGUMENTS.
!              ----------
TYPE(TERDI)       ,INTENT(IN)    :: YDERDI
TYPE(TECMIP)      ,INTENT(INOUT) :: YDECMIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KINDAT 
INTEGER(KIND=JPIB),INTENT(IN)    :: KMINUT 
!     -----------------------------------------------------------------

!*       0.2   LOCAL ARRAYS.
!              -------------
REAL(KIND=JPRB),ALLOCATABLE, SAVE :: ZYTIME(:), ZMDAY(:)


INTEGER(KIND=JPIM) :: IDY, IM1, IM2, IMN, JK, JL, JI, JM, NINDAT
INTEGER(KIND=JPIM) :: IH0,IJ0,IM0,IA0,IDD,ISS, IHR,IMIN,ISC,IYR,ILMOIS(12)

INTEGER(KIND=JPIM) :: I, IUNIT, IDIR, IFIL
LOGICAL            :: LLIS_OPEN
CHARACTER(LEN=132) :: CLSKIP_LINE

CHARACTER (LEN = 300) ::  CLFN

REAL(KIND=JPRB) :: ZTIMI, ZXTIME
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: NLON1, NLAT1 ,NLV1

#include "abor1.intfb.h"
#include "updcalsec.intfb.h"
#include "suecozv.intfb.h"

#include "fcttim.func.h"      
#include "netcdf.inc"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDECOZV',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------


!! Ugly trick to keep this stuff local, but to avoid doing it at each time step. 
!! Should be improvable once IFS time-handling is cleaned up. O Marsden
IF ( .NOT. ALLOCATED(ZMDAY) ) THEN

  ALLOCATE(ZMDAY(12))
  ALLOCATE(ZYTIME(12))
  
  ZMDAY = (/&
   & 31._JPRB,    59.25_JPRB,  90.25_JPRB, 120.25_JPRB, 151.25_JPRB, 181.25_JPRB&
   & ,  212.25_JPRB, 243.25_JPRB, 273.25_JPRB, 304.25_JPRB, 334.25_JPRB, 365.25_JPRB&
   & /)  
  ZYTIME= (/&
   & 22320._JPRB,  64980._JPRB, 107640._JPRB, 151560._JPRB, 195480._JPRB, 239400._JPRB&
   & , 283320._JPRB, 327960._JPRB, 371880._JPRB, 415800._JPRB, 459720._JPRB, 503640._JPRB&
   & /)  

ENDIF


!*         1.     TIME INDEX WITHIN OZONE CLIMATOLOGY
!                 -----------------------------------

!ECEARTH: Time interpolation of climatology for any type of calendar. (See also suecaec.F90!)

IA0=NCCAA(KINDAT)
IM0=NMM(KINDAT)
IJ0=NDD(KINDAT)
IH0=0
IDD=KMINUT/1440_JPIB ! The number of days since KINDAT
ISS=MODULO(KMINUT,1440_JPIB)*60_JPIB ! The number of seconds since the start of the current day
CALL UPDCALSEC(IH0,IJ0,IM0,IA0,IDD,ISS, IHR,IMIN,ISC,IDY,IMN,IYR,ILMOIS,-1)

ZXTIME=IDY-1+(60*(60*IHR+IMIN)+ISC)/RDAY - 0.5*ILMOIS(IMN) ! Number of days relative to center of the current month
IF (ZXTIME < 0.0_JPRB) THEN
  IM1=MODULO(IMN-2,12)+1
  IM2=IMN
  ZXTIME=ZXTIME + 0.5*(ILMOIS(IM1) + ILMOIS(IM2)) ! Adjust relative time to the center of the previous month
ELSE
  IM1=IMN
  IM2=MODULO(IMN,12)+1
ENDIF
ZTIMI=ZXTIME/(0.5*(ILMOIS(IM1) + ILMOIS(IM2))) ! Compute interpolation weight between the centers of the two months

IF (IMN == 1.AND.IM1 == 12) THEN
  IM1=0
ELSEIF (IMN == 12.AND.IM2 == 1) THEN
  IM2=13
ENDIF

!*         1.     TIME INDEX WITHIN OZONE CLIMATOLOGY
!                 -----------------------

! Apply fixed year override BEFORE range checks to avoid spurious warnings
IF (NCMIPFIXYR>0) IYR=NCMIPFIXYR

! SET TIME INTERVAL
IF (YDECMIP%NO3CMIP == 7) THEN ! CMIP7

    ! Note: These are the same as for CMIP6
    ! at least in FZJ-CMIP-ozone-1-0 dataset 
    NLON1 = NLON1_CMIP7
    NLAT1 = NLAT1_CMIP7
    NLV1  = NLV1_CMIP7  

    IF (IYR < 1850) THEN
        WRITE(NULOUT,*) "UPDECOZV: IYR < 1850 not supported in CMIP7 ozone v2.0. Will repeat 1850 for now..."
        IYR=MAX(IYR,1850)
        WRITE(NULOUT,*) "UPDECOZV: IYR = ",IYR 
    ENDIF
    IF (IYR > 2022) THEN
        WRITE(NULOUT,*) "UPDECOZV: IYR > 2022 not supported in CMIP7 (yet). Will repeat 2022 for now..."
        IYR=MIN(IYR,2022)
        WRITE(NULOUT,*) "UPDECOZV: IYR = ",IYR 
    ENDIF

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
  CALL ABOR1('UPDECOZV: Value of NO3CMIP not supported')
ENDIF

!! rerun setup of ozv *if* current year is not that stored in YDECMIP
! Joakim: Also check if we want to read NCMIPFIXYR, but it is not what we have read before
IF (YDECMIP%NCURRYR /= IYR .OR. (NCMIPFIXYR>0 .AND. YDECMIP%NCURRYR /= NCMIPFIXYR)) THEN
  ! Joakim: SUECOZV is called with KINDAT which is the start date
  ! so if model starts in 1990, we will always read 1990 even 
  ! if current year is 2050. 
  ! So we need to make a new date...
  !CALL SUECOZV(YDECMIP,KINDAT) ! original code
  NINDAT=IYR*10000+IMN*100+IDY ! 8-digit YYYYMMDD
  CALL SUECOZV(YDECMIP,NINDAT)
ENDIF


! check that SUECOZV has been called during the initialization phase
IF (.NOT. ( ALLOCATED(YDECMIP%RSINC1)  .AND. ALLOCATED(YDECMIP%ROZT1)    .AND. ALLOCATED(YDECMIP%RPROC1) .AND. &
 &          ALLOCATED(YDECMIP%RLATCLI) .AND. ALLOCATED(YDECMIP%RLONCLI)) .AND. ALLOCATED(YDECMIP%ZOZCL) ) THEN
  CALL ABOR1('Calling UPDEXOZV without first having called SUECOZV ')
ENDIF
  

!*         2.0    TIME INTERPOLATED FIELD
!                 -----------------------

!*( Field is also transformed in kg/kg! )
DO JK=1,NLV1
  DO JL=1,NLAT1
    DO JI=1,NLON1
      YDECMIP%ROZT1(JI,JL,JK)=YDERDI%RO3 * (YDECMIP%ZOZCL(JI,JL,JK,IM2)&
        & +ZTIMI*(YDECMIP%ZOZCL(JI,JL,JK,IM1)-YDECMIP%ZOZCL(JI,JL,JK,IM2)))
    ENDDO
  ENDDO
ENDDO

DO JL=1,NLAT1
  DO JI=1,NLON1
    YDECMIP%ROZT1(JI,JL, 0)=0.0_JPRB
    YDECMIP%ROZT1(JI,JL,NLV1+1)=YDECMIP%ROZT1(JI,JL,NLV1)
  ENDDO
ENDDO


IF (LHOOK) CALL DR_HOOK('UPDECOZV',1,ZHOOK_HANDLE)

!-----------------------------------------------------------------------------

END SUBROUTINE UPDECOZV


