! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE WROUTSPGB(YDMODEL,YDGEOMETRY,YDDIM,PSPFLD,KSHMAX,KGRIBIO,CDLEVTYPE,KRESOL,CDFNSH,YDFPFIELDS,PTSTEP,YDQTYPE)

!**** *WROUTSPGB* _ GRIB codes and writes out spectral fields

!     Purpose.
!     --------
!     Write out spectral fields in GRIB

!**   Interface.
!     ----------
!        *CALL* *WROUTSPGB*(...)

!        Explicit arguments :    
!        --------------------

!        Implicit arguments :      The state variables of the model
!        --------------------

!     Method.
!     -------
!        See documentation
!        - spectral part would only work if input-resol = output resol
!          for vertical and horizontal resolution

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*
!        Original : 01-12-14 Adapted from WR****

!     Modifications.
!     --------------
!       P.Towers : 02-11-13 Fixes for fewer writers than gatherers
!       P.Towers : 03-04-23 Added ISETFIELDCOUNTFDB logic
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!        G. Radnoti:   25-Aug-2004 treatment of tstep<0 for LTWINC
!        M.Hamrud      01-Dec-2005 Generalized IO scheme
!       R. El Khatib : 23-Oct-2008 use suofname
!       R. El Khatib : 20-Aug-2012 optional argument KRESOL and consequences
!       G. Carver :    22-May-2013 fixed write mode for non-FDB I/O
!       K. Yessad (July 2014): Move some variables.
!       R. El Khatib : 17-Mar-2016 CFPATH passed to suofname
!     ------------------------------------------------------------------

USE TYPE_MODEL, ONLY: MODEL
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB, JPIB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : NCONF
USE YOMCT3   , ONLY : NSTEP
USE YOMLUN   , ONLY : NULOUT
USE YOMVAR   , ONLY : MCGLVEC, MBGVEC,&
 & LTRREF, LTWANA, LTWGRA, LTWINC, LTWBGV, LTWCGL  
USE YOMOPH0  , ONLY : CFNAN, CFNBGV, CFNCGL, CFNGR, CFNINSH, CFNRF, LINC
USE YOMMP0   , ONLY : NOUTTYPE
USE YOMRIP0  , ONLY : NINDAT, NSSSSS
USE IOSTREAM_MIX , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 & CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST, Y_IOSTREAM_FDB,&
 & CLOSE_IOREQUEST
USE ALGORITHM_STATE_MOD, ONLY : GET_NSIM4D
USE TYPE_FPFIELDS, ONLY : TFPFIELDS
USE TYPE_FPRQDYNS, ONLY : TYPE_FPRQDYN

#ifdef WITH_XIOS
USE YOMCT0   , ONLY : LXIOS
USE CXIOS    , ONLY : XIOS_PUT
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)       , INTENT(INOUT) :: YDMODEL
TYPE (GEOMETRY)   , INTENT(IN) :: YDGEOMETRY
TYPE(TDIM)        , INTENT(IN) :: YDDIM
INTEGER(KIND=JPIM),PARAMETER :: ISHOUR=3600
INTEGER(KIND=JPIM),PARAMETER :: IHDAY=24

INTEGER(KIND=JPIM), INTENT(IN) :: KSHMAX 
REAL(KIND=JPRB)   , INTENT(IN) :: PSPFLD(:,:) 
INTEGER(KIND=JPIM), INTENT(IN) :: KGRIBIO(4,KSHMAX) 
CHARACTER(LEN=1)  , INTENT(IN) :: CDLEVTYPE 
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
CHARACTER(LEN=210), INTENT(IN) :: CDFNSH
REAL(KIND=JPRB)   , INTENT(IN), OPTIONAL :: PTSTEP
TYPE(TYPE_FPRQDYN), INTENT(IN), OPTIONAL :: YDQTYPE
INTEGER(KIND=JPIM)             :: ILEN

!     ------------------------------------------------------------------

INTEGER(KIND=JPIB) :: IINC,IMTS
INTEGER(KIND=JPIM) :: IHOUR
INTEGER(KIND=JPIB) :: ISEC
INTEGER(KIND=JPIM) :: IH0,IJ0,IM0,IA0,IDD,ISS,IHR,IMIN,ISC,IJOUR,IMOIS,IAN,ILMOIS(12)
INTEGER(KIND=JPIM) :: IYM
CHARACTER :: CLFNSH*210,CLFNRF*30,CLFANA*30,CLR1*5,CLR2*4,CLMODE*1
CHARACTER :: CLFINC*30
CHARACTER :: CLFGRA*30
CHARACTER :: CLFCGL*30
CHARACTER :: CLFBGV*30
CHARACTER :: CLEVT*2

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE(TFPFIELDS)     :: YDFPFIELDS
TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "suofname.intfb.h"
#include "fcttim.func.h"
#include "updcalsec.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WROUTSPGB',0,ZHOOK_HANDLE)
ASSOCIATE(YRRIP=>YDMODEL%YRML_GCONF%YRRIP)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,NSPEC2=>YDGEOMETRY%YRDIM%NSPEC2, &
 & NSTOP=>YRRIP%NSTOP, TSTEP=>YRRIP%TSTEP, NYMLPPSP=>YRRIP%NYMLPPSP)
!     ------------------------------------------------------------------

IA0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
IJ0=NDD(NINDAT)
IH0=NSSSSS/ISHOUR
IF (NSTEP <= 0) NYMLPPSP=0
ISEC=NINT(NSTEP*TSTEP,JPIB)
IHOUR=INT(ISEC/ISHOUR,JPIM)
IDD=IHOUR/IHDAY
ISS=MODULO(IHOUR,IHDAY)*ISHOUR
CALL UPDCALSEC(IH0,IJ0,IM0,IA0,IDD,ISS,IHR,IMIN,ISC,IJOUR,IMOIS,IAN,ILMOIS,NULOUT)
IYM=IAN*100+IMOIS

IF(CDLEVTYPE == 'm') THEN
  CLEVT='ML'
ELSEIF(CDLEVTYPE == 'p') THEN
  CLEVT='PL'
ELSEIF(CDLEVTYPE == 'v') THEN
  CLEVT='PV'
ELSEIF(CDLEVTYPE == 't') THEN
  CLEVT='TH'
ELSE
  CALL ABOR1('WROUTSPGB:UNKNOWN LEVEL TYPE')
ENDIF

CLMODE = 'a'

IF(NOUTTYPE /= 2.AND. KSHMAX > 0 )THEN

!      Open output file

  IF (LTRREF) THEN
    
    CLR1 = CFNRF(1:5)
    CLR2 = CFNRF(6:9)
    IF (NSTEP < 0) THEN
      WRITE(CLFNRF,'(A5,A4,I7.6)') CLR1,CLR2,NSTEP
    ELSE
      IF (LINC) THEN
        IINC = NINT(REAL(NSTEP,JPRB)*TSTEP/3600._JPRB,JPIB)
      ELSE
        IINC = NSTEP
      ENDIF
      WRITE(CLFNRF,'(A5,A4,''+'',I6.6)') CLR1,CLR2,IINC
    ENDIF
    CLFNSH=CLFNRF

  ELSEIF(LTWANA) THEN
      
    CLFANA = CFNAN
    IF (NSTEP < 0) THEN
      CLFANA(10:10) = '-'
      IINC = ABS(NSTEP)
    ELSEIF (NCONF/100 == 1) THEN
      IMTS = NINT(REAL(NSTEP,JPRB)*TSTEP/60._JPRB,JPIB)
      IINC = (IMTS/60)*100+MOD(IMTS,60)
    ELSE
      IINC = NSTEP
    ENDIF
    WRITE(CLFANA(7:9),'(I3.3)') GET_NSIM4D()
    WRITE(CLFANA(11:),'(I6.6)') IINC
    CLFNSH=CLFANA
      
  ELSEIF(LTWINC) THEN
    
    CLFINC = CFNINSH
! Output files labelled with time step in configuration 501 (TL evolution)
    IF (NCONF == 501) THEN
      IINC = NSTEP
    ELSE
      IMTS = NINT(REAL(NSTEP,JPRB)*TSTEP/60._JPRB,JPIB)
      IINC = (IMTS/60)*100+MOD(IMTS,60)
    ENDIF
    WRITE(CLFINC(9:11),'(I3.3)') GET_NSIM4D()
    IF (IINC >= 0) THEN
      WRITE(CLFINC(13:),'(I6.6)') IINC
    ELSE
      WRITE(CLFINC(12:12),'(A1)') '-'
      WRITE(CLFINC(13:),'(I6.6)') -IINC
    ENDIF
    CLFNSH=CLFINC
    
  ELSEIF(LTWGRA) THEN
      
    CLFGRA = CFNGR
    IF (NSTEP < 0) THEN
      CLFGRA(10:10) = '-'
      IINC = ABS(NSTEP)
    ELSEIF (NCONF/100 == 1) THEN
      IMTS = NINT(REAL(NSTEP,JPRB)*TSTEP/60._JPRB,JPIB)
      IINC = (IMTS/60)*100+MOD(IMTS,60)
    ELSEIF (IABS(NCONF)/100 == 8) THEN
      CLFGRA(10:10) = '-'
      IINC = NSTOP-NSTEP
    ELSE
      IINC = NSTEP
    ENDIF
    WRITE(CLFGRA(7:9),'(I3.3)') GET_NSIM4D()
    WRITE(CLFGRA(11:),'(I6.6)') IINC
    CLFNSH=CLFGRA
    
  ELSEIF(LTWBGV) THEN
    
    CLFBGV = CFNBGV
    WRITE(CLFBGV(7:),'(I3.3)') MBGVEC
    CLFNSH=CLFBGV
      
  ELSEIF(LTWCGL) THEN
    
    CLFCGL = CFNCGL
    WRITE(CLFCGL(7:),'(I3.3)') MCGLVEC
    CLFNSH=CLFCGL
    
  ELSE

    CALL SUOFNAME(KSTEP=NSTEP,KDIGITS=6,CDLABEL='ICMSH',CDOFNAME=CLFNSH)
    IF(NSTEP>0) THEN
      ILEN=LEN(TRIM(CLFNSH))-10
      IF ((NYMLPPSP == 0).OR.(IYM == NYMLPPSP).OR. &
    ((IYM /= NYMLPPSP).AND.((IJOUR /= 1).OR.(IHR /= 0)))) THEN
        WRITE(CLFNSH(ILEN:),'("+",I6)') IYM
      ELSE
        WRITE(CLFNSH(ILEN:),'("+",I6)') NYMLPPSP
      ENDIF
    ENDIF
  ENDIF

ENDIF

CLMODE='a'

IF ((NSTEP <= 0).OR.((NYMLPPSP /= 0).AND. &
 & (IYM /= NYMLPPSP).AND.(IJOUR == 1).AND.(IHR == 0))) THEN
  NYMLPPSP = NYMLPPSP
ELSE
  NYMLPPSP = IYM
ENDIF

CALL SETUP_IOREQUEST(YL_IOREQUEST,'SPECTRAL_FIELDS',LDGRIB=.TRUE.,&
 & KGRIB2D=KGRIBIO(1,:),KLEVS2D=KGRIBIO(2,:),KBSET2D=KGRIBIO(3,:),&
 & CDLEVTYPE=CLEVT,KRESOL=KRESOL,PTSTEP=PTSTEP)

IF(NOUTTYPE == 2) THEN
  CALL IO_PUT(Y_IOSTREAM_FDB,YL_IOREQUEST,PR2=PSPFLD)
ELSE
#ifdef WITH_XIOS
  IF (LXIOS) THEN
    ! XIOS_FPOS extra logging
    IF(PRESENT(YDQTYPE)) THEN
      WRITE(NULOUT, '(''XIOSFPOS: ENTERING XIOS FIELD SENDING PROCEDURES FOR SP FIELDS WTIH YDQTYPE'')')
      WRITE(UNIT=NULOUT,FMT='('' NFPISCAG in WROUTSPGB= '',I4)') YDQTYPE%NFPISCAG
      CALL XIOS_PUT(YL_IOREQUEST,YDFPFIELDS,YDQTYPE=YDQTYPE,CMODE=CLMODE,PR2=PSPFLD)
    ELSE
      WRITE(NULOUT, '(''XIOSFPOS: SKIPPTING XIOS FIELD SENDING PROCEDURES FOR SP FIELDS WITHOUT YDQTYPE'')')
      !CALL XIOS_PUT(YL_IOREQUEST,YDFPFIELDS,CMODE=CLMODE,PR2=PSPFLD)
    ENDIF
  ELSE
    CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFNSH),CDMODE='a',KIOMASTER=1)
    CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR2=PSPFLD)
    CALL CLOSE_IOSTREAM(YL_IOSTREAM)
  ENDIF
#else
  CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFNSH),CDMODE='a',KIOMASTER=1)
  CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR2=PSPFLD)
  CALL CLOSE_IOSTREAM(YL_IOSTREAM)
#endif
ENDIF
 
CALL CLOSE_IOREQUEST(YL_IOREQUEST)

END ASSOCIATE
END ASSOCIATE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WROUTSPGB',1,ZHOOK_HANDLE)
END SUBROUTINE WROUTSPGB
