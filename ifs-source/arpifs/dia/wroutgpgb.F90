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

SUBROUTINE WROUTGPGB(YDMODEL,YDGEOMETRY,PGGFLD,KGPTOT,KGGMAX,KGRIBIO,CDLEVTYPE,CDMODE,KRESOL,CDFN,PTSTEP)

!**** *WROUTGPGB* - GRIB codes and writes out grid-point fields

!     Purpose.
!     --------
!     Write out gridpoint fields in GRIB

!**   Interface.
!     ----------
!        *CALL* *WROUTGPGB*(...)

!        Explicit arguments : PGGFLD  - fields to be written out   
!        -------------------- KGPTOT  - number of gridpoints on your proc.
!                             KGGMAX  - number of fields to write out
!                             KGRIBIO - help array for GRIB coding
!                             CDLEVTYPE - type of level for fields
!                             CDMODE - mode to open files 
!                             KRESOL - resolution tag
!                             CDFN - output file name
!                             PTSTEP - model time step

!        Implicit arguments :      
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    - PRGRIBENC - modify GRIB headers, output resolution
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*
!       Original : 01-12-14 Adapted from WR****

!     Modifications.
!     --------------
!       J.Hague(IBM)/M.Hamrud 15-03-02  - Improvements in message passing
!       P.Towers 13-11-02 - Fixes for fewer writers than gatherers
!       John Hague : 13-01-03 WAIT if NPROC=1
!       P.Towers 23-04-03 - Added ISETFIELDCOUNTFDB logic
!       M.Hamrud      01-Oct-2003 CY28 Cleaning
!       M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!       M.Hamrud      01-Dec-2005 Generalized IO scheme
!       R. El Khatib : 01-Mar-2012 LFPOS => LECFPOS
!       T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!       K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      Ph. Lopez 24-Feb-2017 Added option to write out trajectory GP fields 
!                            in TL evolution experiments (conf 501 with LDTLEVOL=T). 
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPIB, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LECFPOS, NFPOS, LXIOS
USE YOMFPC   , ONLY : TNAMFPOBJ
USE YOMCT3   , ONLY : NSTEP
USE YOMGRIB  , ONLY : NSTEPLPP  
USE YOMOPH0  , ONLY : CFNBGV, CFNGG, CFNINGG, LINC
USE YOMMP0   , ONLY : NOUTTYPE, NPROC
USE YOMRIP0  , ONLY : NINDAT, NSSSSS
USE YOMVAR   , ONLY : LTWINC, MBGVEC, LTWBGV
USE ALGORITHM_STATE_MOD  , ONLY : GET_NSIM4D
USE IOSTREAM_MIX , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, IO_PUT,&
 & CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST, Y_IOSTREAM_FDB,&
 & CLOSE_IOREQUEST
USE MPL_MODULE , ONLY : MPL_ALLREDUCE
USE TYPE_FPFIELDS, ONLY : TFPFIELDS
USE YOMLUN , ONLY : NULOUT

#ifdef WITH_XIOS
USE YOMCT0   , ONLY : LXIOS
USE CXIOS    , ONLY : XIOS_PUT
#endif

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),PARAMETER :: ISHOUR=3600
INTEGER(KIND=JPIM),PARAMETER :: IHDAY=24

TYPE(MODEL)       , INTENT(IN) :: YDMODEL
TYPE(GEOMETRY)    , INTENT(IN) :: YDGEOMETRY
TYPE (TFPFIELDS)  , INTENT(IN) :: YDFPFIELDS
INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT 
INTEGER(KIND=JPIM), INTENT(IN) :: KGGMAX 
REAL(KIND=JPRB)   , INTENT(IN) :: PGGFLD(KGPTOT,KGGMAX) 
INTEGER(KIND=JPIM), INTENT(IN) :: KGRIBIO(2,KGGMAX) 
CHARACTER(LEN=1)  , INTENT(IN) :: CDLEVTYPE 
CHARACTER(LEN=1)  , INTENT(IN) :: CDMODE
INTEGER(KIND=JPIM), INTENT(IN) :: KRESOL
CHARACTER(LEN=*)  , INTENT(IN) :: CDFN
REAL(KIND=JPRB)   , INTENT(IN), OPTIONAL :: PTSTEP

!     ------------------------------------------------------------------

INTEGER(KIND=JPIB) :: IINC, IMTS
INTEGER(KIND=JPIM) :: IHOUR
INTEGER(KIND=JPIB) :: ISEC
INTEGER(KIND=JPIM) :: IH0,IJ0,IM0,IA0,IDD,ISS,IHR,IMIN,ISC,IJOUR,IMOIS,IAN,ILMOIS(12)
INTEGER(KIND=JPIM) :: IYM
INTEGER(KIND=JPIM) :: NYMLPPGP
CHARACTER :: CLEVT*3
CHARACTER :: CLR1*5,CLR2*4,CLR1INC*4
CHARACTER :: CLR3*5,CLR4*5,CLR5*5  !KPP CLR5
CHARACTER :: CLFNGG*30
CHARACTER :: CLFNUA*30
CHARACTER :: CLFNOC*30  !KPP
CHARACTER :: CLFN*30
CHARACTER :: CLFBGV*30
CHARACTER :: CLMODE*1

TYPE(TNAMFPOBJ) :: YLNAMFPOBJ
TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST
REAL(KIND=JPRB) :: ZGGFLD(KGPTOT,KGGMAX,1)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "fcttim.func.h"
#include "updcalsec.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WROUTGPGB',0,ZHOOK_HANDLE)

ASSOCIATE(YRRIP=>YDMODEL%YRML_GCONF%YRRIP)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,TSTEP=>YRRIP%TSTEP, NYMLPPGP_ASSOC=>YRRIP%NYMLPPGP, &
 & CFPFMT=>YLNAMFPOBJ%CFPFMT)
!     ------------------------------------------------------------------

NYMLPPGP = NYMLPPGP_ASSOC

!  PREPARATIONS
IA0=NCCAA(NINDAT)
IM0=NMM(NINDAT)
IJ0=NDD(NINDAT)
IH0=NSSSSS/ISHOUR
IF (NSTEP <= 0) NYMLPPGP=0
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
ELSEIF(CDLEVTYPE == 's') THEN
  CLEVT='SFC'
ELSEIF(CDLEVTYPE == 'o') THEN  !KPP
  CLEVT='OML'                  !KPP
ELSE
  CALL ABOR1('WROUTGPGB:UNKNOWN LEVEL TYPE')
ENDIF


! CREATE FILENAMES AND/OR PREPARE FOR WRITING

IF(NOUTTYPE /= 2) THEN  ! Not writing to FDB

  CLFNGG = ' '
  CLR1 = CFNGG(1:5)
  IF (LTWINC) CLR1INC = CFNINGG(1:4)
  CLR2 = CFNGG(6:9)
  CLR3 = 'ICMFP'
  CLR4 = 'ICMUA'
  CLR5 = 'ICMOC'   !KPP
  IF (LINC) THEN
    IINC = NINT(REAL(NSTEP,JPRB)*TSTEP/3600._JPRB,JPIB)
  ELSE
    IINC = NSTEP
  ENDIF

  IF (NSTEP < 0) THEN
    IF (LECFPOS) THEN
      WRITE(CLFNGG,'(A5,A4,I5.4)') CLR1,CLR2,NSTEP
      WRITE(CLFNUA,'(A5,A4,I5.4)') CLR4,CLR2,NSTEP
    ELSE
      WRITE(CLFNGG,'(A5,A4,I5.4)') CLR1,CLR2,NSTEP
    ENDIF
    WRITE(CLFNOC,'(A5,A4,I5.4)') CLR5,CLR2,NSTEP      !KPP  
  ELSE
    IF(LTWINC) THEN
      IMTS = NINT(REAL(NSTEP,JPRB)*TSTEP/60._JPRB,JPIB)
      IINC = (IMTS/60)*100+MOD(IMTS,60)
      WRITE(CLFNGG,'(A4,A4,I3.3,''+'',I6.6)') CLR1INC,CLR2,GET_NSIM4D(),IINC
    ELSE
      WRITE(CLFNGG,'(A5,A4,''+'',I10.10)') CLR1,CLR2,IINC
    ENDIF
    WRITE(CLFNUA,'(A5,A4,''+'',I10.10)') CLR4,CLR2,IINC
    WRITE(CLFNOC,'(A5,A4,''+'',I10.10)') CLR5,CLR2,IINC
  ENDIF

  IF( LECFPOS .AND. NFPOS == 2 )  THEN
    IF(CLEVT == 'SFC') THEN
      IF (NSTEP == 0) THEN
        WRITE(CLFN,'(A10,A6)') CLFNGG(1:10),'000000'
      ELSE
        IF ((NYMLPPGP == 0).OR.(IYM == NYMLPPGP).OR. &
 &         ((IYM /= NYMLPPGP).AND.((IJOUR /= 1).OR.(IHR /= 0)))) THEN
          WRITE(CLFN,'(A10,I6)') CLFNGG(1:10),IYM
        ELSE
          WRITE(CLFN,'(A10,I6)') CLFNGG(1:10),NYMLPPGP
        ENDIF
      ENDIF
    ELSE
      IF (NSTEP == 0) THEN
        WRITE(CLFN,'(A10,A6)') CLFNUA(1:10),'000000'
      ELSE
        IF ((NYMLPPGP == 0).OR.(IYM == NYMLPPGP).OR. &
 &         ((IYM /= NYMLPPGP).AND.((IJOUR /= 1).OR.(IHR /= 0)))) THEN
          WRITE(CLFN,'(A10,I6)') CLFNUA(1:10),IYM
        ELSE
          WRITE(CLFN,'(A10,I6)') CLFNUA(1:10),NYMLPPGP
        ENDIF
      ENDIF
    ENDIF
  ELSE
    CLFN = CLFNGG
  ENDIF

! CLFN = CLFNGG


  IF(CLEVT == 'OML') THEN    !KPP
    CLFN = CLFNOC            !KPP
  ENDIF                      !KPP

ENDIF

IF(LTWBGV) THEN
  CLFBGV=CFNBGV
  WRITE(CLFBGV(7:),'(I3.3)') MBGVEC
  WRITE(CLFBGV(1:2),'(A2)') 'UA'
  CLFN=CLFBGV
ENDIF

IF ((NSTEP <= 0).OR.(NYMLPPGP == 0).OR. &
 & ((IYM /= NYMLPPGP).AND.((IJOUR /= 1).OR.(IHR /= 0)))) THEN
  CLMODE=CDMODE
ELSE
  CLMODE='a'
ENDIF

IF ((NSTEP <= 0).OR.((NYMLPPGP /= 0).AND. &
 & (IYM /= NYMLPPGP).AND.(IJOUR == 1).AND.(IHR == 0))) THEN
  NYMLPPGP = NYMLPPGP
ELSE
  NYMLPPGP = IYM
ENDIF

! TODO: Evaluate what this does. In 43r3 we set this switch locally. In 48r1 its and ASSOCIATE that can affect other places, and was
! not there to begin with. I leave it commented out for now - Jan St.
! CFPFMT == 'GAUSS' 

CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,KRESOL=KRESOL,&
 & KGRIB2D=KGRIBIO(1,:),KLEVS2D=KGRIBIO(2,:),CDLEVTYPE=CLEVT,PTSTEP=PTSTEP)
ZGGFLD(:,:,1)=PGGFLD(:,:)
IF(NOUTTYPE == 2) THEN
  CALL IO_PUT(Y_IOSTREAM_FDB,YL_IOREQUEST,PR3=ZGGFLD)
ELSE
#ifdef WITH_XIOS
  IF (LXIOS) THEN
    ! XIOS_FPOS extra logging
    WRITE(NULOUT, '(''XIOSFPOS: ENTERING XIOS FIELD SENDING PROCEDURES FOR GP FIELDS'')')
    CALL XIOS_PUT(YL_IOREQUEST,YDFPFIELDS,CMODE=CLMODE,PR3=ZGGFLD)
  ELSE
    !CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CLFN),CDMODE=CLMODE,KIOMASTER=1)
    CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CLFN),CDMODE='a',KIOMASTER=1)
    CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR3=ZGGFLD)
    CALL CLOSE_IOSTREAM(YL_IOSTREAM)
  ENDIF
#else
  CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CDFN),CDMODE=CDMODE,KIOMASTER=1)
  CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR3=ZGGFLD)
  CALL CLOSE_IOSTREAM(YL_IOSTREAM)
#endif
ENDIF
 
CALL CLOSE_IOREQUEST(YL_IOREQUEST)

!Make sure everyone knows when a surface parameter was last postprocessed
IF(NPROC > 1 .AND. CLEVT == 'SFC') THEN
  CALL MPL_ALLREDUCE(NSTEPLPP(:,2),'MAX',CDSTRING='WROUTGPGB:')
ENDIF

END ASSOCIATE
END ASSOCIATE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WROUTGPGB',1,ZHOOK_HANDLE)
END SUBROUTINE WROUTGPGB
