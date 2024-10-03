MODULE suxios

USE PARKIND1, ONLY : JPIM, JPRB
USE yomxios
USE xios

IMPLICIT NONE
PRIVATE

PUBLIC :: suxios_ini, suxios_fin, suxios_ini_ctxt, suxios_fin_ctxt, suxios_ctxt, suxios_namfpc_l, suxios_namfpc_obj, &
 & suxios_namfpc_sci, suxios_namct0a, suxios_namct0b, suxios_pextra_fields

CONTAINS

SUBROUTINE suxios_ini

  USE MPL_MODULE,   ONLY : LMPLUSERCOMM, MPLUSERCOMM

  CALL GSTATS(2032,0)
  !$OMP SINGLE

  ! Enabling the usage of MPLUSERCOMM communicator, rather than MPI_COMM_WORLD 
  LMPLUSERCOMM = .TRUE.

  ! Initialization of XIOS and definition of the MPLUSERCOMM communicator to be used by IFS
  CALL xios_initialize(model_name,return_comm=MPLUSERCOMM)

  !$OMP END SINGLE
  CALL GSTATS(2032,1)

END SUBROUTINE suxios_ini

SUBROUTINE suxios_fin

  !$OMP SINGLE

  ! Finalization of XIOS and MPI
  CALL GSTATS(2033,0)
  CALL xios_finalize()
  CALL GSTATS(2033,1)

  !$OMP END SINGLE

END SUBROUTINE suxios_fin

SUBROUTINE suxios_ini_ctxt

  USE MPL_MODULE,   ONLY : MPLUSERCOMM

  !$OMP SINGLE

  ! Context initialization
  CALL GSTATS(2035,0)
  CALL xios_context_initialize(ifs_context, MPLUSERCOMM)
  CALL xios_get_handle(ifs_context, context_handle)
  CALL xios_set_current_context(context_handle)
  CALL GSTATS(2035,1)

  !$OMP END SINGLE

END SUBROUTINE suxios_ini_ctxt

SUBROUTINE suxios_fin_ctxt

  !$OMP SINGLE

  ! Deallocating XIOS buffers if needed
  IF (LOPT_SEND) THEN
    IF (LSINGLE_PREC_SEND) THEN
      IF (ALLOCATED(SFCFLDBUF_SP)) DEALLOCATE(SFCFLDBUF_SP)
      IF (ALLOCATED(MLFLDBUF_SP)) DEALLOCATE(MLFLDBUF_SP)
      IF (ALLOCATED(PLFLDBUF_SP)) DEALLOCATE(PLFLDBUF_SP)
      IF (ALLOCATED(TLFLDBUF_SP)) DEALLOCATE(TLFLDBUF_SP)
      IF (ALLOCATED(VLFLDBUF_SP)) DEALLOCATE(VLFLDBUF_SP)
    ELSE
      IF (ALLOCATED(SFCFLDBUF_DP)) DEALLOCATE(SFCFLDBUF_DP)
      IF (ALLOCATED(MLFLDBUF_DP)) DEALLOCATE(MLFLDBUF_DP)
      IF (ALLOCATED(PLFLDBUF_DP)) DEALLOCATE(PLFLDBUF_DP)
      IF (ALLOCATED(TLFLDBUF_DP)) DEALLOCATE(TLFLDBUF_DP)
      IF (ALLOCATED(VLFLDBUF_DP)) DEALLOCATE(VLFLDBUF_DP)
    END IF
  END IF

  ! Finalization of XIOS context
  CALL xios_context_finalize()

  !$OMP END SINGLE

END SUBROUTINE suxios_fin_ctxt

SUBROUTINE suxios_ctxt(YDGEOMETRY, YDMODEL)

  USE GEOMETRY_MOD, ONLY : GEOMETRY
  USE TYPE_MODEL, ONLY : MODEL

  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
  TYPE(MODEL), INTENT(IN) :: YDMODEL

  CALL GSTATS(2034,0)
  !$OMP SINGLE

  ! Date setting
  CALL GSTATS(2036,0)
  CALL ifs_xios_set_calendar(YDMODEL)
  CALL GSTATS(2036,1)

  ! Definition of axes
  CALL GSTATS(2037,0)
  CALL ifs_xios_set_axis(YDGEOMETRY)
  CALL GSTATS(2037,1)

  ! Definition of domains
  CALL GSTATS(2038,0)
  CALL ifs_xios_set_domain(YDGEOMETRY)
  CALL GSTATS(2038,1)

  ! Definition of type of communication
  CALL ifs_xios_set_type_communication

  ! Close context definition
  CALL GSTATS(2039,0)
  CALL xios_close_context_definition()
  CALL GSTATS(2039,1)

  !$OMP END SINGLE
  CALL GSTATS(2034,1)

END SUBROUTINE suxios_ctxt

SUBROUTINE ifs_xios_set_calendar(YDMODEL)

  USE TYPE_MODEL, ONLY : MODEL
  USE YOMRIP0, ONLY : NINDAT, NSSSSS
  USE YOMLUN , ONLY : NULOUT, NULRCF
  USE YOMIOS , ONLY : CFRCF
  USE YOMRES , ONLY : CSTEP

  TYPE(MODEL), INTENT(IN) :: YDMODEL

  INTEGER(KIND=JPIM) :: year, month, day, hours, minutes, seconds
  LOGICAL :: lexist
  CHARACTER (LEN=20) :: time_origin_str, start_date_str, time_step_str, duration_from_origin_str

  ! Variables not necessary for XIOS, but for reading the NAMRCF namelist
  CHARACTER (LEN=14) :: CTIME
  REAL(KIND=JPRB) :: GMASS0, GMASSI
  INTEGER(KIND=JPIM) :: NSTEPLPP(5,2)

#include "namrcf.nam.h"

  ASSOCIATE (TSTEP => YDMODEL%YRML_GCONF%YRRIP%TSTEP)

  ! Not necessary, the type of calendar is set up in the iodef.xml file
  !CALL xios_define_calendar(type="Gregorian")

  ! The current date is: start_date + NSTEP*TSTEP
  time_step%second = TSTEP
  CALL xios_set_timestep(time_step)

  ! Set up time origin of the simulation
  year = NINDAT/10000
  month = MOD(NINDAT/100, 100)
  day = MOD(NINDAT, 100)
  hours = NSSSSS/3600
  minutes = (NSSSSS - hours*3600)/60
  seconds = NSSSSS - hours*3600 - minutes*60

  time_origin = xios_date(year, month, day, hours, minutes, seconds)

  ! Time origin of the time axis. It will appear as meta-data attached to the time axis in the output file
  CALL xios_set_time_origin(time_origin=time_origin)

  ! Set up start date of the current restart
  INQUIRE(FILE=CFRCF,EXIST=lexist)

  IF (lexist) THEN
    ! It is a restart
    OPEN(NULRCF,FILE=CFRCF,DELIM="QUOTE")
    READ(NULRCF,NAMRCF)
    CLOSE(NULRCF,STATUS='KEEP')

    READ(CSTEP, *) nstep_from_origin
    
    duration_from_origin%second = REAL(nstep_from_origin,JPRB)*TSTEP
    start_date = time_origin + duration_from_origin
  ELSE
    ! It is not a restart
    start_date = time_origin
    nstep_from_origin = 0
  END IF

  ! Start date of the simulation for the current context
  CALL xios_set_start_date(start_date=start_date)

  CALL xios_date_convert_to_string(time_origin, time_origin_str)
  CALL xios_date_convert_to_string(start_date, start_date_str)
  CALL xios_duration_convert_to_string(time_step, time_step_str)
  CALL xios_duration_convert_to_string(duration_from_origin, duration_from_origin_str)

  WRITE(NULOUT,*) 'XIOSFPOS: TIME_STEP IS ',time_step_str
  WRITE(NULOUT,*) 'XIOSFPOS: TIME_ORIGIN IS ',time_origin_str
  WRITE(NULOUT,*) 'XIOSFPOS: START_TIME IS ',start_date_str
  WRITE(NULOUT, '('' XIOSFPOS: NUMBER OF TIME STEPS FROM TIME_ORIGIN IS'',I8)') nstep_from_origin
  WRITE(NULOUT,*) 'XIOSFPOS: DURATION FROM TIME_ORIGIN IS ',duration_from_origin_str

  END ASSOCIATE

END SUBROUTINE ifs_xios_set_calendar

SUBROUTINE ifs_xios_set_axis(YDGEOMETRY)

  USE GEOMETRY_MOD, ONLY : GEOMETRY
  ! XIOS_FPOS extra logging
  USE YOMLUN , ONLY : NULOUT
  
  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY

  INTEGER(KIND=JPIM) :: i
  REAL(KIND=JPRB) :: j
  REAL(KIND=JPRB), ALLOCATABLE :: zML(:)
  REAL(KIND=8),    ALLOCATABLE :: transfer_value_1d(:)

  ! Definition of model levels axis
  ALLOCATE(zML(YDGEOMETRY%YRDIMV%NFLEVG))
  ALLOCATE(transfer_value_1d(YDGEOMETRY%YRDIMV%NFLEVG))
  zML = 0.0_JPRB
  transfer_value_1d = 0.0

  j = 1.0
  DO i = 1, YDGEOMETRY%YRDIMV%NFLEVG
    zML(i) = j
    j = j + 1.0
  END DO

  ! Output all model levels
  transfer_value_1d = REAL(zML,KIND=8)
  CALL xios_set_axis_attr(model_axis_name, n_glo=YDGEOMETRY%YRDIMV%NFLEVG, value=transfer_value_1d, unit="-", positive="up")
  zML=REAL(transfer_value_1d, KIND=JPRB)
  DEALLOCATE(zML)
  DEALLOCATE(transfer_value_1d)

END SUBROUTINE ifs_xios_set_axis

SUBROUTINE ifs_xios_set_domain(YDGEOMETRY)

  USE GEOMETRY_MOD, ONLY : GEOMETRY
  USE YOMMP0,       ONLY : MYPROC, MYSETA, MYSETB
  USE YOMCST,       ONLY : RPI

  TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY

  INTEGER(KIND=JPIM) :: i, j, ndglg, ni_glo, ni, nvertex
  REAL(KIND=JPRB) :: zdeltax, zdeltayup, zdeltaydw
  INTEGER(KIND=JPIM), ALLOCATABLE :: i_index(:)
  REAL(KIND=JPRB), ALLOCATABLE :: lonvalue_1d(:), latvalue_1d(:)
  REAL(KIND=JPRB), ALLOCATABLE :: bounds_lon_1d(:,:), bounds_lat_1d(:,:)
  REAL(KIND=JPRB), ALLOCATABLE :: zrgauslat(:)
  REAL(KIND=8), ALLOCATABLE    :: transfer_lat_1d(:)
  REAL(KIND=8), ALLOCATABLE    :: transfer_lon_1d(:)
  REAL(KIND=8), ALLOCATABLE    :: transfer_lat_2d(:,:)
  REAL(KIND=8), ALLOCATABLE    :: transfer_lon_2d(:,:)

  ni_glo = YDGEOMETRY%YRGEM%NGPTOTG
  ni = YDGEOMETRY%YRGEM%NGPTOT
  ndglg = YDGEOMETRY%YRDIM%NDGLG
  nvertex = 4

  ALLOCATE(i_index(ni))
  ALLOCATE(lonvalue_1d(ni), latvalue_1d(ni))
  ALLOCATE(transfer_lon_1d(ni), transfer_lat_1d(ni))
  lonvalue_1d = 0.0_JPRB
  latvalue_1d = 0.0_JPRB
  transfer_lon_1d = 0.0
  transfer_lat_1d = 0.0
  ALLOCATE(bounds_lon_1d(nvertex, ni), bounds_lat_1d(nvertex, ni)) 
  ALLOCATE(transfer_lon_2d(nvertex, ni), transfer_lat_2d(nvertex, ni))
  bounds_lon_1d = 0.0_JPRB
  bounds_lat_1d = 0.0_JPRB
  transfer_lon_2d = 0.0
  transfer_lat_2d = 0.0
  ALLOCATE (zrgauslat(0:ndglg+1))
  zrgauslat = 0.0_JPRB

  !
  !* Local domain data
  !
  ! Re-write a more efficient version by using NDGLG, NLOENG, NSTA and NONL
  j = 0
  DO i = 1, ni_glo
    IF (YDGEOMETRY%YRMP%NGLOBALPROC(i) == MYPROC) THEN
      j = j + 1
      ! XIOS requires indexing from 0
      i_index(j) = i - 1
    END IF
  END DO

  zrgauslat(0)       = 90.0_JPRB
  zrgauslat(ndglg+1) = -90.0_JPRB
  zrgauslat(1:ndglg) = ASIN(YDGEOMETRY%YRCSGLEG%RMU(1:ndglg))*(180.0_JPRB/(RPI))

  DO i = 1, ni
    !
    !* Longitudes and latitudes for local domain grid-points (from radians to degrees)
    !
    latvalue_1d(i) = REAL(YDGEOMETRY%YRGSGEOM_NB%GELAT(i)*(180.0_JPRB/RPI),JPRB)
    lonvalue_1d(i) = REAL(YDGEOMETRY%YRGSGEOM_NB%GELAM(i)*(180.0_JPRB/RPI),JPRB)

    !
    !* Cells' boundaries for local domain grid-points
    !
    zdeltax =   0.5_JPRB*360.0_JPRB/REAL(YDGEOMETRY%YRGEM%NLOENG(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i)))
    zdeltayup = 0.5_JPRB*(zrgauslat(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i) - 1) - zrgauslat(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i)))
    zdeltaydw = 0.5_JPRB*(zrgauslat(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i) + 1) - zrgauslat(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i)))

    IF (zrgauslat(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i) - 1) == 90.0_JPRB) zdeltayup = 2.0_JPRB*zdeltayup
    IF (zrgauslat(YDGEOMETRY%YRGSGEOM_NB%NGPLAT(i) + 1) == -90.0_JPRB) zdeltaydw = 2.0_JPRB*zdeltaydw

    bounds_lon_1d(1,i) = lonvalue_1d(i) + zdeltax
    bounds_lat_1d(1,i) = latvalue_1d(i) + zdeltaydw
    bounds_lon_1d(2,i) = lonvalue_1d(i) + zdeltax
    bounds_lat_1d(2,i) = latvalue_1d(i) + zdeltayup
    bounds_lon_1d(3,i) = lonvalue_1d(i) - zdeltax
    bounds_lat_1d(3,i) = latvalue_1d(i) + zdeltayup
    bounds_lon_1d(4,i) = lonvalue_1d(i) - zdeltax
    bounds_lat_1d(4,i) = latvalue_1d(i) + zdeltaydw
  END DO

  ! Define global domain and local domain index
  CALL xios_set_domain_attr(gaussian_domain_name, type='gaussian', ni_glo=ni_glo, ibegin=0, ni=ni, i_index=i_index)
  ! Define local domain data
  CALL xios_set_domain_attr(gaussian_domain_name, data_dim=1, data_ibegin=0, data_ni=ni)
  ! Define longitudes and latitudes for grid-point cells
  transfer_lon_1d = REAL(lonvalue_1d,KIND=8)
  transfer_lat_1d = REAL(latvalue_1d,KIND=8)
  CALL xios_set_domain_attr(gaussian_domain_name, lonvalue_1d=transfer_lon_1d, latvalue_1d=transfer_lat_1d)
  lonvalue_1d = REAL(transfer_lon_1d,KIND=JPRB)
  latvalue_1d = REAL(transfer_lat_1d,KIND=JPRB)
  ! Define cell's boundaries
  transfer_lon_2d = REAL(bounds_lon_1d,KIND=8)
  transfer_lat_2d = REAL(bounds_lat_1d,KIND=8)
  CALL xios_set_domain_attr(gaussian_domain_name, nvertex=nvertex, bounds_lon_1d=transfer_lon_2d, bounds_lat_1d=transfer_lat_2d)
  bounds_lon_1d = REAL(transfer_lon_2d,KIND=JPRB)
  bounds_lat_1d = REAL(transfer_lat_2d,KIND=JPRB)
  DEALLOCATE(i_index)
  DEALLOCATE(lonvalue_1d, latvalue_1d)
  DEALLOCATE(bounds_lon_1d, bounds_lat_1d)
  DEALLOCATE(transfer_lon_1d, transfer_lat_1d)
  DEALLOCATE(transfer_lon_2d, transfer_lat_2d)
  DEALLOCATE(zrgauslat)

END SUBROUTINE ifs_xios_set_domain

SUBROUTINE ifs_xios_set_type_communication

  ! XIOS_FPOS extra logging
  USE YOMLUN, ONLY : NULOUT

  ! Setting whether using delayed (optimized) send or not
  IF (xios_getvar(lopt_send_var_name, LOPT_SEND)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: LOPT_SEND IS'',L2)') LOPT_SEND
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS WARNING: LOPT_SEND IS NOT DEFINED IN IODEF.XML. BY DEFAULT IS FALSE'')')
  END IF
  ! Setting whether using single precision to send data to XIOS
  IF (xios_getvar(lsingle_prec_send_var_name, LSINGLE_PREC_SEND)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: LSINGLE_PREC_SEND IS'',L2)') LSINGLE_PREC_SEND
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS WARNING: LSINGLE_PREC_SEND IS NOT DEFINED IN IODEF.XML. BY DEFAULT IS FALSE'')')
  END IF

END SUBROUTINE ifs_xios_set_type_communication



SUBROUTINE suxios_namfpc_sci(YDNAMFPSCI)

  USE YOMFPC, ONLY : TNAMFPSCI, LTRACEFP
  USE YOMLUN, ONLY : NULOUT

  TYPE(TNAMFPSCI), INTENT(OUT)  :: YDNAMFPSCI

  ASSOCIATE(NFITP=>YDNAMFPSCI%NFITP, NFITT=>YDNAMFPSCI%NFITT, NFITV=>YDNAMFPSCI%NFITV, &
    & NFPCLI=>YDNAMFPSCI%NFPCLI, LFPQ=>YDNAMFPSCI%LFPQ, RFPCORR=>YDNAMFPSCI%RFPCORR)

  ! Setting spectral fitting and other FullPos variables
  IF (xios_getvar(nfitp_var_name, NFITP)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFITP IS'',I4)') NFITP
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFITP IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(nfitt_var_name, NFITT)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFITT IS'',I4)') NFITT
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFITT IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(nfitv_var_name, NFITV)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFITV IS'',I4)') NFITV
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFITV IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(nfpcli_var_name, NFPCLI)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFPCLI IS'',I4)') NFPCLI
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFPCLI IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(lfpq_var_name, LFPQ)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: LFPQ IS'',L2)') LFPQ
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: LFPQ IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(ltracefp_var_name, LTRACEFP)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: LTRACEFP IS'',L2)') LTRACEFP
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: LTRACEFP IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(rfpcorr_var_name, RFPCORR)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: RFPCORR IS'',F8.1)') RFPCORR
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: RFPCORR IS NOT DEFINED'')')
  END IF
  YDNAMFPSCI%NFITI=0
  YDNAMFPSCI%NFITS=2
  YDNAMFPSCI%LFPRH100=.TRUE.


  END ASSOCIATE

END SUBROUTINE suxios_namfpc_sci



SUBROUTINE suxios_namfpc_obj(YDNAMFPOBJ)

  USE YOMFPC, ONLY : TNAMFPOBJ
  USE YOMLUN, ONLY : NULOUT

  TYPE(TNAMFPOBJ), INTENT(OUT)  :: YDNAMFPOBJ

  ASSOCIATE(CFPFMT=>YDNAMFPOBJ%CFPFMT, NFRFPOS=>YDNAMFPOBJ%NFRFPOS)

  CFPFMT = 'MODEL'
  ! For some reason IFS has two variables that are described to do the same thing:
  ! NFRFPOS and NFRPOS both control the FULLPOS output frequency.
  ! We write NFRFPOS from XIOS NFRPOS.
  IF (xios_getvar(nfrpos_var_name, NFRFPOS)) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFRFPOS IS'',I6)') NFRFPOS
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFRFPOS IS NOT DEFINED'')')
  END IF

  END ASSOCIATE

END SUBROUTINE suxios_namfpc_obj



SUBROUTINE suxios_namfpc_l(YDGEOMETRY,YDNAMFPL)

  USE GEOMETRY_MOD, ONLY : GEOMETRY
  USE YOMFPC,       ONLY : TNAMFPL
  USE PARFPOS,      ONLY : JPOS3DF, JPOS3S, JPOS3P, JPOS3TH, JPOS3PV, JPOS3H, JPOSSGP, JPOS2DF
  USE YOMLUN,       ONLY : NULOUT

  TYPE(GEOMETRY) , INTENT(IN)   :: YDGEOMETRY
  TYPE(TNAMFPL)  , INTENT(OUT)  :: YDNAMFPL

  INTEGER(KIND=JPIM) :: n_glo_ml, n_glo_pl, n_glo_th, n_glo_pv, n_glo_hl, i
  CHARACTER(LEN=16) :: cgrb
  LOGICAL :: lpost = .false.
  REAL(KIND=8),ALLOCATABLE :: transfer_value_2d(:)

  !$OMP SINGLE

  ASSOCIATE(RFP3P=>YDNAMFPL%RFP3P, RFP3H=>YDNAMFPL%RFP3H, RFP3TH=>YDNAMFPL%RFP3TH, &
    & RFP3PV=>YDNAMFPL%RFP3PV, NRFP3S=>YDNAMFPL%NRFP3S, NFP3DFP=>YDNAMFPL%NFP3DFP, &
    & NFP3DFH=>YDNAMFPL%NFP3DFH, NFP3DFT=>YDNAMFPL%NFP3DFT, NFP3DFV=>YDNAMFPL%NFP3DFV, &
    & NFP3DFS=>YDNAMFPL%NFP3DFS, MFP3DFP=>YDNAMFPL%MFP3DFP, &
    & MFP3DFH=>YDNAMFPL%MFP3DFH, MFP3DFT=>YDNAMFPL%MFP3DFT, MFP3DFV=>YDNAMFPL%MFP3DFV, &
    & MFP3DFS=>YDNAMFPL%MFP3DFS, CFP2DF=>YDNAMFPL%CFP2DF, MFP2DF=>YDNAMFPL%MFP2DF, &
    & NFP2DF=>YDNAMFPL%NFP2DF, MFPPHY=>YDNAMFPL%MFPPHY, NFPPHY=>YDNAMFPL%NFPPHY)


  
  !! RFP3I=>YDNAMFPL%RFP3I, RFP3F=>YDNAMFPL%RFP3F, NFP3DFI=>YDNAMFPL%NFP3DFI, 
  !! NFP3DFF=>YDNAMFPL%NFP3DFF, MFP3DFI=>YDNAMFPL%MFP3DFI, MFP3DFF=>YDNAMFPL%MFP3DFF
  !! CFP3DF=>YDNAMFPL%CFP3DF, CFPCFU=>YDNAMFPL%CFPCFU, CFPXFU=>YDNAMFPL%CFPXFU
  !! CFPPHY=>YDNAMFPL%CFPPHY

  ! Set variables of 3D fields
  ! Model levels
  NFP3DFS = 0
  DO i = 1, JPOS3DF 
    MFP3DFS(i) = 0
  END DO
  !CALL xios_get_axis_attr(model_axis_name, n_glo=n_glo_ml)
  n_glo_ml = YDGEOMETRY%YRDIMV%NFLEVG
  IF (n_glo_ml > 0) THEN
    DO i = 1, N3DFLD
      cgrb = C3DFLD(i)
      IF (xios_is_valid_field(TRIM(cgrb))) THEN
        IF (xios_field_is_active(TRIM(cgrb))) THEN
          NFP3DFS = NFP3DFS + 1
          MFP3DFS(NFP3DFS) = IGRB3DFLD(i)
        END IF
      ELSE
        WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' MODEL LEVELS FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
      END IF
    END DO

    ! Add PEXTRA fields to be output in model levels
    DO i = 1, NVEXTR
      MFP3DFS(NFP3DFS+i) = NVEXTRAGB(i)
    END DO
    NFP3DFS = NFP3DFS + NVEXTR

    ! Set number of full model levels for PEXTRA fields
    IF (NVEXTR > 0) NCEXTR = YDGEOMETRY%YRDIMV%NFLEVG

    DO i = 1, JPOS3S
      NRFP3S(i) = -9
    END DO
    !If used, convert returned floating point values to integer
    !CALL xios_get_axis_attr(model_axis_name, value=NRFP3S(1:n_glo_ml))
    DO i = 1, n_glo_ml
      NRFP3S(i) = i
    END DO
  END IF

  IF (n_glo_ml > 0 .and. NFP3DFS > 0) THEN
    ! Allocating XIOS buffer
    IF (LOPT_SEND) THEN
      IF (LSINGLE_PREC_SEND) THEN
        ALLOCATE(MLFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_ml,NFP3DFS))
#if defined init_alloc_zero
        MLFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        MLFLDBUF_SP = HUGE(MLFLDBUF_SP)
#endif
      ELSE
        ALLOCATE(MLFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_ml,NFP3DFS))
#if defined init_alloc_zero
        MLFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        MLFLDBUF_DP = HUGE(MLFLDBUF_DP)
#endif
      END IF
    END IF
    WRITE(NULOUT, '(''XIOSFPOS: NFP3DFS IS'',I4)') NFP3DFS
    DO i = 1, NFP3DFS
      WRITE(NULOUT, '(''XIOSFPOS: MFP3DFS('',I4,'') IS'',I8)') i, MFP3DFS(i)
    END DO
    WRITE(NULOUT, '(''XIOSFPOS: NRFP3S LENGTH IS'',I4)') n_glo_ml
    DO i = 1, n_glo_ml
       WRITE(NULOUT, '(''XIOSFPOS: NRFP3S('',I4,'') IS'',I4)') i, NRFP3S(i)
    END DO
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: MODEL LEVELS NOT USED'')')
  END IF

  ! Pressure levels
  NFP3DFP = 0
  DO i = 1, JPOS3DF
    MFP3DFP(i) = 0
  END DO
  CALL xios_get_axis_attr(pressure_axis_name, n_glo=n_glo_pl)
  IF (n_glo_pl > 0) THEN
    DO i = 1, N3DFLD
      IF (TRIM(C3DFLD(i)) == TRIM('pres')) CYCLE
      cgrb = TRIM(C3DFLD(i))//'_pl'
      IF (xios_is_valid_field(TRIM(cgrb))) THEN
        IF (xios_field_is_active(TRIM(cgrb))) THEN
          NFP3DFP = NFP3DFP + 1
          MFP3DFP(NFP3DFP) = IGRB3DFLD(i)
        END IF
      ELSE
        WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' PRESSURE LEVELS FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
      END IF
    END DO
    DO i = 1, JPOS3P
      RFP3P(i) = -9._JPRB
    END DO
    ALLOCATE(transfer_value_2d(1:n_glo_pl))
    transfer_value_2d = REAL(RFP3P, KIND=8)
    CALL xios_get_axis_attr(pressure_axis_name, value=transfer_value_2d(1:n_glo_pl))
    RFP3P(1:n_glo_pl) = REAL(transfer_value_2d, KIND=JPRB)
    DEALLOCATE(transfer_value_2d)
  END IF

  IF (n_glo_pl > 0 .and. NFP3DFP > 0) THEN
    lpost = .true.
    ! Allocating XIOS buffer
    IF (LOPT_SEND) THEN
      IF (LSINGLE_PREC_SEND) THEN
        ALLOCATE(PLFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_pl,NFP3DFP))
#if defined init_alloc_zero
        PLFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        PLFLDBUF_SP = HUGE(PLFLDBUF_SP)
#endif
      ELSE
        ALLOCATE(PLFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_pl,NFP3DFP))
#if defined init_alloc_zero
        PLFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        PLFLDBUF_DP = HUGE(PLFLDBUF_SP)
#endif
      END IF
    END IF
    WRITE(NULOUT, '(''XIOSFPOS: NFP3DFP IS'',I4)') NFP3DFP
    DO i = 1, NFP3DFP
      WRITE(NULOUT, '(''XIOSFPOS: MFP3DFP('',I4,'') IS'',I8)') i, MFP3DFP(i)
    END DO
    WRITE(NULOUT, '(''XIOSFPOS: RFP3P LENGTH IS'',I4)') n_glo_pl
    DO i = 1, n_glo_pl
      WRITE(NULOUT, '(''XIOSFPOS: RFP3P('',I4,'') IS'',F10.1)') i, RFP3P(i)
    END DO
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: PRESSURE LEVELS NOT USED'')')
  END IF

  ! Theta levels
  NFP3DFT = 0
  DO i = 1, JPOS3DF
    MFP3DFT(i) = 0
  END DO
  CALL xios_get_axis_attr(theta_axis_name, n_glo=n_glo_th)
  IF (n_glo_th > 0) THEN
    DO i = 1, N3DFLD
      IF (TRIM(C3DFLD(i)) == TRIM('pt')) CYCLE
      cgrb = TRIM(C3DFLD(i))//'_th'
      IF (xios_is_valid_field(TRIM(cgrb))) THEN
        IF (xios_field_is_active(TRIM(cgrb))) THEN
          NFP3DFT = NFP3DFT + 1
          MFP3DFT(NFP3DFT) = IGRB3DFLD(i)
        END IF
      ELSE
        WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' THETA LEVELS FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
      END IF
    END DO
    DO i = 1, JPOS3TH
      RFP3TH(i) = -9._JPRB
    END DO
    ALLOCATE(transfer_value_2d(1:n_glo_th))
    CALL xios_get_axis_attr(theta_axis_name, value=transfer_value_2d(1:n_glo_th))
    RFP3TH(1:n_glo_th)=REAL(transfer_value_2d, KIND=JPRB)
    DEALLOCATE(transfer_value_2d)
  END IF

  IF (n_glo_th > 0 .and. NFP3DFT > 0) THEN
    lpost = .true.
    ! Allocating XIOS buffer
    IF (LOPT_SEND) THEN
      IF (LSINGLE_PREC_SEND) THEN
        ALLOCATE(TLFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_th,NFP3DFT))
#if defined init_alloc_zero
        TLFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        TLFLDBUF_SP = HUGE(TLFLDBUF_SP)
#endif
      ELSE
        ALLOCATE(TLFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_th,NFP3DFT))
#if defined init_alloc_zero
        TLFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        TLFLDBUF_DP = HUGE(TLFLDBUF_SP)
#endif
      END IF
    END IF
    WRITE(NULOUT, '(''XIOSFPOS: NFP3DFT IS'',I4)') NFP3DFT
    DO i = 1, NFP3DFT
      WRITE(NULOUT, '(''XIOSFPOS: MFP3DFT('',I4,'') IS'',I8)') i, MFP3DFT(i)
    END DO
    WRITE(NULOUT, '(''XIOSFPOS: RFP3TH LENGTH IS'',I4)') n_glo_th
    DO i = 1, n_glo_th
      WRITE(NULOUT, '(''XIOSFPOS: RFP3TH('',I4,'') IS'',F8.1)') i, RFP3TH(i)
    END DO
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: THETA LEVELS NOT USED'')')
  END IF

  ! PV levels
  NFP3DFV = 0
  DO i = 1, JPOS3DF
    MFP3DFV(i) = 0
  END DO
  CALL xios_get_axis_attr(pv_axis_name, n_glo=n_glo_pv)
  IF (n_glo_pv > 0) THEN
    DO i = 1, N3DFLD
      IF (TRIM(C3DFLD(i)) == TRIM('pv')) CYCLE
      cgrb = TRIM(C3DFLD(i))//'_pv'
      IF (xios_is_valid_field(TRIM(cgrb))) THEN
        IF (xios_field_is_active(TRIM(cgrb))) THEN
          NFP3DFV = NFP3DFV + 1
          MFP3DFV(NFP3DFV) = IGRB3DFLD(i)
        END IF
      ELSE
        WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' PV LEVELS FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
      END IF
    END DO
    DO i = 1, JPOS3PV
      RFP3PV(i) = 9999*1.E6_JPRB
    END DO
    ALLOCATE(transfer_value_2d(1:n_glo_pv))
    CALL xios_get_axis_attr(pv_axis_name, value=transfer_value_2d(1:n_glo_pv))
    RFP3PV(1:n_glo_pv)=REAL(transfer_value_2d,KIND=JPRB)
    DEALLOCATE(transfer_value_2d)
  END IF

  IF (n_glo_pv > 0 .and. NFP3DFV > 0) THEN
    lpost = .true.
    ! Allocating XIOS buffer
    IF (LOPT_SEND) THEN
      IF (LSINGLE_PREC_SEND) THEN
        ALLOCATE(VLFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_pv,NFP3DFV))
#if defined init_alloc_zero
        VLFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        VLFLDBUF_SP = HUGE(VLFLDBUF_SP)
#endif
      ELSE
        ALLOCATE(VLFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,n_glo_pv,NFP3DFV))
#if defined init_alloc_zero
        VLFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        VLFLDBUF_DP = HUGE(VLFLDBUF_SP)
#endif
      END IF
    END IF
    WRITE(NULOUT, '(''XIOSFPOS: NFP3DFV IS'',I4)') NFP3DFV
    DO i = 1, NFP3DFV
      WRITE(NULOUT, '(''XIOSFPOS: MFP3DFV('',I4,'') IS'',I8)') i, MFP3DFV(i)
    END DO
    WRITE(NULOUT, '(''XIOSFPOS: RFP3PV LENGTH IS'',I4)') n_glo_pv
    DO i = 1, n_glo_pv
      WRITE(NULOUT, '(''XIOSFPOS: RFP3PV('',I4,'') IS'',F12.8)') i, RFP3PV(i)
    END DO
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: PV LEVELS NOT USED'')')
  END IF

  ! Height levels
  !NFP3DFH = 0
  !DO i = 1, JPOS3DF
  !  MFP3DFH(i) = 0
  !END DO
  !CALL xios_get_axis_attr(height_axis_name, n_glo=n_glo_hl)
  !IF (n_glo_hl > 0) THEN
  !  DO i = 1, N3DFLD
  !    IF (TRIM(C3DFLD(i)) == TRIM('z')) CYCLE
  !    cgrb = TRIM(C3DFLD(i))//'_hl'
  !    IF (xios_is_valid_field(TRIM(cgrb))) THEN
  !      IF (xios_field_is_active(TRIM(cgrb))) THEN
  !        NFP3DFH = NFP3DFH + 1
  !        MFP3DFH(NFP3DFH) = IGRB3DFLD(i)
  !      END IF
  !    ELSE
  !      WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' HEIGHT LEVELS FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
  !    END IF
  !  END DO
  !  DO i = 1, JPOS3H
  !    RFP3H(i) = -9._JPRB
  !  END DO
  !  CALL xios_get_axis_attr(height_axis_name, value=RFP3H(1:n_glo_hl))
  !END IF

  !IF (n_glo_hl > 0 .and. NFP3DFH > 0) THEN
  !  lpost = .true.
  !  WRITE(NULOUT, '(''XIOSFPOS: NFP3DFH IS'',I4)') NFP3DFH
  !  DO i = 1, NFP3DFH
  !    WRITE(NULOUT, '(''XIOSFPOS: MFP3DFH('',I4,'') IS'',I8)') i, MFP3DFH(i)
  !  END DO
  !  WRITE(NULOUT, '(''XIOSFPOS: RFP3H LENGTH IS'',I4)') n_glo_hl
  !  DO i = 1, n_glo_hl
  !    WRITE(NULOUT, '(''XIOSFPOS: RFP3H('',I4,'') IS'',F8.1)') i, RFP3H(i)
  !  END DO
  !ELSE
  !  WRITE(NULOUT, '(''XIOSFPOS: HEIGHT LEVELS NOT USED'')')
  !END IF

  ! Set 2D physical fields
  NFPPHY = 0
  DO i = 1, JPOSSGP
    MFPPHY(i) = 0
  END DO
  DO i = 1, NSFCFLD
    cgrb = CSFCFLD(i)
    IF (TRIM(cgrb) == TRIM('lnsp')) CYCLE
    IF (xios_is_valid_field(TRIM(cgrb))) THEN
      IF (xios_field_is_active(TRIM(cgrb))) THEN
        NFPPHY = NFPPHY + 1
        MFPPHY(NFPPHY) = IGRBSFCFLD(i)
      END IF
    ELSE
      WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' SURFACE PHYSICAL FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
    END IF
  END DO

  IF (NFPPHY > 0) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFPPHY IS'',I4)') NFPPHY
    DO i = 1, NFPPHY
      WRITE(NULOUT, '(''XIOSFPOS: MFPPHY('',I4,'') IS'',I8)') i, MFPPHY(i)
    END DO
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: SURFACE PHYSICAL FIELDS NOT USED'')')
  END IF

  ! Set 2D dynamical fields
  NFP2DF = 0
  DO i = 1, JPOS2DF
    MFP2DF(i) = 0
  END DO
  ! "It is strongly recommended these fields are always enabled in the output in order for the 
  ! FullPos post-processing to work correctly"
  IF (lpost) THEN
    NFP2DF = 3
    MFP2DF(1) = 129
    MFP2DF(2) = 134
    MFP2DF(3) = 152
  END IF
  ! Logarithm of surface pressure (lnsp <-> 152) not available in FullPos as a grid-point field,
  ! so it necessary to treat it separately from the rest of surface fields. This means to transform
  ! it from spectral space to grid-point space before sending it to XIOS
  IF (xios_is_valid_field(TRIM('lnsp'))) THEN
    IF (xios_field_is_active(TRIM('lnsp')) .and. (.not. lpost)) THEN
      NFP2DF = 1
      MFP2DF(1) = 152
    END IF
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' SURFACE DYNAMICAL FIELD NOT DEFINED IN IODEF.XML'')') TRIM('lnsp')
  END IF

  IF (NFP2DF > 0) THEN
    WRITE(NULOUT, '(''XIOSFPOS: NFP2DF IS'',I4)') NFP2DF
    DO i = 1, NFP2DF
      WRITE(NULOUT, '(''XIOSFPOS: MFP2DF('',I4,'') IS'',I8)') i, MFP2DF(i)
    END DO
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: SURFACE DYNAMICAL FIELDS NOT USED'')')
  END IF

  ! Allocating XIOS buffer for surface fields
  IF (LOPT_SEND) THEN
    IF (LSINGLE_PREC_SEND) THEN
      IF (NFPPHY > 0 .and. NFP2DF > 0) THEN
        ALLOCATE(SFCFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,NFPPHY+1))
#if defined init_alloc_zero
        SFCFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        SFCFLDBUF_SP = HUGE(SFCFLDBUF_SP)
#endif
      ELSE IF (NFPPHY > 0) THEN
        ALLOCATE(SFCFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,NFPPHY))
#if defined init_alloc_zero
        SFCFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        SFCFLDBUF_SP = HUGE(SFCFLDBUF_SP)
#endif
      ELSE IF (NFP2DF > 0) THEN
        ALLOCATE(SFCFLDBUF_SP(YDGEOMETRY%YRGEM%NGPTOT,1))
#if defined init_alloc_zero
        SFCFLDBUF_SP = 0.0_JPRB
#elif defined init_alloc_huge
        SFCFLDBUF_SP = HUGE(SFCFLDBUF_SP)
#endif
      END IF
    ELSE
      IF (NFPPHY > 0 .and. NFP2DF > 0) THEN
        ALLOCATE(SFCFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,NFPPHY+1))
#if defined init_alloc_zero
        SFCFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        SFCFLDBUF_DP = HUGE(SFCFLDBUF_DP)
#endif
      ELSE IF (NFPPHY > 0) THEN
        ALLOCATE(SFCFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,NFPPHY))
#if defined init_alloc_zero
        SFCFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        SFCFLDBUF_DP = HUGE(SFCFLDBUF_DP)
#endif
      ELSE IF (NFP2DF > 0) THEN
        ALLOCATE(SFCFLDBUF_DP(YDGEOMETRY%YRGEM%NGPTOT,1))
#if defined init_alloc_zero
        SFCFLDBUF_DP = 0.0_JPRB
#elif defined init_alloc_huge
        SFCFLDBUF_DP = HUGE(SFCFLDBUF_DP)
#endif
      END IF
    END IF
  END IF

  ! Ensuring we use 'MODEL' as the format of the output files


  END ASSOCIATE

  !$OMP END SINGLE

END SUBROUTINE suxios_namfpc_l

SUBROUTINE suxios_namct0a

  USE YOMCT0, ONLY : JPNPST, NFPOS, NPOSTS, NHISTS

  INTEGER(KIND=JPIM) :: i

  !$OMP SINGLE

  ! Ensure FullPos is enabled
  NFPOS = 2

  ! Ensure we use regular output
  DO i = 0, JPNPST
    NPOSTS(i) = 0
    NHISTS(i) = 0
  END DO
  
  !$OMP END SINGLE

END SUBROUTINE suxios_namct0a

SUBROUTINE suxios_namct0b(YDMODEL)

  USE TYPE_MODEL, ONLY : MODEL
  USE YOMCT0, ONLY : NFRPOS, NFRHIS
  USE YOMARG, ONLY : NSUPERSEDE
  ! XIOS_FPOS extra logging
  USE YOMLUN, ONLY : NULOUT

  TYPE(MODEL), INTENT(IN) :: YDMODEL

  REAL(KIND=JPRB) :: ZUNIT

  !$OMP SINGLE
  ASSOCIATE (TSTEP => YDMODEL%YRML_GCONF%YRRIP%TSTEP)

  ! Setting FullPos output frequency
  ZUNIT=3600._JPRB
  IF (xios_getvar(nfrpos_var_name, NFRPOS)) THEN
    IF ((NSUPERSEDE > 0).AND.(TSTEP > 0.0_JPRB).AND.(NFRPOS < 0)) THEN
      NFRPOS=NINT((REAL(-NFRPOS,JPRB)*ZUNIT)/TSTEP)
    END IF
    WRITE(NULOUT, '(''XIOSFPOS: NFRPOS IS'',I4)') NFRPOS
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFRPOS IS NOT DEFINED'')')
  END IF
  IF (xios_getvar(nfrhis_var_name, NFRHIS)) THEN
    IF ((NSUPERSEDE > 0).AND.(TSTEP > 0.0_JPRB).AND.(NFRHIS < 0)) THEN
      NFRHIS=NINT((REAL(-NFRHIS,JPRB)*ZUNIT)/TSTEP)
    END IF
    WRITE(NULOUT, '(''XIOSFPOS: NFRHIS IS'',I4)') NFRHIS
  ELSE
    WRITE(NULOUT, '(''XIOSFPOS: NFRHIS IS NOT DEFINED'')')
  END IF

  END ASSOCIATE
  !$OMP END SINGLE

END SUBROUTINE suxios_namct0b

SUBROUTINE suxios_pextra_fields(IFS_NAMELIST)

  USE YOEPHY   , ONLY : YREPHY
  USE YOMDPHY  , ONLY : YRDPHY
  USE YOMPHYDS , ONLY : TPHYDS
  ! XIOS_FPOS extra logging
  USE YOMLUN, ONLY : NULOUT

  CHARACTER(LEN=*), INTENT(IN) :: IFS_NAMELIST

  TYPE(TPHYDS) :: YRPHYDS

  INTEGER(KIND=JPIM) :: i
  CHARACTER(LEN=16)  :: cgrb

  !$OMP SINGLE

  IF (IFS_NAMELIST == 'NAEPHY') THEN
    ! Set variables of PEXTRA fields
    NVEXTR = 0
    NVEXTRAGB(:) = -999
    DO i = 1, NPEXTRAFLD
      cgrb = CPEXTRAFLD(i)
      IF (xios_is_valid_field(TRIM(cgrb))) THEN
        IF (xios_field_is_active(TRIM(cgrb))) THEN
          NVEXTR = NVEXTR + 1
          NVEXTRAGB(NVEXTR) = IGRBPEXTRAFLD(i)
        END IF
      ELSE
        WRITE(NULOUT, '(''XIOSFPOS WARNING:'',A16,'' PEXTRA FIELD NOT DEFINED IN IODEF.XML'')') TRIM(cgrb)
      END IF
    END DO

    IF (NVEXTR > 0) THEN
      LBUD23 = .TRUE.
      WRITE(NULOUT, '(''XIOSFPOS: PHYSICAL TENDENCIES AND FLUXES OUTPUT (PEXTRA FIELDS) IS ENABLED'')')
    ELSE
      LBUD23 = .FALSE.
      WRITE(NULOUT, '(''XIOSFPOS: PHYSICAL TENDENCIES AND FLUXES OUTPUT (PEXTRA FIELDS) IS NOT ENABLED'')')
    END IF

    ! Enable (or not) computation of physics tendencies (PEXTRA fields)
    YREPHY%LBUD23 = LBUD23
  ELSE IF (IFS_NAMELIST == 'NAMDPHY') THEN
    ! Set number of tendency output fields (PEXTRA fields)
    YRDPHY%NVEXTR = NVEXTR
    WRITE(NULOUT, '(''XIOSFPOS: NVEXTR IS'',I4)') NVEXTR
    ! Set number of full model levels e.g. 60, 91, 137, etc.
    YRDPHY%NCEXTR = NCEXTR
    WRITE(NULOUT, '(''XIOSFPOS: NCEXTR IS'',I4)') NCEXTR
  ELSE IF (IFS_NAMELIST == 'NAMPHYDS') THEN
    ! Define GRIB codes for the tendency fields (PEXTRA fields)
    YRPHYDS%NVEXTRAGB(:) = -999
    DO i = 1, NVEXTR
      YRPHYDS%NVEXTRAGB(i) = NVEXTRAGB(i)
      WRITE(NULOUT, '(''XIOSFPOS: NVEXTRAGB('',I4,'') IS'',I8)') i, NVEXTRAGB(i)
    END DO
  END IF

  !$OMP END SINGLE

END SUBROUTINE suxios_pextra_fields

END MODULE suxios
