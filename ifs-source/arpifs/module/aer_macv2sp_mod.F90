!>
!!
!! @brief Module AER_MACv2SP_MOD: provides anthropogenic aerorol optial properties as a function of lat, lon
!!   height, time, and wavelength
!!
!! Adapted from 
!!
!! MO_SIMPLE_PLUMES by @author Bjorn Stevens & Karsten Peters MPI-M, Hamburg (v1-beta release 2015-12-05; adapted to final v1 release of 2016-06-20)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Based on code originally developed at the MPI by Karsten Peters, Bjorn Stevens, Stephanie Fiedler
!!   and Stefan Kinne with input from Thorsten Mauritsen and Robert Pincus
!!
!! @par Copyright
!!
!! S. Boussetta Mars 2016: Adapted and Modified for integration within IFS  
!
MODULE AER_MACv2SP_MOD

  USE netcdf
  USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
  USE mpi
  USE ECPHYS_AUX_TYPE_MOD, ONLY : AUX_TYPE, KEYS_LOCAL_TYPE
  USE ECE_CMIP           , ONLY : LCMIP6, LCMIP7, CMIP7DATADIR, CMIP6DATADIR
  USE YOMLUN             , ONLY : NULOUT
  !  USE RADIATION_AEROSOL,        ONLY : AEROSOL_TYPE

  IMPLICIT NONE

  TYPE (AUX_TYPE)                :: PAUX

  INTEGER(KIND=JPIM), PARAMETER ::              &
       nplumes   = 9                           ,& !< Number of plumes
       nfeatures = 2                           ,& !< Number of features per plume
       ntimes    = 52                          ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251                            !< Number of years of available forcing

  REAL(KIND=JPRB), PARAMETER    ::              &
       pi      = 2.*ASIN(1.),                   & !< half the unit circle (radians)
       deg2rad = pi/180.                          !< conversion factor from degrees to radians

  LOGICAL, SAVE ::                              &
       sp_initialized = .FALSE.                   !< parameter determining whether input needs to be read

  REAL(KIND=JPRB) ::                            &
       is_biomass     (nplumes)                ,& !< if plumes is mainly from biomass
       plume_lat      (nplumes)                ,& !< latitude where plume maximizes
       plume_lon      (nplumes)                ,& !< longitude where plume maximizes
       beta_a         (nplumes)                ,& !< parameter a for beta function vertical profile
       beta_b         (nplumes)                ,& !< parameter b for beta function vertical profile
       aod_spmx       (nplumes)                ,& !< aod at 550 for simple plume (maximum)
       aod_fmbg       (nplumes)                ,& !< aod at 550 for fine mode background (for twomey effect)
       asy550         (nplumes)                ,& !< asymmetry parameter for plume at 550nm
       ssa550         (nplumes)                ,& !< single scattering albedo for plume at 550nm
       angstrom       (nplumes)                ,& !< Angstrom parameter for plume 
       sig_lon_E      (nfeatures, nplumes)     ,& !< Eastward extent of plume feature
       sig_lon_W      (nfeatures, nplumes)     ,& !< Westward extent of plume feature
       sig_lat_E      (nfeatures, nplumes)     ,& !< Southward extent of plume feature
       sig_lat_W      (nfeatures, nplumes)     ,& !< Northward extent of plume feature
       theta          (nfeatures, nplumes)     ,& !< Rotation angle of feature
       ftr_weight     (nfeatures, nplumes)     ,& !< Feature weights = (nfeatures+1) to account for BB background
       time_weight    (nfeatures, nplumes)     ,& !< Time-weights = (nfeatures+1) to account for BB background
       time_weight_bg (nfeatures, nplumes)     ,& !< Time-weights for natural background
       time_weight_ref (nfeatures, nplumes)    ,& !< Time-weights for the reference year
       year_weight (nyears, nplumes)           ,& !< Yearly weight for plume
       ann_cycle   (nfeatures, ntimes, nplumes)    !< annual cycle for feature

  PUBLIC sp_aop_profile, sp_setup


CONTAINS
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_SETUP:  This subroutine should be called at initialization to read the netcdf data that describes the simple plume
  ! climatology.  The information needs to be either read by each processor or distributed to processors.
  !
  SUBROUTINE sp_setup
    !
    ! ----------

    !USE YOERAD   , ONLY : MAC2SPFIL
    INTEGER(KIND=JPIM)    :: iret, ncid, DimID, VarID, xdmy, IFIL

    CHARACTER (LEN = 300) ::  MAC2SPFN, MAC2SPFIL
    IF (LCMIP7) THEN
        WRITE(MAC2SPFIL,*) TRIM(CMIP7DATADIR)//'/macv2sp/'//'SPv2.1_1850-2023_CMIP7.nc'
    ELSE
        WRITE(MAC2SPFIL,*) TRIM(CMIP6DATADIR)//'/'//'SPv2_1850-2020_r20241218.nc'
    END IF
    !
    ! ---------- 
    !
    IFIL = LEN_TRIM(MAC2SPFIL)
    MAC2SPFN = MAC2SPFIL(1:IFIL)
    
    WRITE(NULOUT, *) "AER_MACv2SP_MOD: OPENING ",MAC2SPFN 
    iret = nf90_open(MAC2SPFN, NF90_NOWRITE, ncid)
    IF (iret /= NF90_NOERR) STOP 'NetCDF File not opened: '//MAC2SPFN 
    !
    ! read dimensions and make sure file conforms to expected size
    !
    iret = nf90_inq_dimid(ncid, "plume_number"  , DimId)
    iret = nf90_inquire_dimension(ncid, DimId, len = xdmy)
    IF (xdmy /= nplumes) STOP 'NetCDF improperly dimensioned-- plume_number'

    iret = nf90_inq_dimid(ncid, "plume_feature", DimId)
    iret = nf90_inquire_dimension(ncid, DimId, len = xdmy)
    IF (xdmy /= nfeatures) STOP 'NetCDF improperly dimensioned-- plume_feature'

    iret = nf90_inq_dimid(ncid, "year_fr"   , DimId)
    iret = nf90_inquire_dimension(ncid, DimID, len = xdmy)
    IF (xdmy /= ntimes) STOP 'NetCDF improperly dimensioned-- year_fr'

    iret = nf90_inq_dimid(ncid, "years"   , DimId)
    iret = nf90_inquire_dimension(ncid, DimID, len = xdmy)
    IF (xdmy /= nyears) STOP 'NetCDF improperly dimensioned-- years'
    !
    ! read variables that define the simple plume climatology
    !
    iret = nf90_inq_varid(ncid, "is_biomass", VarId)
    iret = nf90_get_var(ncid, VarID, is_biomass(:), start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading is_biomass'
    iret = nf90_inq_varid(ncid, "plume_lat", VarId)
    iret = nf90_get_var(ncid, VarID, plume_lat(:), start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lat'
    iret = nf90_inq_varid(ncid, "plume_lon", VarId)
    iret = nf90_get_var(ncid, VarID, plume_lon(:), start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lon'
    iret = nf90_inq_varid(ncid, "beta_a"   , VarId)
    iret = nf90_get_var(ncid, VarID, beta_a(:)   , start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading beta_a'
    iret = nf90_inq_varid(ncid, "beta_b"   , VarId)
    iret = nf90_get_var(ncid, VarID, beta_b(:)   , start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading beta_b'
    iret = nf90_inq_varid(ncid, "aod_spmx" , VarId)
    iret = nf90_get_var(ncid, VarID, aod_spmx(:)  , start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading aod_spmx'
    iret = nf90_inq_varid(ncid, "aod_fmbg" , VarId)
    iret = nf90_get_var(ncid, VarID, aod_fmbg(:)  , start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading aod_fmbg'
    iret = nf90_inq_varid(ncid, "ssa550"   , VarId)
    iret = nf90_get_var(ncid, VarID, ssa550(:)  , start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading ssa550'
    iret = nf90_inq_varid(ncid, "asy550"   , VarId)
    iret = nf90_get_var(ncid, VarID, asy550(:)  , start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading asy550'
    iret = nf90_inq_varid(ncid, "angstrom" , VarId)
    iret = nf90_get_var(ncid, VarID, angstrom(:), start=(/1/), count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading angstrom'

    iret = nf90_inq_varid(ncid, "sig_lat_W"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lat_W(:,:)    , start=(/1, 1/), count=(/nfeatures, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lat_W'
    iret = nf90_inq_varid(ncid, "sig_lat_E"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lat_E(:,:)    , start=(/1, 1/), count=(/nfeatures, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lat_E'
    iret = nf90_inq_varid(ncid, "sig_lon_E"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lon_E(:,:)    , start=(/1, 1/), count=(/nfeatures, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lon_E'
    iret = nf90_inq_varid(ncid, "sig_lon_W"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lon_W(:,:)    , start=(/1, 1/), count=(/nfeatures, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lon_W'
    iret = nf90_inq_varid(ncid, "theta"         , VarId)
    iret = nf90_get_var(ncid, VarID, theta(:,:)        , start=(/1, 1/), count=(/nfeatures, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading theta'
    iret = nf90_inq_varid(ncid, "ftr_weight"    , VarId)
    iret = nf90_get_var(ncid, VarID, ftr_weight(:,:)   , start=(/1, 1/), count=(/nfeatures, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lat'
    iret = nf90_inq_varid(ncid, "year_weight"   , VarId)
    iret = nf90_get_var(ncid, VarID, year_weight(:,:)  , start=(/1, 1/), count=(/nyears, nplumes   /))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading year_weight'
    iret = nf90_inq_varid(ncid, "ann_cycle"     , VarId)
    iret = nf90_get_var(ncid, VarID, ann_cycle(:,:,:)  , start=(/1, 1, 1/), count=(/nfeatures, ntimes, nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading ann_cycle'

    iret = nf90_close(ncid)

    sp_initialized = .TRUE.
    RETURN
  END SUBROUTINE sp_setup

  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight.
  !
  SUBROUTINE set_time_weight(year_fr)
    !
    ! ----------
    !
    REAL(KIND=JPRB), INTENT(IN) ::  &
         year_fr           !< Fractional Year (1850.0-2100.99)

    INTEGER(KIND=JPIM) ::  &
         idate(8)         ,& !< integer array of system clock date yyyy, mm, dd
         iyear            ,& !< Integer year values between 1 and 156 (1850-2100) 
         iweek            ,& !< Integer index (between 1 and ntimes); for ntimes = 52 this corresponds to weeks (roughly)
         iplume              ! plume number

    ! Choice of the reference year for the CDNC scale factors.
    ! For the moment, we agreed to set it to the pre-industrial control year 1850, 
    ! but in principle other choices are also possible
    ! (see discussion issue #221 on the EC-Earth 3 portal).
    INTEGER(KIND=JPIM), PARAMETER:: iyear_ref = 1   !< Integer year value for 1850
    !INTEGER(KIND=JPIM), PARAMETER:: iyear_ref = 156 !< Integer year value for 2005
                                          !< Note that 2005 is also the reference year in MAC-SP, 
                                          !< for which all year_weight factors are equal to 1.
    !
    ! ----------
    !
    iyear = FLOOR(year_fr) - 1850+1
    iweek = FLOOR((year_fr-FLOOR(year_fr)) * ntimes) + 1

    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
    DO iplume = 1, nplumes
      time_weight(1, iplume) = year_weight(iyear, iplume) * ann_cycle(1, iweek, iplume)
      time_weight(2, iplume) = year_weight(iyear, iplume) * ann_cycle(2, iweek, iplume)
      time_weight_ref(1, iplume) = year_weight(iyear_ref, iplume) * ann_cycle(1, iweek, iplume)
      time_weight_ref(2, iplume) = year_weight(iyear_ref, iplume) * ann_cycle(2, iweek, iplume)
      time_weight_bg(1, iplume) = ann_cycle(1, iweek, iplume)
      time_weight_bg(2, iplume) = ann_cycle(2, iweek, iplume)
    END DO
    !
    ! check if code is being used in beta version and stop if it is after June 1, 2016
    !
    !CALL DATE_AND_TIME(VALUES = idate)
    !IF (idate(1) >= 2016 .AND. idate(2) >= 10) THEN
    !  PRINT '(A102)', 'System clock says thde date is after June 2016 at which time this beta release should be depreciated.'
    !  PRINT '(A85)' , 'Terminating; please check with developers for a finalized CMIP6 version of the code.'
    !  STOP
    !END IF
  
    RETURN
  END SUBROUTINE set_time_weight
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical properties based on the
  ! the simple plume fit to the MPI Aerosol Climatology (Version 2).  It sums over nplumes to provide a profile of aerosol
  ! optical properties on a host models vertical grid.
  !
  SUBROUTINE sp_aop_profile (KLEV, KIDIA, KFDIA, KLON, LAMDA, LON, LAT, YEAR_FR, PGEOH, CDNC_FACTOR, AOD_PROF, SSA_PROF, ASY_PROF, ZTAOD)
    !
    ! ----------
    !
    USE YOMCT3   , ONLY : NSTEP
    INTEGER(KIND=JPIM), INTENT(IN)  :: &
         KLEV,                         & !< number of levels
         KLON                            !< number of points
    INTEGER(KIND=JPIM), INTENT(IN)  :: KIDIA
    INTEGER(KIND=JPIM), INTENT(IN)  :: KFDIA

    REAL(KIND=JPRB), INTENT(IN)     :: LAMDA                         !< wavelength
    REAL(KIND=JPRB), INTENT(IN)     :: YEAR_FR                       !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
    REAL(KIND=JPRB), INTENT(IN)     :: LON(KLON)                     !< longitude 
    REAL(KIND=JPRB), INTENT(IN)     :: LAT(KLON)                     !< latitude
    REAL(KIND=JPRB), INTENT(IN)     :: PGEOH (KLON, KLEV+1)          !< geopotential height of the model at half levels above sea-level (m)

    REAL(KIND=JPRB), INTENT(OUT)          :: &
         CDNC_FACTOR(KLON)      , & !< multiplicative change factor for CDNC (1.0 implies no change)
                                    !< defined with respect to the reference year
         AOD_PROF(KLON, KLEV)    , & !< profile of aerosol optical depth
         SSA_PROF(KLON, KLEV)    , & !< profile of single scattering albedo
         ASY_PROF(KLON, KLEV)        !< profile of asymmetry parameter

    REAL(KIND=JPRB), INTENT(OUT), OPTIONAL  :: ZTAOD(KLON)    

    INTEGER(KIND=JPIM)           :: iplume, JL, JK

    REAL(KIND=JPRB)              ::  &
         CDNC_FACTOR_MAC(KLON)   , & !< multiplicative change factor for CDNC (1.0 implies no change)
                                     !< defined with respect to the pre-industrial background
                                     !< (in the original MACv2-SP_v1 code, 
                                     !<  this variable is confusingly called dNovrN)
         CDNC_FACTOR_REF(KLON)   , & !< multiplicative change factor for CDNC (1.0 implies no change)
                                     !< for the reference year 
                                     !< defined with respect to the pre-industrial background
         eta(KLON, KLEV+1),         & !< normalized height (by 15 km)
!        z_beta(KLON, KLEV),        & !< profile for scaling column optical depth
         z (KLON, KLEV+1),          & !< geopotential height of the model at half levels above sea-level (m)
         dz(KLON, KLEV),            & !< level thickness (difference between half levels)
         prof(KLON, KLEV),          & !< scaled profile (by beta function)
         beta_sum(KLON),           & !< vertical sum of beta function
         ssa,                      & !< aerosol single-scattering albedo
         asy,                      & !< aerosol asymmetry factor
         cw_an(KLON),              & !< column weight for simple plume (anthropogenic) aod at 550 nm
         cw_bg(KLON),              & !< column weight for fine-mode background aod at 550 nm
         cw_an_ref(KLON),          & !< column weight for simple plume (anthropogenic) aod at 550 nm
                                     !< for reference year
         caod_sp(KLON),            & !< column simple plume (anthropogenic) aod at 550 nm
         caod_bg(KLON),            & !< column fine-mode background aod at 550 nm
         caod_sp_ref(KLON),        & !< column simple plume (anthropogenic) aod at 550 nm
                                     !< for reference year
         a_plume1,                 & !< gaussian longitude factor for feature 1
         a_plume2,                 & !< gaussian longitude factor for feature 2
         b_plume1,                 & !< gaussian latitude factor for feature 1
         b_plume2,                 & !< gaussian latitude factor for feature 2
         delta_lat,                & !< latitude offset
         delta_lon,                & !< longitude offset
         lon1,                     & !< rotated longitude for feature 1
         lat1,                     & !< rotated latitude for feature 2
         lon2,                     & !< rotated longitude for feature 1
         lat2,                     & !< rotated latitude for feature 2
         f1,                       & !< contribution from feature 1
         f2,                       & !< contribution from feature 2
         f3,                       & !< contribution from feature 1 in natural background of Twomey effect
         f4,                       & !< contribution from feature 2 in natural background of Twomey effect
         f5,                       & !< contribution from feature 1 for reference year
         f6,                       & !< contribution from feature 2 for reference year
         aod_550,                  & !< aerosol optical depth at 550nm
         aod_lmd,                  & !< aerosol optical depth at input wavelength
         lfactor,                  & !< factor to compute wavelength dependence of optical properties
         GINV,                     &   !< inverse of gravity constant
         BTOT(KLON)  
    !
    ! ----------
    !
    ! initialize input data (by calling setup at first instance)
    !

    IF (.NOT.sp_initialized) CALL sp_setup
    !
    ! get time weights
    !
    CALL set_time_weight(YEAR_FR)
    !
    ! Compute z from model geopotential height
    !PGEOH (JL, KLEV+1)=POROG(JL) 


    GINV = 0.1019716_JPRB  
    DO JK = 1, KLEV+1
      DO JL = KIDIA, KFDIA
       z(JL, JK)= PGEOH (JL, KLEV-JK+2)*GINV               ! Devide by gravity constant and reverse vertical order of z
       eta(JL, JK) = MAX(0.0_JPRB, MIN(1.0_JPRB, z(JL, JK)/15000._JPRB))  ! original scaling of MACv2 is up to 15km (to check with developers)
      END DO
    END DO

    ! initialize variables, including output
    !
    DO JK = 1, KLEV
      DO JL = KIDIA, KFDIA
        AOD_PROF(JL, JK) = 0.0_JPRB
        SSA_PROF(JL, JK) = 0.0_JPRB
        ASY_PROF(JL, JK) = 0.0_JPRB
!        z_beta(JL, JK)   = MERGE(1.0, 0.0, z(JL, JK) >= oro(JL))
!        when using the IFS geopotential height z-beta is always equal to 1.0 as z in this case is terrain following
!        dz(JL, JK)= eta(JL, JK+1)- eta(JL, JK)
        dz(JL, JK)= z(JL, JK+1)- z(JL, JK)

      END DO
    END DO
    DO JL = KIDIA, KFDIA
      CDNC_FACTOR_MAC(JL) = 1.0_JPRB
      caod_sp(JL)  = 0.00_JPRB
      caod_bg(JL)  = 0.02_JPRB

      CDNC_FACTOR_REF(JL) = 1.0_JPRB
      caod_sp_ref(JL) = 0.00_JPRB
      CDNC_FACTOR(JL) = 1.0_JPRB
      
      IF(PRESENT(ZTAOD)) THEN
          ZTAOD(JL) = 0.0_JPRB
      ENDIF
    END DO
    !
    ! sum contribution from plumes to construct composite profiles of aerosol optical properties
    !
    DO iplume = 1, nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      DO JL = KIDIA, KFDIA
        beta_sum(JL) = 0._JPRB
      END DO
      DO JK = 1, KLEV
        DO JL = KIDIA, KFDIA
          prof(JL, JK)   = (eta(JL, JK)**(beta_a(iplume)-1.) * (1.-eta(JL, JK))**(beta_b(iplume)-1.))*dz(JL, JK)
          beta_sum(JL) = beta_sum(JL) + prof(JL, JK)
        END DO
      END DO
      DO JK = 1, KLEV
        DO JL = KIDIA, KFDIA
          prof(JL, JK)   = prof(JL, JK)/ beta_sum(JL)
        END DO
      END DO
      !
      ! calculate plume weights
      !
      DO JL = KIDIA, KFDIA
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat   = lat(JL) - plume_lat(iplume)
        delta_lon   = lon(JL) - plume_lon(iplume)
        delta_lon   = MERGE ( delta_lon-SIGN(360._JPRB,delta_lon), delta_lon, ABS(delta_lon) > 180._JPRB )
        a_plume1  = 0.5_JPRB / (MERGE(sig_lon_E(1, iplume), sig_lon_W(1, iplume), delta_lon > 0)**2)
        b_plume1  = 0.5_JPRB / (MERGE(sig_lat_E(1, iplume), sig_lat_W(1, iplume), delta_lon > 0)**2)
        a_plume2  = 0.5_JPRB / (MERGE(sig_lon_E(2, iplume), sig_lon_W(2, iplume), delta_lon > 0)**2)
        b_plume2  = 0.5_JPRB / (MERGE(sig_lat_E(2, iplume), sig_lat_W(2, iplume), delta_lon > 0)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        !
        lon1 =   COS(theta(1, iplume))*(delta_lon) + SIN(theta(1, iplume))*(delta_lat)
        lat1 = - SIN(theta(1, iplume))*(delta_lon) + COS(theta(1, iplume))*(delta_lat)
        lon2 =   COS(theta(2, iplume))*(delta_lon) + SIN(theta(2, iplume))*(delta_lat)
        lat2 = - SIN(theta(2, iplume))*(delta_lon) + COS(theta(2, iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode background aerosol (cw_bg)
        !
        f1 = time_weight(1, iplume) * ftr_weight(1, iplume) * EXP(-1.* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f2 = time_weight(2, iplume) * ftr_weight(2, iplume) * EXP(-1.* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))


        ! In the beta version of MAC-SP the same scale factors f1 and f2 
        ! are used both for cw_an and cw_bg.
        ! Since the natural background should not depend on the year index, 
        ! this is probably a bug.
        ! I have reported the issue to Stephanie Fielder, 
        ! who confirmed that "the implementation is not elegant and
        ! we are working on it."
        ! The issue will probably be fixed in a next release.
        ! Until then, I believe the correct implementation should be
        ! to define a new set of scale factors, time_weight_bg, 
        ! which are independent of the year, 
        ! and recalculate f1 and f2 as follows:
        f3 = time_weight_bg(1, iplume) * ftr_weight(1, iplume) * EXP(-1.* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f4 = time_weight_bg(2, iplume) * ftr_weight(2, iplume) * EXP(-1.* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))

        cw_an(JL) = f1*aod_spmx(iplume) + f2*aod_spmx(iplume)
        cw_bg(JL) = f3*aod_fmbg(iplume) + f4*aod_fmbg(iplume)

        ! When using the standard IFS calculation of CDNC, 
        ! the anthropogenic part also needs to be calculated
        ! for the reference year:
        f5 = time_weight_ref(1, iplume) * ftr_weight(1, iplume) * EXP(-1.* (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2))))
        f6 = time_weight_ref(2, iplume) * ftr_weight(2, iplume) * EXP(-1.* (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2))))
        
        cw_an_ref(JL) = f5*aod_spmx(iplume) + f6*aod_spmx(iplume)
      END DO
        !
        ! calculate wavelength-dependent scattering properties
        !
        lfactor   = MIN(1.0_JPRB, 700.0_JPRB/LAMDA)
        ssa = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
        asy =  asy550(iplume) * SQRT(lfactor)
 !     END DO
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the angstrom parameter.
    !
      lfactor = EXP(-angstrom(iplume) * LOG(LAMDA/550.0_JPRB))
      DO JK = 1, KLEV
        DO JL = KIDIA, KFDIA
          aod_550          = prof(JL, JK)     * cw_an(JL)
          aod_lmd          = aod_550          * lfactor
          caod_sp(JL)    = caod_sp(JL)    + prof(JL, JK) * cw_an(JL)
          caod_bg(JL)    = caod_bg(JL)    + prof(JL, JK) * cw_bg(JL)
          caod_sp_ref(JL) = caod_sp_ref(JL) + prof(JL, JK) * cw_an_ref(JL)
          AOD_PROF(JL, JK) = AOD_PROF(JL, JK) + aod_lmd
          SSA_PROF(JL, JK) = SSA_PROF(JL, JK) + aod_lmd*ssa
          ASY_PROF(JL, JK) = ASY_PROF(JL, JK) + aod_lmd*ssa*asy
        END DO
      END DO

! DIAGNOSTIC OUTPUTS FOR DEBUGGING
!   ZTAOD(:)=0.
!   BTOT(:)=0.
!     DO JK = 1, KLEV
!      DO JL = KIDIA, KFDIA
!       IF (NSTEP == 4 .and. LON(JL) <= 22.5 .and. LON(JL) > 22.4 .and. LAT(JL) <= 45.7 .and. LAT(JL) >= 45.6 .and. LAMDA == 530. .and. iplume == 1 ) THEN
!          ZTAOD(JL)=ZTAOD(JL)+AOD_PROF(JL, JK)
!          BTOT(JL)=BTOT(JL)+PROF(JL, JK)
!         print*,  JK, z(JL, JK), ' ZTAOD= ', ZTAOD(JL), ' BTOT= ',BTOT(JL)
!       END IF
!      END DO
!     END DO
! END DIAGNOSTIC OUTPUTS FOR DEBUGGING
    END DO ! iplume

    !
    ! complete optical depth weighting
    !
    DO JK = 1, KLEV
      DO JL = KIDIA, KFDIA
        ASY_PROF(JL, JK) = MERGE(ASY_PROF(JL, JK)/SSA_PROF(JL, JK), 0.0_JPRB, SSA_PROF(JL, JK) > TINY(1._JPRB))
        SSA_PROF(JL, JK) = MERGE(SSA_PROF(JL, JK)/AOD_PROF(JL, JK), 1.0_JPRB, AOD_PROF(JL, JK) > TINY(1._JPRB))
        AOD_PROF(JL, JK) = MAX(AOD_PROF(JL, JK), 0.0_JPRB)
      END DO
    END DO
    !
    ! calculate effective radius normalization (divisor) factor
    !
    DO JL = KIDIA, KFDIA
      CDNC_FACTOR_MAC(JL) = LOG((1000.0_JPRB * (caod_sp(JL) + caod_bg(JL))) + 3.0_JPRB)/LOG((1000.0_JPRB*caod_bg(JL)) + 3.0_JPRB)
      CDNC_FACTOR_REF(JL) = LOG((1000.0_JPRB * (caod_sp_ref(JL) + caod_bg(JL))) + 3.0_JPRB)/LOG((1000.0_JPRB*caod_bg(JL)) + 3.0_JPRB)

      ! In the standard IFS code CDNC is calculated 
      ! based on a parameterization which uses CCN values, 
      ! which can be assumed to be representative for some reference year. 
      ! When this parameterization is used, 
      ! an additional scale factor should be included as follows:
      CDNC_FACTOR(JL) = CDNC_FACTOR_MAC(JL) / CDNC_FACTOR_REF(JL)

      ! When pre-industrial CDNC would be calculated directly, 
      ! based on pre-industrial aerosol fields, 
      ! the reference year should be set to 1850, 
      ! and CDNC_FACTOR_REF will be equal to 1.

      ! Calculate total aod at 550nm
      IF(PRESENT(ZTAOD)) THEN
          IF (LAMDA == 550._JPRB) THEN
            ZTAOD(JL) = SUM(AOD_PROF(JL, :))
          ENDIF
      ENDIF
    END DO

    RETURN
  END SUBROUTINE sp_aop_profile

END MODULE AER_MACv2SP_MOD
