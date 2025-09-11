MODULE YOMORB

  USE PARKIND1, ONLY: JPIS, JPIT, JPIM, JPIB, JPRB, JPRD
  USE YOMHOOK, ONLY: LHOOK, DR_HOOK, JPHOOK
  USE YOMRIP0, ONLY: NINDAT, NSSSSS
  USE YOMCST, ONLY: RDAY, RPI
  USE YOMLUN, ONLY: NULOUT, NULNAM

  IMPLICIT NONE

  PRIVATE

  LOGICAL :: LCORBMD

  CHARACTER(LEN=20)  :: ORBMODE ! orb_mode
  INTEGER(KIND=JPIB) :: ORBIY   ! orb_iyear

  REAL(KIND=JPRB) :: ORBECCEN   ! orbital eccentricity
  REAL(KIND=JPRB) :: ORBOBLIQ   ! obliquity in degrees
  REAL(KIND=JPRB) :: ORBMVELP   ! moving vernal equinox long

  REAL(KIND=JPRB) :: ORBORLIQR  ! Earths obliquity in radians
  REAL(KIND=JPRB) :: ORBLAMBM0  ! Mean long of perihelion at the vernal equinox (radians)
  REAL(KIND=JPRB) :: ORBMVELPP  ! moving vernal equinox longitude of perihelion plus pi (radians)
  REAL(KIND=JPRB) :: ORBCALDAY
  REAL(KIND=JPRB) :: ORBCALDAYM

  INTEGER(KIND=JPIB) :: ORBCYEAR  !last year mark

  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_ECCEN_MIN = 0.0_JPRB ! min value for eccen
  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_ECCEN_MAX = 0.1_JPRB ! max value for eccen
  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_OBLIQ_MIN = -90.0_JPRB ! min value for obliq
  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_OBLIQ_MAX = +90.0_JPRB ! max value for obliq
  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_MVELP_MIN = 0.0_JPRB ! min value for mvelp
  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_MVELP_MAX = 360.0_JPRB ! max value for mvelp

  REAL(KIND=JPRB), PARAMETER :: SHR_ORB_UNDEF_REAL = 1.e36_JPRB  !undefined real
  INTEGER(KIND=JPIB), PARAMETER :: SHR_ORB_UNDEF_INT = 2000000000   !undefined int

  ! Public interface of YOMORB
  PUBLIC :: LCORBMD
  PUBLIC :: ORBCALDAY
  PUBLIC :: ORBCALDAYM
  PUBLIC :: SHR_ORB_UNDEF_REAL
  PUBLIC :: SHR_ORB_UNDEF_INT
  PUBLIC :: SUORB
  PUBLIC :: SHR_ORB_PARAMS
  PUBLIC :: SHR_ORB_DECL
  PUBLIC :: SHR_ORB_PRINT
  PUBLIC :: SHR_ORB_UPD

#include "abor1.intfb.h"

CONTAINS

!====================================================================
  SUBROUTINE SHR_ORB_PARAMS(iyear_AD, eccen, obliq, mvelp, obliqr, lambm0, mvelpp)

    !----------------------------- Arguments ------------------------------------
    integer(KIND=JPIB), intent(in)    :: iyear_AD  ! Year to calculate orbit for
    real(KIND=JPRB), intent(inout) :: eccen     ! orbital eccentricity
    real(KIND=JPRB), intent(inout) :: obliq     ! obliquity in degrees
    real(KIND=JPRB), intent(inout) :: mvelp     ! moving vernal equinox long
    real(KIND=JPRB), intent(out)   :: obliqr    ! Earths obliquity in rad
    real(KIND=JPRB), intent(out)   :: lambm0    ! Mean long of perihelion at vernal equinox (radians)
    real(KIND=JPRB), intent(out)   :: mvelpp    ! moving vernal equinox long of perihelion plus pi (rad)
    !------------------------------ Parameters ----------------------------------
    integer(KIND=JPIT), parameter :: poblen = 47 ! # of elements in series wrt obliquity
    integer(KIND=JPIT), parameter :: pecclen = 19 ! # of elements in series wrt eccentricity
    integer(KIND=JPIT), parameter :: pmvelen = 78 ! # of elements in series wrt vernal equinox
    real(KIND=JPRB), parameter :: psecdeg = 1.0_JPRB/3600.0_JPRB ! arc sec to deg conversion
    real(KIND=JPRB) :: degrad             ! degree to radian conversion factor
    real(KIND=JPRB) :: yb4_1950AD         ! number of years before 1950 AD

    character(len=*), parameter :: subname = '(shr_orb_params)'

    ! Cosine series data for computation of obliquity: amplitude (arc seconds),
    ! rate (arc seconds/year), phase (degrees).

    real(KIND=JPRB), parameter :: obamp(poblen) =  (/  & ! amplitudes for obliquity cos series
        -2462.2214466_JPRB, -857.3232075_JPRB, -629.3231835_JPRB,  &
         -414.2804924_JPRB, -311.7632587_JPRB,  308.9408604_JPRB,  &
         -162.5533601_JPRB, -116.1077911_JPRB,  101.1189923_JPRB,  &
          -67.6856209_JPRB,   24.9079067_JPRB,   22.5811241_JPRB,  &
          -21.1648355_JPRB,  -15.6549876_JPRB,   15.3936813_JPRB,  &
           14.6660938_JPRB,  -11.7273029_JPRB,   10.2742696_JPRB,  &
            6.4914588_JPRB,    5.8539148_JPRB,   -5.4872205_JPRB,  &
           -5.4290191_JPRB,    5.1609570_JPRB,    5.0786314_JPRB,  &
           -4.0735782_JPRB,    3.7227167_JPRB,    3.3971932_JPRB,  &
           -2.8347004_JPRB,   -2.6550721_JPRB,   -2.5717867_JPRB,  &
           -2.4712188_JPRB,    2.4625410_JPRB,    2.2464112_JPRB,  &
           -2.0755511_JPRB,   -1.9713669_JPRB,   -1.8813061_JPRB,  &
           -1.8468785_JPRB,    1.8186742_JPRB,    1.7601888_JPRB,  &
           -1.5428851_JPRB,    1.4738838_JPRB,   -1.4593669_JPRB,  &
            1.4192259_JPRB,   -1.1818980_JPRB,    1.1756474_JPRB,  &
           -1.1316126_JPRB,    1.0896928_JPRB  &
    /)

    real(KIND=JPRB), parameter :: obrate(poblen) = (/  & ! rates for obliquity cosine series
           31.609974_JPRB, 32.620504_JPRB, 24.172203_JPRB,  &
           31.983787_JPRB, 44.828336_JPRB, 30.973257_JPRB,  &
           43.668246_JPRB, 32.246691_JPRB, 30.599444_JPRB,  &
           42.681324_JPRB, 43.836462_JPRB, 47.439436_JPRB,  &
           63.219948_JPRB, 64.230478_JPRB,  1.010530_JPRB,  &
            7.437771_JPRB, 55.782177_JPRB,  0.373813_JPRB,  &
           13.218362_JPRB, 62.583231_JPRB, 63.593761_JPRB,  &
           76.438310_JPRB, 45.815258_JPRB,  8.448301_JPRB,  &
           56.792707_JPRB, 49.747842_JPRB, 12.058272_JPRB,  &
           75.278220_JPRB, 65.241008_JPRB, 64.604291_JPRB,  &
            1.647247_JPRB,  7.811584_JPRB, 12.207832_JPRB,  &
           63.856665_JPRB, 56.155990_JPRB, 77.448840_JPRB,  &
            6.801054_JPRB, 62.209418_JPRB, 20.656133_JPRB,  &
           48.344406_JPRB, 55.145460_JPRB, 69.000539_JPRB,  &
           11.071350_JPRB, 74.291298_JPRB, 11.047742_JPRB,  &
            0.636717_JPRB, 12.844549_JPRB  &
    /)

    real(KIND=JPRB), parameter :: obphas(poblen) = (/  & ! phases for obliquity cosine series
          251.9025_JPRB, 280.8325_JPRB, 128.3057_JPRB,  &
          292.7252_JPRB,  15.3747_JPRB, 263.7951_JPRB,  &
          308.4258_JPRB, 240.0099_JPRB, 222.9725_JPRB,  &
          268.7809_JPRB, 316.7998_JPRB, 319.6024_JPRB,  &
          143.8050_JPRB, 172.7351_JPRB,  28.9300_JPRB,  &
          123.5968_JPRB,  20.2082_JPRB,  40.8226_JPRB,  &
          123.4722_JPRB, 155.6977_JPRB, 184.6277_JPRB,  &
          267.2772_JPRB,  55.0196_JPRB, 152.5268_JPRB,  &
           49.1382_JPRB, 204.6609_JPRB,  56.5233_JPRB,  &
          200.3284_JPRB, 201.6651_JPRB, 213.5577_JPRB,  &
           17.0374_JPRB, 164.4194_JPRB,  94.5422_JPRB,  &
          131.9124_JPRB,  61.0309_JPRB, 296.2073_JPRB,  &
          135.4894_JPRB, 114.8750_JPRB, 247.0691_JPRB,  &
          256.6114_JPRB,  32.1008_JPRB, 143.6804_JPRB,  &
           16.8784_JPRB, 160.6835_JPRB,  27.5932_JPRB,  &
          348.1074_JPRB,  82.6496_JPRB  &
    /)

    ! Cosine/sine series data for computation of eccentricity and fixed vernal
    ! equinox longitude of perihelion (fvelp): amplitude,
    ! rate (arc seconds/year), phase (degrees).

    real(KIND=JPRB), parameter :: ecamp(pecclen) = (/  & ! ampl for eccen/fvelp cos/sin series
            0.01860798_JPRB,  0.01627522_JPRB, -0.01300660_JPRB,  &
            0.00988829_JPRB, -0.00336700_JPRB,  0.00333077_JPRB,  &
           -0.00235400_JPRB,  0.00140015_JPRB,  0.00100700_JPRB,  &
            0.00085700_JPRB,  0.00064990_JPRB,  0.00059900_JPRB,  &
            0.00037800_JPRB, -0.00033700_JPRB,  0.00027600_JPRB,  &
            0.00018200_JPRB, -0.00017400_JPRB, -0.00012400_JPRB,  &
            0.00001250_JPRB  &
    /)

    real(KIND=JPRB), parameter :: ecrate(pecclen) = (/  & ! rates for eccen/fvelp cos/sin series
            4.2072050_JPRB,  7.3460910_JPRB, 17.8572630_JPRB,  &
           17.2205460_JPRB, 16.8467330_JPRB,  5.1990790_JPRB,  &
           18.2310760_JPRB, 26.2167580_JPRB,  6.3591690_JPRB,  &
           16.2100160_JPRB,  3.0651810_JPRB, 16.5838290_JPRB,  &
           18.4939800_JPRB,  6.1909530_JPRB, 18.8677930_JPRB,  &
           17.4255670_JPRB,  6.1860010_JPRB, 18.4174410_JPRB,  &
            0.6678630_JPRB  &
    /)

    real(KIND=JPRB), parameter :: ecphas(pecclen) = (/  & ! phases for eccen/fvelp cos/sin series
           28.620089_JPRB, 193.788772_JPRB, 308.307024_JPRB,  &
          320.199637_JPRB, 279.376984_JPRB,  87.195000_JPRB,  &
          349.129677_JPRB, 128.443387_JPRB, 154.143880_JPRB,  &
          291.269597_JPRB, 114.860583_JPRB, 332.092251_JPRB,  &
          296.414411_JPRB, 145.769910_JPRB, 337.237063_JPRB,  &
          152.092288_JPRB, 126.839891_JPRB, 210.667199_JPRB,  &
           72.108838_JPRB  &
    /)

    ! Sine series data for computation of moving vernal equinox longitude of
    ! perihelion: amplitude (arc seconds), rate (arc sec/year), phase (degrees).

    real(KIND=JPRB), parameter :: mvamp(pmvelen) = (/  & ! amplitudes for mvelp sine series
         7391.0225890_JPRB, 2555.1526947_JPRB, 2022.7629188_JPRB,  &
        -1973.6517951_JPRB, 1240.2321818_JPRB,  953.8679112_JPRB,  &
         -931.7537108_JPRB,  872.3795383_JPRB,  606.3544732_JPRB,  &
         -496.0274038_JPRB,  456.9608039_JPRB,  346.9462320_JPRB,  &
         -305.8412902_JPRB,  249.6173246_JPRB, -199.1027200_JPRB,  &
          191.0560889_JPRB, -175.2936572_JPRB,  165.9068833_JPRB,  &
          161.1285917_JPRB,  139.7878093_JPRB, -133.5228399_JPRB,  &
          117.0673811_JPRB,  104.6907281_JPRB,   95.3227476_JPRB,  &
           86.7824524_JPRB,   86.0857729_JPRB,   70.5893698_JPRB,  &
          -69.9719343_JPRB,  -62.5817473_JPRB,   61.5450059_JPRB,  &
          -57.9364011_JPRB,   57.1899832_JPRB,  -57.0236109_JPRB,  &
          -54.2119253_JPRB,   53.2834147_JPRB,   52.1223575_JPRB,  &
          -49.0059908_JPRB,  -48.3118757_JPRB,  -45.4191685_JPRB,  &
          -42.2357920_JPRB,  -34.7971099_JPRB,   34.4623613_JPRB,  &
          -33.8356643_JPRB,   33.6689362_JPRB,  -31.2521586_JPRB,  &
          -30.8798701_JPRB,   28.4640769_JPRB,  -27.1960802_JPRB,  &
           27.0860736_JPRB,  -26.3437456_JPRB,   24.7253740_JPRB,  &
           24.6732126_JPRB,   24.4272733_JPRB,   24.0127327_JPRB,  &
           21.7150294_JPRB,  -21.5375347_JPRB,   18.1148363_JPRB,  &
          -16.9603104_JPRB,  -16.1765215_JPRB,   15.5567653_JPRB,  &
           15.4846529_JPRB,   15.2150632_JPRB,   14.5047426_JPRB,  &
          -14.3873316_JPRB,   13.1351419_JPRB,   12.8776311_JPRB,  &
           11.9867234_JPRB,   11.9385578_JPRB,   11.7030822_JPRB,  &
           11.6018181_JPRB,  -11.2617293_JPRB,  -10.4664199_JPRB,  &
           10.4333970_JPRB,  -10.2377466_JPRB,   10.1934446_JPRB,  &
          -10.1280191_JPRB,   10.0289441_JPRB,  -10.0034259_JPRB  &
    /)

    real(KIND=JPRB), parameter :: mvrate(pmvelen) = (/  & ! rates for mvelp sine series
           31.609974_JPRB, 32.620504_JPRB, 24.172203_JPRB,  &
            0.636717_JPRB, 31.983787_JPRB,  3.138886_JPRB,  &
           30.973257_JPRB, 44.828336_JPRB,  0.991874_JPRB,  &
            0.373813_JPRB, 43.668246_JPRB, 32.246691_JPRB,  &
           30.599444_JPRB,  2.147012_JPRB, 10.511172_JPRB,  &
           42.681324_JPRB, 13.650058_JPRB,  0.986922_JPRB,  &
            9.874455_JPRB, 13.013341_JPRB,  0.262904_JPRB,  &
            0.004952_JPRB,  1.142024_JPRB, 63.219948_JPRB,  &
            0.205021_JPRB,  2.151964_JPRB, 64.230478_JPRB,  &
           43.836462_JPRB, 47.439436_JPRB,  1.384343_JPRB,  &
            7.437771_JPRB, 18.829299_JPRB,  9.500642_JPRB,  &
            0.431696_JPRB,  1.160090_JPRB, 55.782177_JPRB,  &
           12.639528_JPRB,  1.155138_JPRB,  0.168216_JPRB,  &
            1.647247_JPRB, 10.884985_JPRB,  5.610937_JPRB,  &
           12.658184_JPRB,  1.010530_JPRB,  1.983748_JPRB,  &
           14.023871_JPRB,  0.560178_JPRB,  1.273434_JPRB,  &
           12.021467_JPRB, 62.583231_JPRB, 63.593761_JPRB,  &
           76.438310_JPRB,  4.280910_JPRB, 13.218362_JPRB,  &
           17.818769_JPRB,  8.359495_JPRB, 56.792707_JPRB,  &
            8.448301_JPRB,  1.978796_JPRB,  8.863925_JPRB,  &
            0.186365_JPRB,  8.996212_JPRB,  6.771027_JPRB,  &
           45.815258_JPRB, 12.002811_JPRB, 75.278220_JPRB,  &
           65.241008_JPRB, 18.870667_JPRB, 22.009553_JPRB,  &
           64.604291_JPRB, 11.498094_JPRB,  0.578834_JPRB,  &
            9.237738_JPRB, 49.747842_JPRB,  2.147012_JPRB,  &
            1.196895_JPRB,  2.133898_JPRB,  0.173168_JPRB  &
    /)

    real(KIND=JPRB), parameter :: mvphas(pmvelen) = (/  & ! phases for mvelp sine series
          251.9025_JPRB, 280.8325_JPRB, 128.3057_JPRB,  &
          348.1074_JPRB, 292.7252_JPRB, 165.1686_JPRB,  &
          263.7951_JPRB,  15.3747_JPRB,  58.5749_JPRB,  &
           40.8226_JPRB, 308.4258_JPRB, 240.0099_JPRB,  &
          222.9725_JPRB, 106.5937_JPRB, 114.5182_JPRB,  &
          268.7809_JPRB, 279.6869_JPRB,  39.6448_JPRB,  &
          126.4108_JPRB, 291.5795_JPRB, 307.2848_JPRB,  &
           18.9300_JPRB, 273.7596_JPRB, 143.8050_JPRB,  &
          191.8927_JPRB, 125.5237_JPRB, 172.7351_JPRB,  &
          316.7998_JPRB, 319.6024_JPRB,  69.7526_JPRB,  &
          123.5968_JPRB, 217.6432_JPRB,  85.5882_JPRB,  &
          156.2147_JPRB,  66.9489_JPRB,  20.2082_JPRB,  &
          250.7568_JPRB,  48.0188_JPRB,   8.3739_JPRB,  &
           17.0374_JPRB, 155.3409_JPRB,  94.1709_JPRB,  &
          221.1120_JPRB,  28.9300_JPRB, 117.1498_JPRB,  &
          320.5095_JPRB, 262.3602_JPRB, 336.2148_JPRB,  &
          233.0046_JPRB, 155.6977_JPRB, 184.6277_JPRB,  &
          267.2772_JPRB,  78.9281_JPRB, 123.4722_JPRB,  &
          188.7132_JPRB, 180.1364_JPRB,  49.1382_JPRB,  &
          152.5268_JPRB,  98.2198_JPRB,  97.4808_JPRB,  &
          221.5376_JPRB, 168.2438_JPRB, 161.1199_JPRB,  &
           55.0196_JPRB, 262.6495_JPRB, 200.3284_JPRB,  &
          201.6651_JPRB, 294.6547_JPRB,  99.8233_JPRB,  &
          213.5577_JPRB, 154.1631_JPRB, 232.7153_JPRB,  &
          138.3034_JPRB, 204.6609_JPRB, 106.5938_JPRB,  &
          250.4676_JPRB, 332.3345_JPRB,  27.3039_JPRB  &
    /)

    !---------------------------Local variables----------------------------------
    integer(KIND=JPIT) :: i
    real(KIND=JPRB) :: obsum   ! Obliquity series summation
    real(KIND=JPRB) :: cossum  ! Cos series summation for eccentricity/fvelp
    real(KIND=JPRB) :: sinsum  ! Sin series summation for eccentricity/fvelp
    real(KIND=JPRB) :: fvelp   ! Fixed vernal equinox long of perihelion
    real(KIND=JPRB) :: mvsum   ! mvelp series summation
    real(KIND=JPRB) :: beta    ! Intermediate argument for lambm0
    real(KIND=JPRB) :: years   ! Years to time of interest ( pos <=> future)
    real(KIND=JPRB) :: eccen2  ! eccentricity squared
    real(KIND=JPRB) :: eccen3  ! eccentricity cubed

    !-------------------------- Formats -----------------------------------------
    character(len=*), parameter :: F00 = "('(shr_orb_params) ',4a)"
    character(len=*), parameter :: F01 = "('(shr_orb_params) ',a,i9)"
    character(len=*), parameter :: F02 = "('(shr_orb_params) ',a,f6.3)"
    character(len=*), parameter :: F03 = "('(shr_orb_params) ',a,es14.6)"

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_ORB_PARAMS', 0, ZHOOK_HANDLE)
    !----------------------------------------------------------------------------
    ! radinp and algorithms below will need a degree to radian conversion factor
    ! Check for flag to use input orbit parameters

    degrad = RPI/180.0_JPRB

    IF (iyear_AD == SHR_ORB_UNDEF_INT) THEN

      ! Check input obliq, eccen, and mvelp to ensure reasonable

      if (obliq == SHR_ORB_UNDEF_REAL) then
        write (NULOUT, F00) trim(subname)//' Have to specify orbital parameters:'
        write (NULOUT, F00) 'Either set: iyear_AD, OR [obliq, eccen, and mvelp]:'
        write (NULOUT, F00) 'iyear_AD is the year to simulate orbit for (ie. 1950): '
        write (NULOUT, F00) 'obliq, eccen, mvelp specify the orbit directly:'
        write (NULOUT, F00) 'The AMIP II settings (for a 1995 orbit) are: '
        write (NULOUT, F00) ' obliq =  23.4441'
        write (NULOUT, F00) ' eccen =   0.016715'
        write (NULOUT, F00) ' mvelp = 102.7'
        CALL ABOR1(subname//' ERROR: unreasonable obliq')
      end if
      if ((obliq < SHR_ORB_OBLIQ_MIN) .or. (obliq > SHR_ORB_OBLIQ_MAX)) then
        write (NULOUT, F03) 'Input obliquity unreasonable: ', obliq
        call ABOR1(subname//' ERROR: unreasonable obliq')
      end if
      if ((eccen < SHR_ORB_ECCEN_MIN) .or. (eccen > SHR_ORB_ECCEN_MAX)) then
        write (NULOUT, F03) 'Input eccentricity unreasonable: ', eccen
        call ABOR1(subname//' ERROR: unreasonable eccen')
      end if
      if ((mvelp < SHR_ORB_MVELP_MIN) .or. (mvelp > SHR_ORB_MVELP_MAX)) then
        write (NULOUT, F03) 'Input mvelp unreasonable: ', mvelp
        call ABOR1(subname//' ERROR: unreasonable mvelp')
      end if
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

    ELSE  ! Otherwise calculate based on years before present

      yb4_1950AD = 1950.0_JPRB - real(iyear_AD, JPRB)
      if (abs(yb4_1950AD) .gt. 1000000.0_JPRB) then
        write (NULOUT, F00) 'orbit only valid for years+-1000000'
        write (NULOUT, F00) 'Relative to 1950 AD'
        write (NULOUT, F03) '# of years before 1950: ', yb4_1950AD
        write (NULOUT, F01) 'Year to simulate was  : ', iyear_AD
        call ABOR1(subname//' ERROR: unreasonable year')
      end if

      years = -yb4_1950AD

      obsum = 0.0_JPRB
      do i = 1, poblen
        obsum = obsum  &
            + obamp(i)*psecdeg*cos(  &
                (obrate(i)*psecdeg*years + obphas(i))*degrad  &
            )
      end do
      obliq = 23.320556_JPRB + obsum

      cossum = 0.0_JPRB
      do i = 1, pecclen
        cossum = cossum + ecamp(i)*cos((ecrate(i)*psecdeg*years + ecphas(i))*degrad)
      end do

      sinsum = 0.0_JPRB
      do i = 1, pecclen
        sinsum = sinsum + ecamp(i)*sin((ecrate(i)*psecdeg*years + ecphas(i))*degrad)
      end do

      ! Use summations to calculate eccentricity

      eccen2 = cossum*cossum + sinsum*sinsum
      eccen = sqrt(eccen2)
      eccen3 = eccen2*eccen

      ! A series of cases for fvelp, which is in radians.

      if (abs(cossum) .le. 1.0E-8_JPRB) then
        if (sinsum .eq. 0.0_JPRB) then
          fvelp = 0.0_JPRB
        else if (sinsum .lt. 0.0_JPRB) then
          fvelp = 1.5_JPRB*RPI
        else if (sinsum .gt. 0.0_JPRB) then
          fvelp = .5_JPRB*RPI
        end if
      else if (cossum .lt. 0.0_JPRB) then
        fvelp = atan(sinsum/cossum) + RPI
      else if (cossum .gt. 0.0_JPRB) then
        if (sinsum .lt. 0.0_JPRB) then
          fvelp = atan(sinsum/cossum) + 2.0_JPRB*RPI
        else
          fvelp = atan(sinsum/cossum)
        end if
      end if

      mvsum = 0.0_JPRB
      do i = 1, pmvelen
        mvsum = mvsum  &
            + mvamp(i)*psecdeg*sin(  &
                (mvrate(i)*psecdeg*years + mvphas(i))*degrad  &
            )
      end do
      mvelp = fvelp/degrad + 50.439273_JPRB*psecdeg*years + 3.392506_JPRB + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp .lt. 0.0_JPRB)
        mvelp = mvelp + 360.0_JPRB
      end do
      do while (mvelp .ge. 360.0_JPRB)
        mvelp = mvelp - 360.0_JPRB
      end do

    END IF  ! end of test on whether to calculate or use input orbital params

    ! Orbit needs the obliquity in radians

    obliqr = obliq*degrad

    mvelpp = (mvelp + 180.0_JPRB)*degrad

    ! Set up an argument used several times in lambm0 calculation ahead.

    beta = sqrt(1.0_JPRB - eccen2)

    ! The mean longitude at the vernal equinox (lambda m nought in Berger
    ! 1978; in radians) is calculated from the following formula given in
    ! Berger 1978.  At the vernal equinox the true longitude (lambda in Berger
    ! 1978) is 0.

    lambm0 =  &
        2.0_JPRB*(  &
            (0.5_JPRB*eccen + 0.125_JPRB*eccen3)*(1.0_JPRB + beta)*sin(mvelpp)  &
            - 0.25_JPRB*eccen2*(0.5_JPRB + beta)*sin(2.0_JPRB*mvelpp)  &
            + 0.125_JPRB*eccen3*(1.0_JPRB/3.0_JPRB + beta)*sin(3.0_JPRB*mvelpp)  &
        )

    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_ORB_PARAMS', 1, ZHOOK_HANDLE)

  END SUBROUTINE SHR_ORB_PARAMS

!====================================================================
  SUBROUTINE SHR_ORB_DECL(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf)

    !------------------------------Arguments--------------------------------
    real(KIND=JPRB), intent(in)  :: calday ! Calendar day, including fraction
    real(KIND=JPRB), intent(in)  :: eccen  ! Eccentricity
    real(KIND=JPRB), intent(in)  :: obliqr ! Earths obliquity in radians
    real(KIND=JPRB), intent(in)  :: lambm0 ! Mean long of perihelion at the vernal equinox (radians)
    real(KIND=JPRB), intent(in)  :: mvelpp ! moving vernal equinox longitude of perihelion plus pi (radians)
    real(KIND=JPRB), intent(out) :: delta  ! Solar declination angle in rad
    real(KIND=JPRB), intent(out) :: eccf   ! Earth-sun distance factor (ie. (1/r)**2)

    !---------------------------Local variables-----------------------------
    real(KIND=JPRB), parameter :: dayspy = 365.25_JPRB  ! days per year
    real(KIND=JPRB), parameter :: ve = 80.5_JPRB   ! Calday of vernal equinox
    ! assumes Jan 1 = calday 1

    real(KIND=JPRB) ::   lambm  ! Lambda m, mean long of perihelion (rad)
    real(KIND=JPRB) ::   lmm    ! Intermediate argument involving lambm
    real(KIND=JPRB) ::   lamb   ! Lambda, the earths long of perihelion
    real(KIND=JPRB) ::   invrho ! Inverse normalized sun/earth distance
    real(KIND=JPRB) ::   sinl   ! Sine of lmm

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_ORB_DECL', 0, ZHOOK_HANDLE)

    lambm = lambm0 + (calday - ve)*2.0_JPRB*RPI/dayspy
    lmm = lambm - mvelpp

    ! The earths true longitude, in radians, is then found from
    ! the formula in Berger 1978:

    sinl = sin(lmm)
    lamb = lambm + eccen*(  &
            2.0_JPRB*sinl + eccen*(  &
                1.25_JPRB*sin(2.0_JPRB*lmm) + eccen*(  &
                    (13.0_JPRB/12.0_JPRB)*sin(3.0_JPRB*lmm) - 0.25_JPRB*sinl  &
                )  &
            )  &
        )

    ! Using the obliquity, eccentricity, moving vernal equinox longitude of
    ! perihelion (plus), and earths true longitude, the declination (delta)
    ! and the normalized earth/sun distance (rho in Berger 1978; actually inverse
    ! rho will be used), and thus the eccentricity factor (eccf), can be
    ! calculated from formulas given in Berger 1978.

    invrho = (1.0_JPRB + eccen*cos(lamb - mvelpp))/(1.0_JPRB - eccen*eccen)

    ! Set solar declination and eccentricity factor

    delta = asin(sin(obliqr)*sin(lamb))
    eccf = invrho*invrho
    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_ORB_DECL', 1, ZHOOK_HANDLE)
    return

  END SUBROUTINE SHR_ORB_DECL

!====================================================================
  SUBROUTINE SHR_ORB_PRINT(iyear_AD, eccen, obliq, mvelp)

    !---------------------------Arguments----------------------------------------
    integer(KIND=JPIB), intent(in) :: iyear_AD ! requested Year (AD)
    real(KIND=JPRB), intent(in) :: eccen    ! eccentricity (unitless, typically 0 to 0.1)
    real(KIND=JPRB), intent(in) :: obliq    ! obliquity (-90 to +90 degrees, typically 22-26)
    real(KIND=JPRB), intent(in) :: mvelp    ! moving vernal equinox at perhel (0 to 360 degrees)
    !-------------------------- Formats -----------------------------------------
    character(len=*), parameter :: F00 = "('(shr_orb_print) ',4a)"
    character(len=*), parameter :: F01 = "('(shr_orb_print) ',a,i9.4)"
    character(len=*), parameter :: F02 = "('(shr_orb_print) ',a,f6.3)"
    character(len=*), parameter :: F03 = "('(shr_orb_print) ',a,es14.6)"
    !----------------------------------------------------------------------------
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_ORB_PRINT', 0, ZHOOK_HANDLE)

    if (iyear_AD .ne. SHR_ORB_UNDEF_INT) then
      if (iyear_AD > 0) then
        write (NULOUT, F01) 'Orbital parameters calculated for year: AD ', iyear_AD
      else
        write (NULOUT, F01) 'Orbital parameters calculated for year: BC ', iyear_AD
      end if
    else if (obliq /= SHR_ORB_UNDEF_REAL) then
      write (NULOUT, F03) 'Orbital parameters: '
      write (NULOUT, F03) 'Obliquity (degree):              ', obliq
      write (NULOUT, F03) 'Eccentricity (unitless):         ', eccen
      write (NULOUT, F03) 'Long. of moving Perhelion (deg): ', mvelp
    else
      write (NULOUT, F03) 'Orbit parameters not set!'
    end if

    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_ORB_PRINT', 1, ZHOOK_HANDLE)
  END SUBROUTINE SHR_ORB_PRINT

!====================================================================
  SUBROUTINE SHR_CAL_DATE2YMD(date, year, month, day)

    integer(KIND=JPIM), intent(in)  :: date             ! coded-date (yyyymmdd)
    integer(KIND=JPIM), intent(out) :: year, month, day   ! calendar year,month,day

    integer(KIND=JPIB) :: tdate   ! temporary date
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_CAL_DATE2YMD', 0, ZHOOK_HANDLE)

    tdate = abs(date)
    year = int(tdate/10000)
    if (date < 0) year = -year
    month = int(mod(tdate, 10000)/100)
    day = mod(tdate, 100)

    IF (LHOOK) CALL DR_HOOK('YOMORB:SHR_CAL_DATE2YMD', 1, ZHOOK_HANDLE)
  END SUBROUTINE SHR_CAL_DATE2YMD

!====================================================================
  SUBROUTINE GET_CURR_CALDAY(YEAR, MONTH, DAY, SEC, CALDAY)

    INTEGER(KIND=JPIM), INTENT(IN)  :: YEAR, MONTH, DAY, SEC
    REAL(KIND=JPRB), INTENT(OUT) :: CALDAY

    integer(KIND=JPIS) A(12), I
    DATA A(1:12) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('YOMORB:GET_CURR_CALDAY', 0, ZHOOK_HANDLE)

    CALDAY = 0.0_JPRB
    IF (  &
      MOD(YEAR, 4) == 0 .AND. MOD(YEAR, 400) /= 100  &
      .AND. MOD(YEAR, 400) /= 200 .AND. MOD(YEAR, 400) /= 300  &
    ) THEN
      A(2) = 29
    ENDIF

    DO I = 1, MONTH - 1
      CALDAY = CALDAY + A(I)
    END DO

    CALDAY = CALDAY + DAY + SEC/86400.0_JPRB

    IF (CALDAY > 366.0_JPRB .AND. CALDAY <= 367.0_JPRB) THEN
      CALDAY = CALDAY - 1.0_JPRB
    END IF

    IF (CALDAY < 1.0_JPRB .OR. CALDAY > 366.0_JPRB) THEN
      CALL ABOR1('YOMORB:GET_CURR_CALDAY: ABOR1 CALLED')
    END IF

    IF (LHOOK) CALL DR_HOOK('YOMORB:GET_CURR_CALDAY', 1, ZHOOK_HANDLE)
  END SUBROUTINE GET_CURR_CALDAY

!====================================================================
  SUBROUTINE SHR_ORB_UPD(ZRSTATI, CALDAY, DELTA, ECCF)

    REAL(KIND=JPRD), INTENT(IN)  :: ZRSTATI
    REAL(KIND=JPRB), INTENT(OUT) :: CALDAY
    REAL(KIND=JPRB), INTENT(OUT) :: DELTA
    REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: ECCF

    INTEGER(KIND=JPIB) :: NYEAR
    INTEGER(KIND=JPIM) :: JROF, IJ0, IM0, IA0, IDINCR, OYY, OMM, ODD, OLMOIS(12)
    INTEGER(KIND=JPIM) :: ISEC
    REAL(KIND=JPRB) :: ZECCF

    CHARACTER(LEN=*), PARAMETER :: F00 = "(' (ORB_PARAMS) ',I6,6F12.8)"

    IJ0 = MOD(NINDAT, 100)
    IM0 = MOD((NINDAT - IJ0)/100, 100)
    IA0 = NINDAT/10000

    IDINCR = (NSSSSS + NINT(ZRSTATI, JPIB))/NINT(RDAY)
    CALL UPDCAL(IJ0, IM0, IA0, IDINCR, ODD, OMM, OYY, OLMOIS, -1)

    IF (TRIM(ORBMODE) == "variable_year") THEN
      NYEAR = ORBIY + (OYY - IA0)
      IF (ORBCYEAR /= NYEAR) THEN
        ORBCYEAR = NYEAR
        CALL SHR_ORB_PARAMS(NYEAR, ORBECCEN, ORBOBLIQ, ORBMVELP, ORBORLIQR, ORBLAMBM0, ORBMVELPP)
        WRITE (NULOUT, F00) NYEAR, ORBECCEN, ORBOBLIQ, ORBMVELP, ORBORLIQR, ORBLAMBM0, ORBMVELPP
      END IF
    END IF

    ISEC = MOD(NSSSSS + NINT(ZRSTATI, JPIB), NINT(RDAY))
    CALL GET_CURR_CALDAY(OYY, OMM, ODD, ISEC, CALDAY)

    CALL SHR_ORB_DECL(CALDAY, ORBECCEN, ORBMVELPP, ORBLAMBM0, ORBORLIQR, DELTA, ZECCF)
    IF (PRESENT(ECCF)) ECCF = ZECCF

  END SUBROUTINE SHR_ORB_UPD

!====================================================================
  SUBROUTINE SUORB(KULOUT)

    IMPLICIT NONE

    INTEGER(KIND=JPIM), INTENT(IN)    :: KULOUT

    INTEGER(KIND=JPIM) :: IJ0, IM0, IA0, IDINCR, OYY, OMM, ODD, OLMOIS(12)
    INTEGER(KIND=JPIM) :: ISTAT
    REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "posname.intfb.h"
#include "fcttim.func.h"
#include "namorb.h"

!-------------------------------------------------
! 1. Default values.
!-------------------------------------------------

    IF (LHOOK) CALL DR_HOOK('SUORB', 0, ZHOOK_HANDLE)

! set default value
    LCORBMD = .FALSE.
    ORBMODE = "fixed_year"
    ORBIY = 1850
    ORBECCEN = 0.0
    ORBOBLIQ = 0.0
    ORBMVELP = 0.0
    ORBCALDAY = 0.0
    ORBCALDAYM = 0.0

! modify default value
    CALL POSNAME(NULNAM, 'NAMORB', ISTAT)
    SELECT CASE (ISTAT)
    CASE (0)
      READ (NULNAM, NAMORB)
    CASE (1)
      WRITE (KULOUT, '(A,L)') 'SUORB : CANNOT LOCATE NAMORB, USING LCORBMD=',LCORBMD
    CASE DEFAULT
      CALL ABOR1('POSNAM:READ ERROR IN NAMELIST FILE')
    END SELECT

!---------------------------------------------------------------
! check orbital mode, reset unused parameters, validate settings
!---------------------------------------------------------------
    IF (LCORBMD) THEN
      IF (TRIM(ORBMODE) == "fixed_year") THEN
        ORBOBLIQ = SHR_ORB_UNDEF_REAL
        ORBECCEN = SHR_ORB_UNDEF_REAL
        ORBMVELP = SHR_ORB_UNDEF_REAL
        IF (ORBIY == SHR_ORB_UNDEF_INT) THEN
          WRITE (KULOUT, *) 'ERROR: invalid settings ORBMODE =', TRIM(ORBMODE)
          WRITE (KULOUT, *) 'ERROR: fixed_year settings = ', ORBIY
          CALL ABOR1('ERROR: invalid settings for ORBMODE')
        END IF
      ELSEIF (TRIM(ORBMODE) == "variable_year") THEN
        ORBOBLIQ = SHR_ORB_UNDEF_REAL
        ORBECCEN = SHR_ORB_UNDEF_REAL
        ORBMVELP = SHR_ORB_UNDEF_REAL
        IF (ORBIY == SHR_ORB_UNDEF_INT) THEN
          WRITE (KULOUT, *) 'ERROR: invalid settings ORBMODE =', TRIM(ORBMODE)
          WRITE (KULOUT, *) 'ERROR: variable_year settings = ', ORBIY
          CALL ABOR1('ERROR: invalid settings for ORBMODE')
        END IF
      ELSEIF (TRIM(ORBMODE) == "fixed_parameters") THEN
        !-- force orb_iyear to undef to make sure shr_orb_params works properly
        ORBIY = SHR_ORB_UNDEF_INT
        IF (ORBECCEN == SHR_ORB_UNDEF_REAL .OR. &
            ORBOBLIQ == SHR_ORB_UNDEF_REAL .OR. &
            ORBMVELP == SHR_ORB_UNDEF_REAL) THEN
          WRITE (KULOUT, *) 'ERROR: invalid settings ORBMODE =', TRIM(ORBMODE)
          WRITE (KULOUT, *) 'ERROR: ORBECCEN = ', ORBECCEN
          WRITE (KULOUT, *) 'ERROR: ORBOBLIQ = ', ORBOBLIQ
          WRITE (KULOUT, *) 'ERROR: ORBMVELP = ', ORBMVELP
          CALL ABOR1('ERROR: invalid settings for ORBMODE')
        END IF
      ELSE
        CALL ABOR1('ERROR: invalid ORBMODE')
      END IF

      CALL SHR_ORB_PARAMS(ORBIY, ORBECCEN, ORBOBLIQ, ORBMVELP, ORBORLIQR, ORBLAMBM0, ORBMVELPP)

      WRITE (KULOUT, *) "========== MODULE YOMORB ==========="
      WRITE (KULOUT, *) "LCORBMD=", LCORBMD
      WRITE (KULOUT, *) "ORBMODE=", ORBMODE, " ORBCYEAR=", ORBCYEAR
      WRITE (KULOUT, *) "ORBIY=", ORBIY
      WRITE (KULOUT, *) "ORBECCEN=", ORBECCEN, " ORBOBLIQ=", ORBOBLIQ
      WRITE (KULOUT, *) "ORBMVELP=", ORBMVELP, " ORBORLIQR=", ORBORLIQR
      WRITE (KULOUT, *) "ORBLAMBM0=", ORBLAMBM0, " ORBMVELPP=", ORBMVELPP
      WRITE (KULOUT, *) "====================================="
    ELSE
      WRITE (KULOUT, *) "========== MODULE YOMORB ==========="
      WRITE (KULOUT, *) "LCORBMD=", LCORBMD
      WRITE (KULOUT, *) "====================================="
    END IF
    IF (LHOOK) CALL DR_HOOK('SUORB', 1, ZHOOK_HANDLE)

  END SUBROUTINE SUORB

!====================================================================

END MODULE YOMORB
