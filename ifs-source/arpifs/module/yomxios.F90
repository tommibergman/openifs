MODULE yomxios

USE PARKIND1, ONLY : JPIM, JPRM, JPRB
USE YOMPHYDS, ONLY : JPVEXTRL
USE xios

IMPLICIT NONE

SAVE

! XIOS context
TYPE(xios_context)  :: context_handle
CHARACTER(LEN=*), PARAMETER :: model_name="OpenIFS"
CHARACTER(LEN=*), PARAMETER :: ifs_context="oifs"

! Variables definition
CHARACTER(LEN=9),  PARAMETER :: lopt_send_var_name="lopt_send"
CHARACTER(LEN=17), PARAMETER :: lsingle_prec_send_var_name="lsingle_prec_send"
CHARACTER(LEN=5),  PARAMETER :: nfitp_var_name="nfitp"
CHARACTER(LEN=5),  PARAMETER :: nfitt_var_name="nfitt"
CHARACTER(LEN=5),  PARAMETER :: nfitv_var_name="nfitv"
CHARACTER(LEN=6),  PARAMETER :: nfpcli_var_name="nfpcli"
CHARACTER(LEN=4),  PARAMETER :: lfpq_var_name="lfpq"
CHARACTER(LEN=8),  PARAMETER :: ltracefp_var_name="ltracefp"
CHARACTER(LEN=7),  PARAMETER :: rfpcorr_var_name="rfpcorr"
CHARACTER(LEN=6),  PARAMETER :: nfrpos_var_name="nfrpos"
CHARACTER(LEN=6),  PARAMETER :: nfrhis_var_name="nfrhis"

! Calendar management
TYPE(xios_duration) :: time_step
TYPE(xios_date)     :: time_origin
TYPE(xios_date)     :: start_date
INTEGER(KIND=JPIM)  :: nstep_from_origin
TYPE(xios_duration) :: duration_from_origin

! Axes definition
CHARACTER(LEN=12), PARAMETER :: model_axis_name="model_levels"
CHARACTER(LEN=15), PARAMETER :: pressure_axis_name="pressure_levels"
CHARACTER(LEN=12), PARAMETER :: theta_axis_name="theta_levels"
CHARACTER(LEN=9),  PARAMETER :: pv_axis_name="pv_levels"
!CHARACTER(LEN=13), PARAMETER :: height_axis_name="height_levels"

! Domains definition
CHARACTER(LEN=16), PARAMETER :: gaussian_domain_name="reduced_gaussian"
CHARACTER(LEN=7),  PARAMETER :: regular_domain_name="regular"

! Surface shortnames and grib codes
INTEGER(KIND=JPIM), PARAMETER :: NSFCFLD=123
CHARACTER (LEN=16),  PARAMETER :: CSFCFLD(NSFCFLD)=&
    &(/ 'sro             ', 'ssro            ', 'ci              ', 'asn             ', 'rsn             ', &
    &   'sstk            ', 'istl1           ', 'istl2           ', 'istl3           ', 'istl4           ', &
    &   'swvl1           ', 'swvl2           ', 'swvl3           ', 'swvl4           ', 'es              ', &
    &   'smlt            ', '10fg            ', 'lspf            ', 'uvb             ', 'par             ', &
    &   'cape            ', 'tclw            ', 'tciw            ', 'mx2t6           ', 'mn2t6           ', &
    &   '10fg6           ', 'emis            ', 'vite            ', 'sz              ', 'sp              ', &
    &   'tcw             ', 'tcwv            ', 'stl1            ', 'stl2            ', 'stl3            ', &
    &   'stl4            ', 'sd              ', 'lsp             ', 'cp              ', 'sf              ', &
    &   'bld             ', 'sshf            ', 'slhf            ', 'chnk            ', 'msl             ', &
    &   'lnsp            ', 'blh             ', 'tcc             ', '10u             ', '10v             ', &
    &   '2t              ', '2d              ', 'ssrd            ', 'lsm             ', 'sr              ', &
    &   'al              ', 'strd            ', 'ssr             ', 'str             ', 'tsr             ', &
    &   'ttr             ', 'ewss            ', 'nsss            ', 'e               ', 'lcc             ', &
    &   'mcc             ', 'hcc             ', 'sund            ', 'lgws            ', 'mgws            ', &
    &   'gwd             ', 'src             ', 'mx2t            ', 'mn2t            ', 'ro              ', &
    &   'tco3            ', 'tsrc            ', 'ttrc            ', 'ssrc            ', 'strc            ', &
    &   'tisr            ', 'vimd            ', 'tp              ', 'iews            ', 'inss            ', &
    &   'ishf            ', 'ie              ', 'lsrh            ', 'skt             ', 'tsn             ', &
    &   'fal             ', 'fsr             ', 'flsr            ', 'cin             ', 'lmlt            ', &
    &   'lmld            ', 'lblt            ', 'ltlt            ', 'lshf            ', 'lict            ', &
    &   'licd            ', 'deg0l           ', 'tcrw            ', 'tcsw            ', '100u            ', &
    &   '100v            ', 'kx              ', 'totalx          ', &
    &   'tvl             ', 'cvl             ', 'lai_lv          ', &
    &   'tvh             ', 'cvh             ', 'lai_hv          ', 'macv2sp_taod550 ', 'aod550          ', &
    &   'aod550ss        ', 'aod550du        ', 'aod550om        ', 'aod550bc        ', 'aod550su        ', &
    &   'macv2sp_cdncf   ', &
    &   'lwcs            ' /)
INTEGER(KIND=JPIM), PARAMETER :: IGRBSFCFLD(NSFCFLD)=&
    &(/                  8,                  9,                 31,                 32,                 33, &
    &                   34,                 35,                 36,                 37,                 38, &
    &                   39,                 40,                 41,                 42,                 44, &
    &                   45,                 49,                 50,                 57,                 58, &
    &                   59,                 78,                 79,                121,                122, &
    &                  123,                124,                125,                129,                134, &
    &                  136,                137,                139,                170,                183, &
    &                  236,                141,                142,                143,                144, &
    &                  145,                146,                147,                148,                151, &
    &                  152,                159,                164,                165,                166, &
    &                  167,                168,                169,                172,                173, &
    &                  174,                175,                176,                177,                178, &
    &                  179,                180,                181,                182,                186, &
    &                  187,                188,                189,                195,                196, &
    &                  197,                198,                201,                202,                205, &
    &                  206,                208,                209,                210,                211, &
    &                  212,                213,                228,                229,                230, &
    &                  231,                232,                234,                235,                238, &
    &                  243,                244,                245,             228001,             228008, &
    &               228009,             228010,             228011,             228012,             228013, &
    &               228014,             228024,             228089,             228090,             228246, &
    &               228247,             260121,             260123, &
    &                   29,                 27,                 66, &
    &                   30,                 28,                 67,             215089,             210207, &
    &               210208,             210209,             210210,             210211,             210212, &   ! AOD for diff. aerosol types
    &               210243, &
    &               228038 /)

! 3D field shortnames and grib codes
INTEGER(KIND=JPIM), PARAMETER :: N3DFLD=20
CHARACTER (LEN=16), PARAMETER :: C3DFLD(N3DFLD)=&
    &(/ 'pt              ', 'mont            ', 'pres            ', 'pv              ', 'crwc            ', &
    &   'cswc            ', 'etadot          ', 'z               ', 't               ', 'u               ', &
    &   'v               ', 'q               ', 'w               ', 'vo              ', 'd               ', &
    &   'r               ', 'o3              ', 'clwc            ', 'ciwc            ', 'cc              '/)
INTEGER(KIND=JPIM), PARAMETER :: IGRB3DFLD(N3DFLD)=&
    &(/                  3,                 53,                 54,                 60,                 75, &
    &                   76,                 77,                129,                130,                131, &
    &                  132,                133,                135,                138,                155, &
    &                  157,                203,                246,                247,                248/)

! 3D PEXTRA fields ids and grib codes
INTEGER(KIND=JPIM), PARAMETER :: NPEXTRAFLD=24
CHARACTER (LEN=16), PARAMETER :: CPEXTRAFLD(NPEXTRAFLD)=&
    &(/ 'pextra_91       ', 'pextra_92       ', 'pextra_93       ', 'pextra_94       ', 'pextra_95       ', &
    &   'pextra_96       ', 'pextra_97       ', 'pextra_98       ', 'pextra_99       ', 'pextra_100      ', &
    &   'pextra_101      ', 'pextra_102      ', 'pextra_103      ', 'pextra_104      ', 'pextra_105      ', &
    &   'pextra_106      ', 'pextra_107      ', 'pextra_108      ', 'pextra_109      ', 'pextra_110      ', &
    &   'pextra_111      ', 'pextra_112      ', 'pextra_113      ', 'pextra_114      '/)
INTEGER(KIND=JPIM), PARAMETER :: IGRBPEXTRAFLD(NPEXTRAFLD)=&
    &(/                 91,                 92,                 93,                 94,                 95, &
    &                   96,                 97,                 98,                 99,                100, &
    &                  101,                102,                103,                104,                105, &
    &                  106,                107,                108,                109,                110, &
    &                  111,                112,                113,                114/)

! PEXTRA setup variables
LOGICAL :: LBUD23=.FALSE.
INTEGER(KIND=JPIM) :: NVEXTR=0
INTEGER(KIND=JPIM) :: NCEXTR=0
INTEGER(KIND=JPIM) :: NVEXTRAGB(JPVEXTRL)

! Buffers & variables to delay I/O and/or enable single precision sends to XIOS
LOGICAL :: LOPT_SEND=.FALSE.
LOGICAL :: LSINGLE_PREC_SEND=.FALSE.
! Surface fields
INTEGER(KIND=JPIM) :: NSFCFLDBUF=0
CHARACTER (LEN=16) :: CSFCFLDBUF(NSFCFLD)
REAL(KIND=JPRM), ALLOCATABLE :: SFCFLDBUF_SP(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: SFCFLDBUF_DP(:,:)
! Model levels
INTEGER(KIND=JPIM) :: NMLFLDBUF=0
CHARACTER (LEN=16) :: CMLFLDBUF(N3DFLD)
REAL(KIND=JPRM), ALLOCATABLE :: MLFLDBUF_SP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: MLFLDBUF_DP(:,:,:)
! Pressure levels
INTEGER(KIND=JPIM) :: NPLFLDBUF=0
CHARACTER (LEN=16) :: CPLFLDBUF(N3DFLD)
REAL(KIND=JPRM), ALLOCATABLE :: PLFLDBUF_SP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: PLFLDBUF_DP(:,:,:)
! Theta levels
INTEGER(KIND=JPIM) :: NTLFLDBUF=0
CHARACTER (LEN=16) :: CTLFLDBUF(N3DFLD)
REAL(KIND=JPRM), ALLOCATABLE :: TLFLDBUF_SP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: TLFLDBUF_DP(:,:,:)
! PV levels
INTEGER(KIND=JPIM) :: NVLFLDBUF=0
CHARACTER (LEN=16) :: CVLFLDBUF(N3DFLD)
REAL(KIND=JPRM), ALLOCATABLE :: VLFLDBUF_SP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: VLFLDBUF_DP(:,:,:)

END MODULE yomxios
