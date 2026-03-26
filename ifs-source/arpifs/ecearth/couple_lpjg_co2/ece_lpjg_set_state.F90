SUBROUTINE ECE_LPJG_SET_STATE(KDIM,YDMODEL,PGFL)

#ifdef WITH_CPLNG2

    USE PARKIND1,   ONLY: JPIM, JPRB
    USE YOMHOOK,    ONLY: LHOOK, DR_HOOK, JPHOOK
    USE YOMPHYDER,  ONLY: DIMENSION_TYPE
    USE TYPE_MODEL, ONLY: MODEL
    USE YOMCST,     ONLY : RMD, RMCO2
    USE CPLNG2
    USE YOM_GRIB_CODES, ONLY: NGRBGHG

    IMPLICIT NONE

    ! Arguments
    TYPE(DIMENSION_TYPE),          INTENT(IN) :: KDIM
    TYPE(MODEL),                   INTENT(IN) :: YDMODEL
    REAL(KIND=JPRB),              INTENT(IN) :: PGFL(KDIM%KLON,KDIM%KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)

    ! Locals
    REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
    INTEGER(KIND=JPIM) :: IL,IE,IG,IH,JL,ICO2_INDEX,IDCO2,JT

#include "update_fields.intfb.h"

    IF (LHOOK) CALL DR_HOOK('ECE_LPJG_SET_STATE',0,ZHOOK_HANDLE)
    ASSOCIATE( &
         & YGHG=>YDMODEL%YRML_GCONF%YGFL%YGHG, & 
         & NGHG=>YDMODEL%YRML_GCONF%YGFL%NGHG)

    ! =========================================================================
    ! *** Pre-compute indices
    ! =========================================================================

    IL = KDIM%KIDIA
    IE = KDIM%KFDIA - KDIM%KIDIA
    IG = KDIM%KSTGLO - 1 + KDIM%KIDIA

    ! =========================================================================
    ! *** CO2 concentration from lowest model level to LPJG
    ! =========================================================================
    
    ! Get the index of CO2_GHG in the GFL array
    IDCO2=-999_JPIM
    DO JT=1,NGHG
      IF (YGHG(JT)%IGRBCODE == NGRBGHG(1) ) IDCO2=JT
    ENDDO
    
    ICO2_INDEX = YGHG(IDCO2)%MP

    ! Send CO2 from lowest model level (level 1) via coupling
    CPLNG2_FLD(CPLNG2_IDX('XCO2Veg'))%D(IG:IG+IE,1,1) = PGFL(IL:IL+IE,KDIM%KLEV,ICO2_INDEX) * RMD / RMCO2 * (10**6)

    END ASSOCIATE
    IF (LHOOK) CALL DR_HOOK('ECE_LPJG_SET_STATE',1,ZHOOK_HANDLE)

#endif

END SUBROUTINE ECE_LPJG_SET_STATE
