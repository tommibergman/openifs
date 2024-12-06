MODULE CPLNG2_DATA_MOD


    USE PARKIND1, ONLY: JPIM, JPRB
    USE CPLNG2_TYPES_MOD, ONLY: CPLNG2_FLD_TYPE

    IMPLICIT NONE

    PRIVATE

    PUBLIC CPLNG2_FLD

    PUBLIC CPLNG2_FLD_TYPE_GRIDPOINT
    PUBLIC CPLNG2_FLD_TYPE_SPECTRAL

    PUBLIC CPLNG2_ADD_FLD
    PUBLIC CPLNG2_ADD_FLD_COMPLETED
    PUBLIC CPLNG2_IDX

    TYPE(CPLNG2_FLD_TYPE), ALLOCATABLE, TARGET :: CPLNG2_FLD(:)

    INTEGER(KIND=JPIM), PARAMETER :: CPLNG2_FLD_TYPE_GRIDPOINT = 0
    INTEGER(KIND=JPIM), PARAMETER :: CPLNG2_FLD_TYPE_SPECTRAL = 1


    ! NUM_FIELDS is for internal bookkeeping
    INTEGER(KIND=JPIM) :: NUM_FIELDS = 0

CONTAINS

    FUNCTION CPLNG2_IDX(fld_name)
        ! arguments
        CHARACTER(LEN=*), INTENT(IN) :: fld_name
        ! result
        INTEGER(KIND=JPIM) :: CPLNG2_IDX
        ! locals
        INTEGER(KIND=JPIM) :: i

        CPLNG2_IDX = i
    END FUNCTION CPLNG2_IDX

    SUBROUTINE CPLNG2_ADD_FLD(name, TYPE, inout, stage, lvl, cat)
        ! parameters
        INTEGER, PARAMETER :: ALLOCATE_CHUNK = 20
        ! arguments
        CHARACTER(LEN=*), INTENT(IN) :: name
        INTEGER(KIND=JPIM), INTENT(IN) :: TYPE
        INTEGER(KIND=JPIM), INTENT(IN) :: inout
        INTEGER(KIND=JPIM), INTENT(IN) :: stage
        INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: lvl
        INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: cat
        ! locals
        INTEGER(KIND=JPIM) :: num_lvl, num_cat
        TYPE(CPLNG2_FLD_TYPE), ALLOCATABLE :: fld_tmp(:)

        ! Set number of levels/categories to default or optional arguments
        num_lvl = 1
        num_cat = 1

        IF (PRESENT(lvl)) num_lvl = lvl
        IF (PRESENT(cat)) num_cat = cat

        ! Make initial allocation, if necessary
        IF (.NOT. ALLOCATED(CPLNG2_FLD)) THEN
            ALLOCATE (CPLNG2_FLD(ALLOCATE_CHUNK))
        END IF

        ! Increase (by reallocation) size of CPLNG2_FLD, if necessary
        IF (NUM_FIELDS == SIZE(CPLNG2_FLD)) THEN
            ALLOCATE (fld_tmp(SIZE(CPLNG2_FLD) + ALLOCATE_CHUNK))
            fld_tmp(1:NUM_FIELDS) = CPLNG2_FLD
            CALL MOVE_ALLOC(fld_tmp, CPLNG2_FLD)
        END IF

        NUM_FIELDS = NUM_FIELDS + 1

        CPLNG2_FLD(NUM_FIELDS)%name = name
        CPLNG2_FLD(NUM_FIELDS)%TYPE = TYPE
        CPLNG2_FLD(NUM_FIELDS)%inout = inout
        CPLNG2_FLD(NUM_FIELDS)%stage = stage
        CPLNG2_FLD(NUM_FIELDS)%num_lvl = num_lvl
        CPLNG2_FLD(NUM_FIELDS)%num_cat = num_cat
    END SUBROUTINE CPLNG2_ADD_FLD

    SUBROUTINE CPLNG2_ADD_FLD_COMPLETED(YDGEOMETRY)
        USE GEOMETRY_MOD, ONLY: GEOMETRY
        USE YOMLUN, ONLY: NULOUT
        USE YOMCT0, ONLY: LXIOS


        ! arguments
        TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
        ! locals
        TYPE(CPLNG2_FLD_TYPE), ALLOCATABLE :: fld_tmp(:)

        INTEGER(KIND=JPIM) :: error
        INTEGER(KIND=JPIM), ALLOCATABLE :: oas_part_spec(:)
        INTEGER(KIND=JPIM) :: oas_part_id, oas_part_id_gp, oas_part_id_sp
        INTEGER(KIND=JPIM) :: oas_nodims(2)
        INTEGER(KIND=JPIM) :: oas_actual_shape(2)

        LOGICAL :: leg_start(YDGEOMETRY%YRGEM%NGPTOT)
        INTEGER(KIND=JPIM) :: i, ip
        INTEGER(KIND=JPIM) :: ilvl, icat
        INTEGER(KIND=JPIM) :: offset
        INTEGER(KIND=JPIM) :: extend
        INTEGER(KIND=JPIM) :: num_segs
        CHARACTER(LEN=128) :: fld_name
        CHARACTER(LEN=3) :: err_str

        INTEGER(KIND=JPIM) :: ispm

        ! Short-cuts for geometry description variables
        ASSOCIATE (NGPTOT => YDGEOMETRY%YRGEM%NGPTOT, &
                   NSMAX => YDGEOMETRY%YRDIM%NSMAX, &
                   NUMP => YDGEOMETRY%YRDIM%NUMP, &
                   NSPEC2 => YDGEOMETRY%YRDIM%NSPEC2, &
                   NGLOBALINDEX => YDGEOMETRY%YRMP%NGLOBALINDEX, &
                   MYMS => YDGEOMETRY%YRLAP%MYMS)

        END ASSOCIATE
    END SUBROUTINE CPLNG2_ADD_FLD_COMPLETED

END MODULE CPLNG2_DATA_MOD
