MODULE CPLNG2_DATA_MOD

    USE MOD_OASIS

    USE PARKIND1, ONLY: JPIM, JPRB
    USE CPLNG2_TYPES_MOD, ONLY: CPLNG2_FLD_TYPE

    IMPLICIT NONE

    PRIVATE

    PUBLIC CPLNG2_FLD

    PUBLIC CPLNG2_FLD_TYPE_GRIDPOINT
    PUBLIC CPLNG2_FLD_TYPE_SPECTRAL

    PUBLIC CPLNG2_FLD_IN
    PUBLIC CPLNG2_FLD_OUT

    PUBLIC CPLNG2_ADD_FLD
    PUBLIC CPLNG2_ADD_FLD_COMPLETED
    PUBLIC CPLNG2_IDX

    TYPE(CPLNG2_FLD_TYPE), ALLOCATABLE, TARGET :: CPLNG2_FLD(:)

    INTEGER(KIND=JPIM), PARAMETER :: CPLNG2_FLD_TYPE_GRIDPOINT = 0
    INTEGER(KIND=JPIM), PARAMETER :: CPLNG2_FLD_TYPE_SPECTRAL = 1

    INTEGER, PARAMETER :: CPLNG2_FLD_IN = OASIS_IN
    INTEGER, PARAMETER :: CPLNG2_FLD_OUT = OASIS_OUT

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

        IF (.NOT. ALLOCATED(CPLNG2_FLD)) THEN
            CALL ABOR1("CPLNG2_IDX: CPLNG2_FLD not allocated upon call")
        END IF

        ! Simple loop searching for the field name
        DO i = 1, SIZE(CPLNG2_FLD)
            IF (TRIM(fld_name) == TRIM(CPLNG2_FLD(i)%name)) EXIT
        END DO

        IF (i > SIZE(CPLNG2_FLD)) THEN
            CALL ABOR1("CPLNG2_IDX: Field name not found: "//TRIM(FLD_NAME))
        END IF

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

        USE MOD_OASIS
#ifdef WITH_XIOS
        USE XIOS
#endif

        ! arguments
        TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
        ! locals
        TYPE(CPLNG2_FLD_TYPE), ALLOCATABLE :: fld_tmp(:)

        ! variables needed in OASIS calls
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

            ! Reallocate CPLNG2_FLD with the correct lengths
            ! (cutting extra lengths added by chunk allocation in CPLNG2_ADD_FLD)
            ALLOCATE (FLD_TMP(NUM_FIELDS), SOURCE=CPLNG2_FLD(1:NUM_FIELDS))
            CALL MOVE_ALLOC(FLD_TMP, CPLNG2_FLD)

            ! -------------------------------------------------------------------------
            ! * SET UP OASIS PARTITION
            ! -------------------------------------------------------------------------

            ! Grid point fields
            IF (ANY(CPLNG2_FLD(:)%TYPE == CPLNG2_FLD_TYPE_GRIDPOINT)) THEN
                ! Set up logical array leg_start that indicats where OASIS segments start
                ! Method: Use the NGLOBALINDEX array from YOMMP that contains a
                !         local-to-global mapping of the grid points. If two consecutive
                !         indices in this array differ by any other value than one, it's the
                !         start of a new segment.
                leg_start(1) = .TRUE.
                leg_start(2:NGPTOT) = NGLOBALINDEX(2:NGPTOT) - NGLOBALINDEX(1:NGPTOT - 1) /= 1

                num_segs = COUNT(leg_start) ! Counts number of segments

                ! For shape and meaning of oas_part_spec, see IG_PARAL(:) in OASIS documentation
                ALLOCATE (oas_part_spec(2 + 2*num_segs))

                oas_part_spec(1) = 3 ! Value 3 indicates an orange partition, which is
                ! an ensemble of segments of the global domain
                oas_part_spec(2) = num_segs

                ! Compute segment offsets and extents and store them in oas_part_spec for later
                ! use in the OASIS_DEF_PARTITION call
                offset = NGLOBALINDEX(1)
                extend = 1
                ip = 3 ! Pointer into oas_part_spec array
                DO i = 2, NGPTOT
                    IF (leg_start(i)) THEN
                        ! OASIS counts offsets starting with zero, hence minus one
                        oas_part_spec(ip) = offset - 1
                        oas_part_spec(ip + 1) = extend
                        ip = ip + 2
                        offset = NGLOBALINDEX(i)
                        extend = 1
                    ELSE
                        extend = extend + 1
                    END IF
                END DO
                oas_part_spec(ip) = offset - 1
                oas_part_spec(ip + 1) = extend

                ! Define partition for OASIS
                CALL OASIS_DEF_PARTITION(oas_part_id_gp, oas_part_spec, error)
                IF (error /= OASIS_OK) THEN
                    WRITE (err_str, '(I3)') error
                    CALL ABOR1("CPLNG2_ADD_FLD_COMPLETED: Error on OASIS_DEF_PARTITION (gridpoint): "//err_str)
                END IF

                DEALLOCATE (oas_part_spec)
            END IF ! Field type == gridpoint

            ! Spectral fields
            IF (ANY(CPLNG2_FLD(:)%TYPE == CPLNG2_FLD_TYPE_SPECTRAL)) THEN
                ! Orange partition: each processor owns a number of (not following) wave numbers.
                !
                !  NSMAX        : Spectral truncation T
                !  NUMP         : Number of wavenumbers locally
                !  MYMS(1:NUMP) : Actual wave numbers   (NOTE: not 0:NUMP as in comment!)

                ! Allocate definition array
                ALLOCATE (oas_part_spec(2 + 2*NUMP))

                ! Indicator for orange partition
                oas_part_spec(1) = 3

                ! Number of segments
                oas_part_spec(2) = NUMP

                ! Global offset and size
                DO i = 1, NUMP
                    ! Global wave number
                    ispm = MYMS(i)
                    ! Offset
                    oas_part_spec(1 + i*2) = ((NSMAX + 1) + (NSMAX + 2 - ispm))*ispm
                    ! Size
                    oas_part_spec(2 + i*2) = (NSMAX + 1 - ispm)*2  ! M:T, Re/Im
                END DO

                ! Define partition for OASIS
                CALL OASIS_DEF_PARTITION(oas_part_id_sp, oas_part_spec, error)
                IF (error /= OASIS_OK) THEN
                    WRITE (err_str, '(I3)') error
                    CALL ABOR1("CPLNG2_ADD_FLD_COMPLETED: Error on OASIS_DEF_PARTITION (spectral): "//err_str)
                END IF

                DEALLOCATE (oas_part_spec)
            END IF ! Field type == spectral

            ! -------------------------------------------------------------------------
            ! * DEFINE COUPLING FIELDS FOR OASIS
            ! -------------------------------------------------------------------------
            oas_nodims = (/1, 1/)

            WRITE (NULOUT, '(/,1X,A,/,1X,A,I4)') &
                'CPLNG2: Coupling field definition.', &
                'CPLNG2: ... Number of fields is: ', SIZE(CPLNG2_FLD)

            DO i = 1, SIZE(CPLNG2_FLD)
                WRITE (NULOUT, '(1X,A,I4,A,I2,A,I2,A,A20)') &
                    'CPLNG2: .... Field ', i, &
                    '(', CPLNG2_FLD(i)%num_lvl, ' lvl,', &
                    CPLNG2_FLD(i)%num_cat, ' cat) Name: ', TRIM(CPLNG2_FLD(i)%name)

                ! Allocate for ID's of all levels/categories of the field
                ALLOCATE (CPLNG2_FLD(i)%id(CPLNG2_FLD(i)%num_lvl, CPLNG2_FLD(i)%num_cat))

                ! Set OASIS partition id and data size according to field type
                ! (gp/spectral). Allocate data member.
                SELECT CASE (CPLNG2_FLD(i)%TYPE)

                CASE (CPLNG2_FLD_TYPE_GRIDPOINT)
                    oas_part_id = oas_part_id_gp
                    oas_actual_shape = (/1, NGPTOT/)
                    ALLOCATE (CPLNG2_FLD(i)%d(NGPTOT, CPLNG2_FLD(i)%num_lvl, CPLNG2_FLD(i)%num_cat))

                CASE (CPLNG2_FLD_TYPE_SPECTRAL)
                    oas_part_id = oas_part_id_sp
                    oas_actual_shape = (/1, NSPEC2/)
                    ALLOCATE (CPLNG2_FLD(i)%d(NSPEC2, CPLNG2_FLD(i)%num_lvl, CPLNG2_FLD(i)%num_cat))

                CASE DEFAULT
                    CALL ABOR1("CPLNG2_ADD_FLD_COMPLETED: Wrong field type for "//TRIM(CPLNG2_FLD(i)%name))

                END SELECT

                ! Initialise to a conspicuous value
                CPLNG2_FLD(i)%d = -HUGE(0.0_JPRB)

                DO ICAT = 1, CPLNG2_FLD(i)%num_cat
                    DO ILVL = 1, CPLNG2_FLD(i)%num_lvl

                        fld_name = TRIM(CPLNG2_FLD(i)%name)
                        IF (CPLNG2_FLD(i)%num_cat > 1) WRITE (fld_name, '(A,".C",I3.3)') TRIM(fld_name), ICAT
                        IF (CPLNG2_FLD(i)%num_lvl > 1) WRITE (fld_name, '(A,".L",I3.3)') TRIM(fld_name), ILVL

                        CALL OASIS_DEF_VAR(CPLNG2_FLD(i)%id(ILVL, ICAT), &
                                           fld_name, &
                                           oas_part_id, &
                                           oas_nodims, &
                                           CPLNG2_FLD(i)%inout, &
                                           oas_actual_shape, &
                                           OASIS_REAL, &
                                           error)
                        WRITE (NULOUT, '(1X,A,I4,A,I3)') 'CPLNG2: .... Field ', i, ' OASIS_DEF_VAR returns ', error

                        IF (error /= OASIS_OK) THEN
                            WRITE (err_str, '(I3)') error
                            CALL ABOR1("CPLNG2_ADD_FLD_COMPLETED: Error in OASIS_DEF_VAR: "//err_str// &
                                       " ("//TRIM(fld_name)//")")
                        END IF
                    END DO
                END DO
            END DO

            ! -------------------------------------------------------------------------
            ! * FINALISE OASIS DEFINITION PHASE
            ! -------------------------------------------------------------------------
#ifdef WITH_XIOS
            CALL xios_oasis_enddef()
#endif
            CALL OASIS_ENDDEF(error)
            IF (error /= OASIS_OK) THEN
                WRITE (err_str, '(I3)') error
                CALL ABOR1("CPLNG2_ADD_FLD_COMPLETED: Error in OASIS_ENDDEF: "//err_str)
            END IF
        END ASSOCIATE
    END SUBROUTINE CPLNG2_ADD_FLD_COMPLETED

END MODULE CPLNG2_DATA_MOD
