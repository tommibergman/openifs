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

SUBROUTINE SUOROG (YDGEOMETRY, PSPOR)

!**** *SUOROG*  - Routine to initialize the grid-point orography fields.

!     Purpose.
!     --------
!           Initialize the grid-point orography fields of the model.

!**   Interface.
!     ----------
!        *CALL* *SUOROG(.....)*

!        Explicit arguments :
!        --------------------
!          INPUT:
!            PSPOR   - spectral orography field

!        Implicit arguments :
!        --------------------
!        The grid-point fields of the model.
!        The boundary condition fields of the model.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      David Dent *ECMWF*
!      Original : 92-05-19

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      R. El Khatib 14-May-2018 move allocations to sugeometry and merge initializations with sueorog
!      J. Bernales   11-Dec-2025 Read ice sheet orography from external file (optional)
!      J. Bernales   04-Mar-2026 Read ice sheet mask from external file (optional)
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM
#ifdef WITH_ATLAS
USE ATLAS_MODULE, ONLY : ATLAS_REAL, ATLAS_FIELD, ATLAS_MESH_NODES, &
 & ATLAS_FUNCTIONSPACE_NODECOLUMNS
#endif
USE TYPE_GEOMETRY , ONLY : GEOMETRY
USE YOMDYNA_STATIC  , ONLY : LGRADGP

! Ice sheet coupling
USE ECEARTH     , ONLY : ECE_CPL_ISMM, ISMRGNPFX, NISMRGNS
USE SURFECE     , ONLY : ECE_LANDICE, SURFECE_SET_LANDICE_LOCAL, SURFECE_RESET_LANDICE
USE YOMCST      , ONLY : RG
USE YOMMP0      , ONLY : MYPROC
USE DISGRID_MOD , ONLY : DISGRID_SEND, DISGRID_RECV
USE NETCDF     , ONLY : nf90_open, NF90_NOWRITE, NF90_NOERR, nf90_strerror, &
 &                       nf90_inq_varid, nf90_get_var, nf90_close

IMPLICIT NONE

TYPE(GEOMETRY)  , INTENT(INOUT)    :: YDGEOMETRY
REAL(KIND=JPRB) , INTENT(INOUT) :: PSPOR(YDGEOMETRY%YRDIM%NSPEC2)

REAL(KIND=JPRB) , ALLOCATABLE :: ZU(:), ZSWORK(:)

REAL(KIND=JPRB) :: ZDIV      (YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZZOROG    (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGL   (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGM   (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGLL  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGMM  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZZOROGLM  (YDGEOMETRY%YRGEM%NGPTOT)

#ifdef WITH_ATLAS
INTEGER(KIND=JPIM),ALLOCATABLE :: IPA(:)
REAL(KIND=JPRB), POINTER :: ZVAR(:,:)
REAL(KIND=JPRB), POINTER :: ZGRAD(:,:,:)
LOGICAL   , POINTER :: LLIS_GHOST(:)
TYPE(ATLAS_FIELD) :: GHOSTFIELD
TYPE(ATLAS_FIELD) :: VARFIELD
TYPE(ATLAS_FIELD) :: GRADFIELD
TYPE(ATLAS_MESH_NODES) :: NODES
TYPE(ATLAS_FUNCTIONSPACE_NODECOLUMNS) :: NODE_COLUMNS
#endif
INTEGER(KIND=JPIM) :: JL,IP,INODE_SIZE
INTEGER(KIND=JPIM) :: JKGLO, IEND, IBL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Ice sheet coupling
INTEGER(KIND=JPIM) :: II
INTEGER :: NCID_IN, VARID_USURF, VARID_PLIT, IRET_NC, IRET_DUMMY
INTEGER :: IACTIVE_CELLS
REAL(KIND=JPRB) :: ZPLIT_SUM
REAL(KIND=JPRB), ALLOCATABLE :: ZGLOB_USURF(:)  ! Global GP buffer
REAL(KIND=JPRB), ALLOCATABLE :: ZGLOB_PLIT(:)   ! Global GP buffer
REAL(KIND=JPRB) :: ZLOC_USURF(YDGEOMETRY%YRGEM%NGPTOT)  ! Local GP buffer
REAL(KIND=JPRB) :: ZLOC_PLIT(YDGEOMETRY%YRGEM%NGPTOT)   ! Combined local mask (accumulated across regions)
REAL(KIND=JPRB) :: ZLOC_PLIT_RGN(YDGEOMETRY%YRGEM%NGPTOT) ! Per-region scatter target
REAL(KIND=JPRB), ALLOCATABLE :: ZLOC_PLIT_B(:,:)  ! Blocked form (NPROMA,NGPBLKS)
INTEGER(KIND=JPIM), PARAMETER :: RPRC = 1_JPIM  ! Root I/O task

#include "fv_gradient.intfb.h"
#include "reespe.intfb.h"
#include "speree.intfb.h"
#include "speuv.intfb.h"
#include "sueorog.intfb.h"

!     ------------------------------------------------------------------

!*       1.    INITIALIZE GRIDPOINT OROGRAPHY.
!              -------------------------------

IF (LHOOK) CALL DR_HOOK('SUOROG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT)

IF (LELAM) THEN

!*       1.1   PLANE GEOMETRY.

  CALL SUEOROG (YDGEOMETRY, PSPOR, ZZOROG, ZZOROGL, ZZOROGM, ZZOROGLL, ZZOROGLM, ZZOROGMM)

ELSE

!*       1.2   SPHERICAL GEOMETRY.

!*       1.2.1   ALLOCATIONS.

!*       1.2.2   INITIALIZE GRIDPOINT OROGRAPHY

  CALL SPEREE(YDGEOMETRY,1,1,PSPOR,ZZOROG)

!*       1.2.2.1   ICE SHEET FIELDS
!
! One region loop serves both features. Each region file
! ({prefix}_pism2ece.nc) can contain 'usurf' and/or 'plit'.
! ECE_CPL_ISMM reads orography (usurf). ECE_LANDICE reads
! the fractional ice sheet mask (plit, 0-1 after regridding).

  IF (ECE_CPL_ISMM .OR. ECE_LANDICE) THEN

    ! Allocate global buffers on root
    IF (MYPROC == RPRC) THEN
      IF (ECE_CPL_ISMM) ALLOCATE(ZGLOB_USURF(YDGEOMETRY%YRGEM%NGPTOTG))
      IF (ECE_LANDICE)  ALLOCATE(ZGLOB_PLIT(YDGEOMETRY%YRGEM%NGPTOTG))
    END IF
    IF (ECE_LANDICE) ZLOC_PLIT(:) = 0._JPRB

    DO II=1,NISMRGNS

      ! --- Open region file ---
      IF (MYPROC == RPRC) THEN
        IRET_NC = nf90_open(TRIM(ISMRGNPFX(II))//"_pism2ece.nc", &
         & NF90_NOWRITE, NCID_IN)
        IF (IRET_NC /= NF90_NOERR) THEN
          WRITE(*,*) '>>> SUOROG: FATAL - cannot open ', &
           & TRIM(ISMRGNPFX(II))//'_pism2ece.nc: ', &
           & TRIM(nf90_strerror(IRET_NC))
          CALL ABOR1('SUOROG: Cannot open ' &
           & //TRIM(ISMRGNPFX(II))//'_pism2ece.nc')
        END IF
      END IF

      ! --- ISMM: read orography (usurf) ---
      IF (ECE_CPL_ISMM) THEN
        ZLOC_USURF(:) = 0._JPRB
        IF (MYPROC == RPRC) THEN
          ZGLOB_USURF(:) = 0._JPRB
          IRET_NC = nf90_inq_varid(NCID_IN, 'usurf_oifs', VARID_USURF)
          IF (IRET_NC == NF90_NOERR) THEN
            IRET_NC = nf90_get_var(NCID_IN, VARID_USURF, ZGLOB_USURF, &
             & start=(/1,1/), count=(/YDGEOMETRY%YRGEM%NGPTOTG,1/))
          END IF
          IF (IRET_NC == NF90_NOERR) THEN
            write(*,*) '>>> SUOROG: read USURF from ', TRIM(ISMRGNPFX(II))
            ZGLOB_USURF(:) = ZGLOB_USURF(:) * RG
          ELSE
            WRITE(*,*) '>>> SUOROG: FATAL - cannot read USURF: ', &
             & IRET_NC, TRIM(nf90_strerror(IRET_NC))
            CALL ABOR1('SUOROG: USURF read failed from ' &
             & //TRIM(ISMRGNPFX(II))//'_pism2ece.nc')
          END IF
        END IF
        IF (MYPROC == RPRC) THEN
          CALL DISGRID_SEND(YDGEOMETRY, 1, ZGLOB_USURF, 1, ZLOC_USURF)
        ELSE
          CALL DISGRID_RECV(YDGEOMETRY, RPRC, 1, ZLOC_USURF, 1)
        END IF
        ! Patch orography where this region has ice sheet data
        ! TODO: usurf > 0 fails for ice below sea level (e.g. bedrock depression)
        WHERE (ZLOC_USURF > 0._JPRB) ZZOROG = ZLOC_USURF
      END IF

      ! --- LANDICE: read fractional ice sheet mask (plit) ---
      IF (ECE_LANDICE) THEN
        ZLOC_PLIT_RGN(:) = 0._JPRB
        IF (MYPROC == RPRC) THEN
          ZGLOB_PLIT(:) = 0._JPRB
          IRET_NC = nf90_inq_varid(NCID_IN, 'plit_oifs', VARID_PLIT)
          IF (IRET_NC == NF90_NOERR) THEN
            IRET_NC = nf90_get_var(NCID_IN, VARID_PLIT, ZGLOB_PLIT, &
             & start=(/1,1/), count=(/YDGEOMETRY%YRGEM%NGPTOTG,1/))
          END IF
          IF (IRET_NC /= NF90_NOERR) THEN
            WRITE(*,*) '>>> SUOROG: FATAL - cannot read PLIT: ', &
             & IRET_NC, TRIM(nf90_strerror(IRET_NC))
            CALL ABOR1('SUOROG: PLIT read failed from ' &
             & //TRIM(ISMRGNPFX(II))//'_pism2ece.nc')
          END IF
          IACTIVE_CELLS = COUNT(ZGLOB_PLIT > 0.5_JPRB)
          ZPLIT_SUM = SUM(ZGLOB_PLIT)
          WRITE(*,*) '>>> SUOROG: Read PLIT from ', TRIM(ISMRGNPFX(II)), &
           & ', range [', MINVAL(ZGLOB_PLIT), ',', MAXVAL(ZGLOB_PLIT), ']'
          WRITE(*,'(A,A,A,I0,A,F10.2)') &
           & '>>> SUOROG: ', TRIM(ISMRGNPFX(II)), &
           & ' PLIT active cells (>0.5): ', IACTIVE_CELLS, &
           & ', sum(plit): ', ZPLIT_SUM
        END IF
        IF (MYPROC == RPRC) THEN
          CALL DISGRID_SEND(YDGEOMETRY, 1, ZGLOB_PLIT, 2, ZLOC_PLIT_RGN)
        ELSE
          CALL DISGRID_RECV(YDGEOMETRY, RPRC, 1, ZLOC_PLIT_RGN, 2)
        END IF
        ! Accumulate across regions (non-overlapping, so no double-count)
        ZLOC_PLIT(:) = ZLOC_PLIT(:) + ZLOC_PLIT_RGN(:)
      END IF

      ! --- Close region file ---
      IF (MYPROC == RPRC) IRET_DUMMY = nf90_close(NCID_IN)

    END DO  ! regions

    ! Post-loop: orography spectral consistency
    IF (ECE_CPL_ISMM) THEN
      CALL REESPE(YDGEOMETRY,1,1,PSPOR,ZZOROG)
      ! TODO: check if extra lapse-rate corrections are needed for ISM orography
      CALL SPEREE(YDGEOMETRY,1,1,PSPOR,ZZOROG)
    END IF

    ! Post-loop: block and publish combined mask
    IF (ECE_LANDICE) THEN
      ALLOCATE(ZLOC_PLIT_B(NPROMA,NGPBLKS))
      ZLOC_PLIT_B(:,:) = 0._JPRB
      DO JKGLO=1,NGPTOT,NPROMA
        IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
        IBL=(JKGLO-1)/NPROMA+1
        ZLOC_PLIT_B(1:IEND,IBL)=MAX(0._JPRB,MIN(1._JPRB, &
         & ZLOC_PLIT(JKGLO:JKGLO+IEND-1)))
      ENDDO
      CALL SURFECE_SET_LANDICE_LOCAL(ZLOC_PLIT_B)
      DEALLOCATE(ZLOC_PLIT_B)
    END IF

    ! De-allocate global buffers
    IF (MYPROC == RPRC) THEN
      IF (ALLOCATED(ZGLOB_USURF)) DEALLOCATE(ZGLOB_USURF)
      IF (ALLOCATED(ZGLOB_PLIT))  DEALLOCATE(ZGLOB_PLIT)
    END IF

  ELSE
    CALL SURFECE_RESET_LANDICE()
  END IF

!*       1.2.3   COMPUTE FIRST ORDER DERIVATIVES OF OROGRAPHY.

  IF(.NOT.LGRADGP) THEN

    ZDIV(:)=0.0_JPRB
    CALL SPEUV(YDGEOMETRY,ZDIV,PSPOR,ZZOROGL,ZZOROGM,1,1,2)

  ELSE
#ifdef WITH_ATLAS
    NODE_COLUMNS = YDGEOMETRY%YRATLAS%FVM%NODE_COLUMNS()
    NODES        = YDGEOMETRY%YRATLAS%MESH%NODES()
    INODE_SIZE   = NODES%SIZE()
    GHOSTFIELD   = NODES%GHOST()
    CALL GHOSTFIELD%DATA(LLIS_GHOST)
    ALLOCATE(IPA(INODE_SIZE))
    IP = 0
    DO JL=1,INODE_SIZE
      IF(.NOT. LLIS_GHOST(JL)) THEN
        IP = IP+1
        IPA(JL) = IP
      ENDIF
    ENDDO
    VARFIELD  = NODE_COLUMNS%CREATE_FIELD(NAME="var_orog", KIND=ATLAS_REAL(JPRB),LEVELS=1)
    GRADFIELD = NODE_COLUMNS%CREATE_FIELD(NAME="grad_orog",KIND=ATLAS_REAL(JPRB),LEVELS=1,VARS=[2])
    CALL VARFIELD%DATA(ZVAR)
    CALL GRADFIELD%DATA(ZGRAD)
    DO JL=1,INODE_SIZE
      IF(.NOT. LLIS_GHOST(JL)) THEN
        IP = IPA(JL)
        ZVAR(1,JL) = ZZOROG(IP)
      ENDIF
    ENDDO
    CALL NODE_COLUMNS%HALO_EXCHANGE(VARFIELD)
    CALL FV_GRADIENT(YDGEOMETRY,ZVAR,ZGRAD)
    DO JL=1,INODE_SIZE
      IF(.NOT. LLIS_GHOST(JL)) THEN
        IP = IPA(JL)
        ZZOROGL(IP) = ZGRAD(1,1,JL)
        ZZOROGM(IP) = ZGRAD(2,1,JL)
      ENDIF
    ENDDO
    CALL VARFIELD%FINAL()
    CALL GRADFIELD%FINAL()
    CALL GHOSTFIELD%FINAL()
    CALL NODES%FINAL()
    CALL NODE_COLUMNS%FINAL()
#endif
  ENDIF

!*       1.2.4   COMPUTE SOME SECOND ORDER DERIVATIVES OF OROGRAPHY FOR NH.

  IF(YDGEOMETRY%LNONHYD_GEOM) THEN
    ALLOCATE(ZSWORK(NSPEC2))
    ALLOCATE(ZU(NGPTOT))
!       * orogll, oroglm.
    CALL REESPE(YDGEOMETRY,1,1,ZSWORK,ZZOROGL)
    CALL SPEUV(YDGEOMETRY,ZDIV,ZSWORK,ZZOROGLL,ZZOROGLM,1,1,2)
!       * orogmm:
    CALL REESPE(YDGEOMETRY,1,1,ZSWORK,ZZOROGM)
    CALL SPEUV(YDGEOMETRY,ZDIV,ZSWORK,ZU,ZZOROGMM,1,1,2)
    DEALLOCATE(ZU)
    DEALLOCATE(ZSWORK)
  ENDIF

ENDIF

!*       1.3    FILL YROROG
!               -----------

! Allocate and copy native NPROMA orography data structure
IF (.NOT. ALLOCATED(YDGEOMETRY%YROROG_B%OROG))  ALLOCATE(YDGEOMETRY%YROROG_B%OROG (NPROMA,NGPBLKS))
IF (.NOT. ALLOCATED(YDGEOMETRY%YROROG_B%OROGL)) ALLOCATE(YDGEOMETRY%YROROG_B%OROGL(NPROMA,NGPBLKS))
IF (.NOT. ALLOCATED(YDGEOMETRY%YROROG_B%OROGM)) ALLOCATE(YDGEOMETRY%YROROG_B%OROGM(NPROMA,NGPBLKS))

! * orography and horizontal gradient of orography:
DO JKGLO=1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  YDGEOMETRY%YROROG(IBL)%OROG (1:IEND)=ZZOROG (JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG(IBL)%OROGL(1:IEND)=ZZOROGL(JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG(IBL)%OROGM(1:IEND)=ZZOROGM(JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG_B%OROG (1:IEND,IBL)=ZZOROG (JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG_B%OROGL(1:IEND,IBL)=ZZOROGL(JKGLO:JKGLO+IEND-1)
  YDGEOMETRY%YROROG_B%OROGM(1:IEND,IBL)=ZZOROGM(JKGLO:JKGLO+IEND-1)
ENDDO

! * horizontal second-order derivatives of orography:
IF (YDGEOMETRY%LNONHYD_GEOM) THEN
  DO JKGLO=1,NGPTOT,NPROMA
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    YDGEOMETRY%YROROG(IBL)%OROGLL(1:IEND)=ZZOROGLL(JKGLO:JKGLO+IEND-1)
    YDGEOMETRY%YROROG(IBL)%OROGLM(1:IEND)=ZZOROGLM(JKGLO:JKGLO+IEND-1)
    YDGEOMETRY%YROROG(IBL)%OROGMM(1:IEND)=ZZOROGMM(JKGLO:JKGLO+IEND-1)
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUOROG',1,ZHOOK_HANDLE)
END SUBROUTINE SUOROG
