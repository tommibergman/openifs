set(OASIS_ARCH_DIR ${CMAKE_SOURCE_DIR}/oasis/arch_ecearth)

find_path(
    OASIS_INCLUDE_DIRECTORIES
    NAMES mod_oasis.mod mct_mod.mod m_mpout.mod remap_bicubic_reduced.mod
    PATHS ${OASIS_ARCH_DIR}/include
    REQUIRED
)

find_library(
    psmile_lib
    NAMES psmile psmile.MPI1
    PATHS ${OASIS_ARCH_DIR}/lib
    REQUIRED
)
find_library(
    mct_lib
    mct
    PATHS ${OASIS_ARCH_DIR}/lib
    REQUIRED
)
find_library(
    mpeu_lib
    mpeu
    PATHS ${OASIS_ARCH_DIR}/lib
    REQUIRED
)
find_library(
    scrip_lib
    scrip
    PATHS ${OASIS_ARCH_DIR}/lib
    REQUIRED
)
set(OASIS_LIBRARIES ${psmile_lib} ${mct_lib} ${mpeu_lib} ${scrip_lib})
