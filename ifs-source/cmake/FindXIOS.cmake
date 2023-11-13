set(XIOS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/xios)

find_path(
    XIOS_INCLUDE_DIRECTORIES
    xios.mod
    PATHS ${XIOS_SOURCE_DIR}/inc
    REQUIRED
)
find_library(
    XIOS_LIBRARIES
    xios
    PATHS ${XIOS_SOURCE_DIR}/lib
    REQUIRED
)
