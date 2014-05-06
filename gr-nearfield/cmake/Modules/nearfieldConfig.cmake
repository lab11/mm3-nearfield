INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_NEARFIELD nearfield)

FIND_PATH(
    NEARFIELD_INCLUDE_DIRS
    NAMES nearfield/api.h
    HINTS $ENV{NEARFIELD_DIR}/include
        ${PC_NEARFIELD_INCLUDEDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/include
          /usr/local/include
          /usr/include
)

FIND_LIBRARY(
    NEARFIELD_LIBRARIES
    NAMES gnuradio-nearfield
    HINTS $ENV{NEARFIELD_DIR}/lib
        ${PC_NEARFIELD_LIBDIR}
    PATHS ${CMAKE_INSTALL_PREFIX}/lib
          ${CMAKE_INSTALL_PREFIX}/lib64
          /usr/local/lib
          /usr/local/lib64
          /usr/lib
          /usr/lib64
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NEARFIELD DEFAULT_MSG NEARFIELD_LIBRARIES NEARFIELD_INCLUDE_DIRS)
MARK_AS_ADVANCED(NEARFIELD_LIBRARIES NEARFIELD_INCLUDE_DIRS)

