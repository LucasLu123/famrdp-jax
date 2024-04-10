##############################################################################
# HO cmake build system
##############################################################################
CMAKE_MINIMUM_REQUIRED( VERSION 3.2 FATAL_ERROR )

#SET(CMAKE_FORTRAN_COMPILER_ENV_VAR gfortran)

#SET(CMAKE_Fortran_COMPILER /usr/bin/gfortran)

#get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

PROJECT(HOSTA Fortran)

ENABLE_LANGUAGE(Fortran)

IF (UNIX)
    MESSAGE(STATUS "current platform: Linux ")
ELSEIF (WIN32)
    MESSAGE(STATUS "current platform: Windows")
ELSE ()
    MESSAGE(STATUS "other platform: ${CMAKE_SYSTEM_NAME}")
ENDIF ()

#SET(CMAKE_CXX_STANDARD_REQUIRED ON)

#SET(CMAKE_CXX_STANDARD 11)
SET(MAKECMDGOALS mpi)
SET(EXEPRE HOSTA)

IF (WIN32)
    SET(FFLAGS "-Dlogo")

    SET(MPI_PATH "C:/Program Files/MPICH2")
    SET(mpilib fmpich2.lib)
ELSE()
    #SET(FFLAGS "-Dlogo -O3 -xHOST  -no-prec-div  -cpp -r8 -no-vec ${OMP_FLAG} ${MPI_FLAG}")
    SET(FFLAGS "-DO3 -cpp")

    SET(MPI_PATH "/usr/local/")
    SET(mpilib "")
ENDIF()

IF(MAKECMDGOALS STREQUAL "seq")
    SET(EXESUF seq)
    SET(FD ifort)
    SET(LD ifort)
    IF (UNIX)
        SET(CMAKE_Fortran_COMPILER /usr/bin/gfortran)
    ENDIF()
    SET(MPI_FLAG "")
    SET(OMP_FLAG "")
ENDIF()

IF(MAKECMDGOALS STREQUAL "mpi")
    SET(EXESUF mpi)
    SET(FD mpif90)
    SET(LD mpif90)
    SET(MPI_FLAG "-DPARALLEL")    
    SET(OMP_FLAG "")

    IF (WIN32)
    ELSE()
        SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90)
    ENDIF()
    INCLUDE_DIRECTORIES(${MPI_PATH}/include)
ENDIF()

IF(MAKECMDGOALS STREQUAL "omp")
    SET(EXESUF omp)
    SET(FD mpif90)
    SET(LD mpif90)
    SET(MPI_FLAG "")
    SET(OMP_FLAG "-DOMP_IMP -openmp")

    IF (WIN32)
    ELSE()
        SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90)
    ENDIF()    
    INCLUDE_DIRECTORIES(${MPI_PATH}/include)
ENDIF()

IF(MAKECMDGOALS STREQUAL "mpiomp")
    SET(EXESUF mpiomp)
    SET(FD mpif90)
    SET(LD mpif90)
    SET(MPI_FLAG "-DPARALLEL")
    SET(OMP_FLAG "-DOMP_IMP -openmp")

    IF (WIN32)
    ELSE()
        SET(CMAKE_Fortran_COMPILER /usr/bin/mpif90)
    ENDIF()    
    INCLUDE_DIRECTORIES(${MPI_PATH}/include)
ENDIF()

SET(SET_OPT_FLAGS1 "-DPARAM_OPT")
SET(SET_OPT_FLAGS2 "-DMemOPT1")
SET(SET_OPT_FLAGS3 "-DMemOPT2")
SET(SET_OPT_FLAGS4 "-DAFLUX_VIS_OPT")

SET(target_name ${EXEPRE}.${EXESUF})

ADD_SUBDIRECTORY (src)