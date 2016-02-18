#******************************************************************************
#
# MantaFlow fluid solver framework
#
# Copyright 2011-2015 Tobias Pfaff, Nils Thuerey
#
# This program is free software, distributed under the terms of the
# GNU General Public License (GPL)
# http://www.gnu.org/licenses
#
# modified by David Ullmann 2016
#
#******************************************************************************

project (Mantaflow)
cmake_minimum_required(VERSION 2.8)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/source/cmake/")
SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
SET(VERBOSE 1)

SET(MANTAVERSION "0.8")



#******************************************************************************
# Default paths
#******************************************************************************
IF(WIN32)
	set(CMAKE_LIBRARY_PATH "C:/Program Files (x86)/Windows Kits/8.0/Lib/win8/um/x64")
	SET(CMAKE_PREFIX_PATH ${WIN_QT_PATH})
ENDIF()

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
    string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${EXECUTABLE_OUTPUT_PATH})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${LIBRARY_OUTPUT_PATH})
	set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${LIBRARY_OUTPUT_PATH})
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )
IF(APPLE)
	IF(NOT CMAKE_PREFIX_PATH)
		SET(CMAKE_PREFIX_PATH "/usr/local/Cellar/qt5/5.5/") #  mac/homebrew
		#SET(CMAKE_PREFIX_PATH "/home/myname/qt/5.5/clang_64") #  other...
	ENDIF()
ENDIF()



#******************************************************************************
# Setup default params
#******************************************************************************
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE "Release")
ELSE()
	MESSAGE("Build-type: '${CMAKE_BUILD_TYPE}'")
ENDIF()

# compilation versions
OPTION(DEBUG "Enable debug compilation" OFF)
OPTION(DOUBLEPRECISION "Compile with double floating point precision" OFF)

# make sure debug settings match...
IF(NOT DEBUG)
	IF(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
		SET(DEBUG 1)
	ENDIF()
ENDIF()
IF(DEBUG)
	SET(CMAKE_BUILD_TYPE "Debug")
	add_definitions ( -DDEBUG )
ENDIF()

MESSAGE(STATUS "")
MESSAGE(STATUS "Options - "
	" -DDEBUG='${DEBUG}' "
	" -DDOUBLEPRECISION='${DOUBLEPRECISION}' "
	)
MESSAGE(STATUS "")



#******************************************************************************
# Setup main project
#******************************************************************************
SET(F_LIBS "" )
SET(F_LIB_PATHS)
SET(F_LINKADD "") # additional linker flags, not a list
SET(SILENT_SOURCES)

# need pre-processing
SET(MANTA_SOURCES
	source/conjugategrad.cpp
	source/edgecollapse.cpp
	source/fastmarch.cpp
	source/fileio.cpp
	source/fluidsolver.cpp
	source/general.cpp
	source/grid.cpp
	source/kernel.cpp
	source/levelset.cpp
	source/mesh.cpp
	source/movingobs.cpp
	source/noisefield.cpp
	source/particle.cpp
	source/shapes.cpp
	source/test.cpp
	source/timing.cpp
	source/turbulencepart.cpp
	source/vortexpart.cpp
	source/vortexsheet.cpp)

SET(MANTA_PLUGIN_SOURCES
	source/plugin/advection.cpp
	source/plugin/extforces.cpp
	source/plugin/flip.cpp
	source/plugin/initplugins.cpp
	source/plugin/kepsilon.cpp
	source/plugin/meshplugins.cpp
	source/plugin/pressure.cpp
	source/plugin/vortexplugins.cpp
	source/plugin/waveletturbulence.cpp
	source/plugin/waves.cpp)

SET(MANTA_UTIL_SOURCES
	source/util/simpleimage.cpp
	source/util/vectorbase.cpp)

SET(MANTA_HEADERS
	include/commonkernels.h
	include/conjugategrad.h
	include/edgecollapse.h
	include/fastmarch.h
	include/fileio.h
	include/fluidsolver.h
	include/general.h
	include/grid.h
	include/kernel.h
	include/levelset.h
	include/mant.h
	include/mesh.h
	include/movingobs.h
	include/noisefield.h
	include/particle.h
	include/pclass.h
	include/shapes.h
	include/timing.h
	include/turbulencepart.h
	include/vortexpart.h
	include/vortexsheet.h)

SET(MANTA_PLUGIN_HEADERS
	include/plugin/advection.h
	include/plugin/extforces.h
	include/plugin/flip.h
	include/plugin/initplugins.h
	include/plugin/pressure.h)

SET(MANTA_UTIL_HEADERS
	include/util/integrator.h
	include/util/interpol.h
	include/util/interpolHigh.h
	include/util/mcubes.h
	include/util/quaternion.h
	include/util/randomstream.h
	include/util/simpleimage.h
	include/util/solvana.h
	include/util/vectorbase.h)

SET(ZLIB_SOURCES
	dependencies/zlib-1.2.8/adler32.c
	dependencies/zlib-1.2.8/compress.c
	dependencies/zlib-1.2.8/crc32.c
	dependencies/zlib-1.2.8/deflate.c
	dependencies/zlib-1.2.8/gzclose.c
	dependencies/zlib-1.2.8/gzlib.c
	dependencies/zlib-1.2.8/gzread.c
	dependencies/zlib-1.2.8/gzwrite.c
	dependencies/zlib-1.2.8/inflate.c
	dependencies/zlib-1.2.8/infback.c
	dependencies/zlib-1.2.8/inftrees.c
	dependencies/zlib-1.2.8/inffast.c
	dependencies/zlib-1.2.8/trees.c
	dependencies/zlib-1.2.8/uncompr.c
	dependencies/zlib-1.2.8/zutil.c)

# include dirs
SET(INCLUDE_PATHS
	include
	include/plugin
	include/util
)

# Multithreading with OpenMP is on by default
add_definitions( -DMT)
add_definitions( -DOPENMP )
if (WIN32)
	add_definitions( /openmp)
else()
	add_definitions(-fopenmp)
	SET(F_LINKADD "${F_LINKADD} -fopenmp ")
endif()

# increase FP precision?
if(DOUBLEPRECISION)
	add_definitions(-DFLOATINGPOINT_PRECISION=2)
endif()



#******************************************************************************
# Setup Zlib static library
#******************************************************************************
if (1)
	# compiler flags
	IF(NOT WIN32)
		IF(DEBUG)
			add_definitions( -DDEBUG )
			add_definitions( -O0 -fno-inline )
			# additional c++ warning
			SET(CMAKE_CXX_FLAGS_DEBUG " -ggdb  -Wall")

			# stricter: no optimizations and inlining, comment out if needed...
			#SET(CMAKE_CXX_FLAGS_DEBUG " -O0 -fno-inline ${CMAKE_CXX_FLAGS_DEBUG}   ")
		ELSE()
			# non-debug, optimized version
			add_definitions( -O3 -Wall)
		ENDIF()
	ELSE()
		# get rid of some MSVC warnings
		add_definitions( /wd4018 /wd4146 /wd4800 )

		# for zlib
		add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
		add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)

		# unsigned to signed conversions
		add_definitions( /wd4267 )
		# double <> single precision
		add_definitions( /wd4244 /wd4305 )
		# disable warnings for unsecure functions
		add_definitions( /D "_CRT_SECURE_NO_WARNINGS" )

		# enable when using Qt creator:
		#add_definitions(/FS)
	ENDIF()
	
	# default: build from source
	SET(ZLIB_LIBRARY_NAME zlib)

	set(ZLIB_ADDFLAGS "-Dverbose=-1")
	if(NOT WIN32)
		set(ZLIB_ADDFLAGS "-Wno-implicit-function-declaration -Dverbose=-1")
	endif()
	add_definitions(${ZLIB_ADDFLAGS})
	
	list(APPEND INCLUDE_PATHS dependencies/zlib-1.2.8)

	add_library(${ZLIB_LIBRARY_NAME} STATIC ${ZLIB_SOURCES})
	
else()
	include(FindZLIB)
	list(APPEND INCLUDE_PATHS ${ZLIB_INCLUDE_DIR})
	list(APPEND F_LIBS ${ZLIB_LIBRARIES})
endif()



#******************************************************************************
# Setup Mantaflow static library
#******************************************************************************
# compiler flags
IF(NOT WIN32)
	IF(DEBUG)
		add_definitions( -DDEBUG )
		add_definitions( -O0 -fno-inline )
		# additional c++ warning
		SET(CMAKE_CXX_FLAGS_DEBUG " -ggdb  -Wall")

		# stricter: no optimizations and inlining, comment out if needed...
		#SET(CMAKE_CXX_FLAGS_DEBUG " -O0 -fno-inline ${CMAKE_CXX_FLAGS_DEBUG}   ")
	ELSE()
		# non-debug, optimized version
		add_definitions( -O3 -Wall)
	ENDIF()
ELSE()
	# get rid of some MSVC warnings
	add_definitions( /wd4018 /wd4146 /wd4800 )

	# for zlib
	add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
	add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)

	# unsigned to signed conversions
	add_definitions( /wd4267 )
	# double <> single precision
	add_definitions( /wd4244 /wd4305 )
	# disable warnings for unsecure functions
	add_definitions( /D "_CRT_SECURE_NO_WARNINGS" )

	# enable when using Qt creator:
	#add_definitions(/FS)
ENDIF()

SET(MANTA_LIBRARY_NAME Mantaflow)
SET(CMAKE_EXE_LINKER_FLAGS_DEBUG   "${CMAKE_EXE_LINKER_FLAGS_DEBUG}   ${F_LINKADD} ")
SET(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${F_LINKADD} ")

include_directories(${INCLUDE_PATHS})

if (WIN32)
	add_library(${MANTA_LIBRARY_NAME} STATIC ${MANTA_SOURCES} ${MANTA_PLUGIN_SOURCES} ${MANTA_UTIL_SOURCES})

	# make nice folders for Visual Studio
	source_group(Plugins FILES ${MANTA_PLUGIN_SOURCES} ${MANTA_PLUGIN_HEADERS})
	source_group(Util FILES ${MANTA_UTIL_SOURCES} ${MANTA_UTIL_HEADERS})
else()
endif()



#******************************************************************************
# Sample executables
#******************************************************************************
SET(TESTS_SOURCES
	source/tests/DropObjMesh.cpp
	source/tests/DropObjMesh_resurfaced.cpp
	source/tests/flip01_simple.cpp
	source/tests/flip02_surface.cpp
	source/tests/flip03_gen.cpp
	source/tests/flip04_adaptDt.cpp
	source/tests/flip05_secOrderBnd.cpp
	source/tests/plume_2d.cpp
	source/tests/plume_adaptDt.cpp) 
 
foreach(testfile ${TESTS_SOURCES})
	STRING(REPLACE ".cpp" "" executable ${testfile}) 
	string(REPLACE "source/tests/" "" testname ${executable})
    add_executable(${testname} ${testfile})
	target_link_libraries(${testname} ${MANTA_LIBRARY_NAME} ${ZLIB_LIBRARY_NAME})
endforeach()
