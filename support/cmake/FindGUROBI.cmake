# Find GUROBI library
#
# The original file was downloaded from
# https://support.gurobi.com/hc/en-us/articles/360039499751-How-do-I-use-CMake-to-build-Gurobi-C-C-projects
# on April 22, 2024
#
# Modified to check GUROBI_INCLUDE_DIRS and some lib subdirs

find_path(GUROBI_INCLUDE_DIRS
		NAMES gurobi_c.h
		HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
		PATH_SUFFIXES include)

find_library(GUROBI_LIBRARY
		NAMES gurobi gurobi120 gurobi110 gurobi100
		HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
		PATH_SUFFIXES lib lib/linux64 lib/osx64 lib/win64)

if(CXX)
		if(MSVC)
				set(MSVC_YEAR "2017")

				if(MT)
						set(M_FLAG "mt")
				else()
						set(M_FLAG "md")
				endif()

				find_library(GUROBI_CXX_LIBRARY
						NAMES gurobi_c++${M_FLAG}${MSVC_YEAR}
						HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
						PATH_SUFFIXES lib)
				find_library(GUROBI_CXX_DEBUG_LIBRARY
						NAMES gurobi_c++${M_FLAG}d${MSVC_YEAR}
						HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
						PATH_SUFFIXES lib)
		else()
				find_library(GUROBI_CXX_LIBRARY
						NAMES gurobi_c++
						HINTS ${GUROBI_DIR} $ENV{GUROBI_HOME}
						PATH_SUFFIXES lib)
				set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_CXX_LIBRARY})
		endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_INCLUDE_DIRS)
