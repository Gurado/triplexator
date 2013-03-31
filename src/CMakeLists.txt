
cmake_minimum_required (VERSION 2.6)
project (triplexator)

if (OPENMP_FOUND)
   # We need a configuration header file.  First, configure this file, then
   # add the target directory to the include path.
   configure_file (${triplexator_SOURCE_DIR}/config.h.in
                   ${triplexator_BINARY_DIR}/config.h)
   
   # Enable OpenMP.
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

elseif (OPENMP_FOUND)
   message(STATUS "Building triplexator without OpenMP support.")
endif (OPENMP_FOUND)

include_directories(${triplexator_BINARY_DIR})
include_directories(${triplexator_SOURCE_DIR})

add_executable(triplexator triplexator.cpp 
				triplexator.h
				triplex.h
				helper.h
				find_index_qgrams.h
				gardener.h
				triplex_alphabet.h
				triplex_functors.h
				triplex_pattern.h
)

# The boost component iostreams is optional and only be used to compress output files on request
if (Boost_IOSTREAMS_FOUND)
	set(Boost_USE_STATIC_LIBS ON)
	include_directories(${Boost_INCLUDE_DIRS})
	target_link_libraries(triplexator ${Boost_IOSTREAMS_LIBRARIES})
	message("Boost_IOSTREAMS_LIBRARY_RELEASE = ${Boost_IOSTREAMS_LIBRARY_RELEASE}")
endif (Boost_IOSTREAMS_FOUND)

if(UNIX AND NOT APPLE)
	target_link_libraries(triplexator rt)
endif()

install(TARGETS triplexator RUNTIME DESTINATION bin COMPONENT applications) 

