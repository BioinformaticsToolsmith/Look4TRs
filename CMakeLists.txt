cmake_minimum_required (VERSION 2.8)
project (Look4TRs)



if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 7.0
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7.0)
        message(FATAL_ERROR "GCC version must be at least 7.0!")
    endif()
endif()

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})


add_library(SOURCES 
${CMAKE_SOURCE_DIR}/src/align/GlobAlignE.cpp
${CMAKE_SOURCE_DIR}/src/align/RepeatAlignE.cpp
${CMAKE_SOURCE_DIR}/src/cluster/DivergencePoint.cpp
${CMAKE_SOURCE_DIR}/src/cluster/Feature.cpp
${CMAKE_SOURCE_DIR}/src/cluster/Progress.cpp
${CMAKE_SOURCE_DIR}/src/cluster/Selector.cpp
${CMAKE_SOURCE_DIR}/src/exception/FileDoesNotExistException.cpp
${CMAKE_SOURCE_DIR}/src/exception/InvalidInputException.cpp
${CMAKE_SOURCE_DIR}/src/exception/InvalidOperationException.cpp
${CMAKE_SOURCE_DIR}/src/exception/InvalidStateException.cpp
${CMAKE_SOURCE_DIR}/src/matrix/GLM.cpp
${CMAKE_SOURCE_DIR}/src/matrix/Matrix.cpp
${CMAKE_SOURCE_DIR}/src/motif/FindMotif.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromListMaker.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/Chromosome.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeOneDigit.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeRandom.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeTR.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeReadTR.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/ChromosomeSelfTR.cpp
${CMAKE_SOURCE_DIR}/src/nonltr/HMM.cpp
${CMAKE_SOURCE_DIR}/src/satellites/Client_Scanner.cpp
${CMAKE_SOURCE_DIR}/src/satellites/Client_Trainer.cpp
${CMAKE_SOURCE_DIR}/src/satellites/Executor.cpp
${CMAKE_SOURCE_DIR}/src/satellites/IClient.cpp
${CMAKE_SOURCE_DIR}/src/satellites/ScorerAdjusted.cpp
${CMAKE_SOURCE_DIR}/src/satellites/ScorerSat.cpp
${CMAKE_SOURCE_DIR}/src/train/Predictor.cpp
${CMAKE_SOURCE_DIR}/src/utility/Location.cpp
${CMAKE_SOURCE_DIR}/src/utility/RegionList.cpp
${CMAKE_SOURCE_DIR}/src/utility/Util.cpp
${CMAKE_SOURCE_DIR}/src/mutate/HandleSeq.cpp
${CMAKE_SOURCE_DIR}/src/mutate/MultiMute.cpp
${CMAKE_SOURCE_DIR}/src/mutate/SingMute.cpp
)


set(CMAKE_CXX_COMPILER g++)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "-fopenmp -O3 -fmessage-length=0")

add_executable(Look4TRs ${CMAKE_SOURCE_DIR}/src/satellites/Look4TRs.cpp)
target_link_libraries(Look4TRs SOURCES)


