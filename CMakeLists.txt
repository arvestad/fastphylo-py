cmake_minimum_required(VERSION 3.15)
project(fastphylo LANGUAGES CXX C)

find_package(pybind11 REQUIRED)

set(FASTPHYLO_DIR src/cpp/fastphylo)

set(FASTPHYLO_SOURCES
    ${FASTPHYLO_DIR}/BitVector.cpp
    ${FASTPHYLO_DIR}/DistanceMatrix.cpp
    ${FASTPHYLO_DIR}/DistanceRow.cpp
    ${FASTPHYLO_DIR}/DNA_b128_String.cpp
    ${FASTPHYLO_DIR}/computeDistance_DNA_b128_String.cpp
    ${FASTPHYLO_DIR}/computeTAMURANEIDistance_DNA_b128_String.cpp
    ${FASTPHYLO_DIR}/Exception.cpp
    ${FASTPHYLO_DIR}/file_utils.cpp
    ${FASTPHYLO_DIR}/FloatDistanceMatrix.cpp
    ${FASTPHYLO_DIR}/InitAndPrintOn_utils.cpp
    ${FASTPHYLO_DIR}/Kimura2parameter.cpp
    ${FASTPHYLO_DIR}/NeighborJoining.cpp
    ${FASTPHYLO_DIR}/Object.cpp
    ${FASTPHYLO_DIR}/Sequence.cpp
    ${FASTPHYLO_DIR}/Sequences2DistanceMatrix.cpp
    ${FASTPHYLO_DIR}/SequenceTree.cpp
    ${FASTPHYLO_DIR}/ambiguity_nucleotide.cpp
    ${FASTPHYLO_DIR}/dna_pairwise_sequence_likelihood.cpp
    ${FASTPHYLO_DIR}/stl_utils.cpp
    ${FASTPHYLO_DIR}/string_compare.cpp
    ${FASTPHYLO_DIR}/TamuraNei.cpp
    ${FASTPHYLO_DIR}/sse2_wrapper.c
    ${FASTPHYLO_DIR}/xml_output_global.cpp
)

pybind11_add_module(_fastphylo
    src/cpp/bindings.cpp
    ${FASTPHYLO_SOURCES}
)

target_compile_features(_fastphylo PRIVATE cxx_std_17)
target_include_directories(_fastphylo PRIVATE ${FASTPHYLO_DIR})

# SSE2 for DNA_b128 on x86 (not needed/available on arm64)
if(NOT CMAKE_OSX_ARCHITECTURES MATCHES "arm64" AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i[3-6]86")
    target_compile_options(_fastphylo PRIVATE -msse2)
endif()

install(TARGETS _fastphylo DESTINATION fastphylo)
