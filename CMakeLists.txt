# 3.18, not pybind11's own 3.15 floor: matches FastPhylo's own
# cmake_minimum_required (src/c++/CMakeLists.txt's find_package(BLAS/
# LAPACK) needs the BLAS::BLAS/LAPACK::LAPACK imported targets, only
# reliable from CMake 3.18 on) - the submodule below re-declares 3.18
# itself regardless, but stating it here too makes the real floor
# obvious without having to open extern/fastphylo to find it.
cmake_minimum_required(VERSION 3.18)
project(fastphylo LANGUAGES CXX C)

find_package(pybind11 REQUIRED)

# fastphylo-py integration plan, Phase 2: FastPhylo is now a git
# submodule (pinned to a tag once one exists past this integration
# work - see the submodule's own detached-HEAD commit in the meantime)
# instead of a vendored, independently-drifting copy of its source
# (the old src/cpp/fastphylo/ - deleted once bindings.cpp is ported,
# see Phase 3/5). OFF/OFF: this extension only needs the fastphylo
# library target, not FastPhylo's own CLI apps or test executables.
set(FASTPHYLO_BUILD_APPS OFF CACHE BOOL "" FORCE)
set(FASTPHYLO_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(extern/fastphylo)

pybind11_add_module(_fastphylo src/cpp/bindings.cpp)

target_compile_features(_fastphylo PRIVATE cxx_std_17)
# fastphylo's PUBLIC target_include_directories (its own include/ root,
# plus Eigen/BLAS/LAPACK's) propagate automatically - no manual include
# path needed here the way the vendored FASTPHYLO_DIR used to require.
target_link_libraries(_fastphylo PRIVATE fastphylo)

install(TARGETS _fastphylo DESTINATION fastphylo)
