cmake_minimum_required(VERSION 3.18)

project(PnetCDF VERSION 1.15.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCSourceCompiles)
include(CheckFunctionExists)
include(CheckIncludeFile)
include(CheckSymbolExists)
include(CheckTypeSize)
include(CTest)
include(GNUInstallDirs)
include(TestBigEndian)

option(PNETCDF_ENABLE_ERANGE_FILL "Use fill value when type conversion returns NC_ERANGE" ON)
option(PNETCDF_ENABLE_NULL_BYTE_HEADER_PADDING "Check null-byte padding in classic headers" OFF)
option(PNETCDF_ENABLE_PROFILING "Enable internal profiling counters" OFF)
option(PNETCDF_ENABLE_DEBUG "Enable PnetCDF debug mode" OFF)
option(PNETCDF_ENABLE_SUBFILING "Enable subfiling support" OFF)
option(PNETCDF_BUILD_PYTHON_COMMUX "Build the thin Python binding for the commux C ABI" ON)
set(PNETCDF_COMM_BACKEND "commux" CACHE STRING "Communication backend: mpi, commux, or auto")
set_property(CACHE PNETCDF_COMM_BACKEND PROPERTY STRINGS mpi commux auto)
set(COMMUX_DIR "" CACHE PATH "Installed Python commux package directory")
set(COMMUX_UCX_INCLUDE_DIR "" CACHE PATH "UCX include directory used by commux")
set(TORCH_LIBRARY_DIR "" CACHE PATH "Torch library directory used by commux")
set(TORCH_INSTALL_PREFIX "" CACHE PATH "Torch package directory used by commux")
set(TORCH_CXX11_ABI "" CACHE STRING "Torch GLIBCXX C++11 ABI setting")
set(Python3_EXECUTABLE "" CACHE FILEPATH "Python used to locate commux")

set(PNETCDF_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
set(PNETCDF_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
set(PNETCDF_INSTALL_PYDIR "${CMAKE_INSTALL_LIBDIR}")
if(SKBUILD)
  set(PNETCDF_INSTALL_INCLUDEDIR "pinc/include")
  set(PNETCDF_INSTALL_LIBDIR "pinc/lib")
  set(PNETCDF_INSTALL_PYDIR "pinc")
endif()

set(PNETCDF_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(PNETCDF_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(PNETCDF_VERSION_SUB ${PROJECT_VERSION_PATCH})
set(PNETCDF_VERSION_PRE "")
string(TIMESTAMP PNETCDF_RELEASE_DATE "%B %d, %Y")

if(PNETCDF_COMM_BACKEND STREQUAL "auto" OR PNETCDF_COMM_BACKEND STREQUAL "commux")
  if(NOT COMMUX_DIR)
    if(NOT Python3_EXECUTABLE)
      find_package(Python3 QUIET COMPONENTS Interpreter)
      if(Python3_Interpreter_FOUND)
        set(Python3_EXECUTABLE "${Python3_EXECUTABLE}")
      endif()
    endif()
    if(Python3_EXECUTABLE)
      execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c "import commux, os; print(os.path.dirname(commux.__file__))"
        OUTPUT_VARIABLE COMMUX_DIR_DETECTED
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
      if(COMMUX_DIR_DETECTED)
        set(COMMUX_DIR "${COMMUX_DIR_DETECTED}" CACHE PATH "Installed Python commux package directory" FORCE)
      endif()
    endif()
  endif()
endif()

if((PNETCDF_COMM_BACKEND STREQUAL "auto" OR PNETCDF_COMM_BACKEND STREQUAL "commux") AND NOT TORCH_LIBRARY_DIR)
  if(Python3_EXECUTABLE)
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import torch, os; root=os.path.dirname(torch.__file__); print(root); print(os.path.join(root, 'lib')); print(int(torch._C._GLIBCXX_USE_CXX11_ABI))"
      OUTPUT_VARIABLE TORCH_INFO_DETECTED
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(TORCH_INFO_DETECTED)
      string(REPLACE "\n" ";" TORCH_INFO_LINES "${TORCH_INFO_DETECTED}")
      list(GET TORCH_INFO_LINES 0 TORCH_INSTALL_PREFIX_DETECTED)
      list(GET TORCH_INFO_LINES 1 TORCH_LIBRARY_DIR_DETECTED)
      list(GET TORCH_INFO_LINES 2 TORCH_CXX11_ABI_DETECTED)
      set(TORCH_INSTALL_PREFIX "${TORCH_INSTALL_PREFIX_DETECTED}" CACHE PATH "Torch package directory used by commux" FORCE)
      set(TORCH_LIBRARY_DIR "${TORCH_LIBRARY_DIR_DETECTED}" CACHE PATH "Torch library directory used by commux" FORCE)
      set(TORCH_CXX11_ABI "${TORCH_CXX11_ABI_DETECTED}" CACHE STRING "Torch GLIBCXX C++11 ABI setting" FORCE)
    endif()
  endif()
endif()

set(ENABLE_COMMUX 0)
set(PNETCDF_COMM_BACKEND_VALUE 1)
if(COMMUX_DIR
   AND EXISTS "${COMMUX_DIR}/include/commux/process_group_ucx.hpp"
   AND EXISTS "${COMMUX_DIR}/lib/libcommux.so")
  set(ENABLE_COMMUX 1)
  if(PNETCDF_COMM_BACKEND STREQUAL "commux")
    set(PNETCDF_COMM_BACKEND_VALUE 2)
  endif()
elseif(PNETCDF_COMM_BACKEND STREQUAL "commux")
  message(FATAL_ERROR "PNETCDF_COMM_BACKEND=commux requested, but COMMUX_DIR is not usable")
endif()

if(PNETCDF_COMM_BACKEND STREQUAL "auto")
  if(ENABLE_COMMUX)
    set(PNETCDF_COMM_BACKEND "commux")
    set(PNETCDF_COMM_BACKEND_VALUE 2)
  else()
    set(PNETCDF_COMM_BACKEND "mpi")
  endif()
endif()

set(ENABLE_FORTRAN 0)
set(ENABLE_CXX 0)
set(ENABLE_ERANGE_FILL 0)
if(PNETCDF_ENABLE_ERANGE_FILL)
  set(ENABLE_ERANGE_FILL 1)
endif()
set(ENABLE_SUBFILING 0)
if(PNETCDF_ENABLE_SUBFILING)
  set(ENABLE_SUBFILING 1)
endif()
set(ENABLE_THREAD_SAFE 0)
set(ENABLE_NETCDF4 0)
set(ENABLE_ADIOS 0)
set(ENABLE_GIO 0)
set(ENABLE_BURST_BUFFER 0)
set(ENABLE_NULL_BYTE_HEADER_PADDING 0)
if(PNETCDF_ENABLE_NULL_BYTE_HEADER_PADDING)
  set(ENABLE_NULL_BYTE_HEADER_PADDING 1)
endif()
set(PNETCDF_DEBUG 0)
if(PNETCDF_ENABLE_DEBUG)
  set(PNETCDF_DEBUG 1)
endif()
set(PNETCDF_PROFILING 0)
if(PNETCDF_ENABLE_PROFILING)
  set(PNETCDF_PROFILING 1)
endif()
set(RELAX_COORD_BOUND 1)
set(IN_PLACE_SWAP 0)

if(PNETCDF_COMM_BACKEND STREQUAL "commux")
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  if(NOT COMMUX_UCX_INCLUDE_DIR)
    find_path(COMMUX_UCX_INCLUDE_DIR
      NAMES ucp/api/ucp.h
      HINTS
        "${COMMUX_DIR}/include"
        "$ENV{UCX_ROOT}/include")
  endif()
  if(NOT COMMUX_UCX_INCLUDE_DIR
     AND EXISTS "/home/chengcli/scix/repos/commux/build/_deps/commux-ucx/install/include/ucp/api/ucp.h")
    set(COMMUX_UCX_INCLUDE_DIR
      "/home/chengcli/scix/repos/commux/build/_deps/commux-ucx/install/include"
      CACHE PATH "UCX include directory used by commux" FORCE)
  endif()
  if(NOT COMMUX_UCX_INCLUDE_DIR)
    message(FATAL_ERROR
      "Could not find ucp/api/ucp.h for commux. Set COMMUX_UCX_INCLUDE_DIR.")
  endif()
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src/include")
  configure_file(src/include/pnetcdf_comm.h.in "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf_comm.h" @ONLY)

  add_library(pnc src/comm/pnc_commux_runtime.cpp)
  target_sources(pnc PRIVATE src/comm/pnc_netcdf_io.cpp)
  target_include_directories(pnc
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src/include>"
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
    PRIVATE
      "${COMMUX_DIR}/include"
      "${COMMUX_UCX_INCLUDE_DIR}"
      "${TORCH_INSTALL_PREFIX}/include"
      "${TORCH_INSTALL_PREFIX}/include/torch/csrc/api/include")
  target_compile_definitions(pnc PRIVATE _GLIBCXX_USE_CXX11_ABI=${TORCH_CXX11_ABI})
  target_link_libraries(pnc PRIVATE
    "${COMMUX_DIR}/lib/libcommux.so"
    "${TORCH_LIBRARY_DIR}/libtorch.so"
    "${TORCH_LIBRARY_DIR}/libtorch_cpu.so"
    "${TORCH_LIBRARY_DIR}/libc10.so"
    netcdf)
  if(EXISTS "${COMMUX_DIR}/lib/libcommux_io.so")
    target_link_libraries(pnc PRIVATE "${COMMUX_DIR}/lib/libcommux_io.so")
  endif()
  if(COMMUX_DIR)
    set(_pnc_commux_runtime_rpath "${COMMUX_DIR}/lib;${TORCH_LIBRARY_DIR}")
    if(SKBUILD)
      # In the wheel, commux and torch are sibling site-packages while libpnc
      # lives in pinc/lib, so reach them via $ORIGIN-relative paths that stay
      # valid regardless of the install prefix.
      set(_pnc_commux_runtime_rpath "$ORIGIN;$ORIGIN/../../commux/lib;$ORIGIN/../../torch/lib")
    endif()
    set_target_properties(pnc PROPERTIES
      BUILD_RPATH "${COMMUX_DIR}/lib;${TORCH_LIBRARY_DIR}"
      INSTALL_RPATH "${_pnc_commux_runtime_rpath}")
  endif()
  set_target_properties(pnc PROPERTIES
    OUTPUT_NAME pnc
    VERSION 0.1.0
    SOVERSION 0)

  add_library(pnetcdf src/comm/pnc_comm_commux.c)
  target_compile_definitions(pnetcdf PRIVATE PNETCDF_COMMUX_ONLY=1)
  target_include_directories(pnetcdf
    PUBLIC
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src/include>"
      "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>"
      "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
  if(COMMUX_DIR)
    target_include_directories(pnetcdf PRIVATE "${COMMUX_DIR}/include")
    target_link_libraries(pnetcdf PRIVATE "${COMMUX_DIR}/lib/libcommux.so")
    if(TORCH_LIBRARY_DIR
       AND EXISTS "${TORCH_LIBRARY_DIR}/libtorch.so"
       AND EXISTS "${TORCH_LIBRARY_DIR}/libtorch_cpu.so"
       AND EXISTS "${TORCH_LIBRARY_DIR}/libc10.so")
      target_link_libraries(pnetcdf PRIVATE
        "${TORCH_LIBRARY_DIR}/libtorch.so"
        "${TORCH_LIBRARY_DIR}/libtorch_cpu.so"
        "${TORCH_LIBRARY_DIR}/libc10.so")
    endif()
    set(_pnc_commux_rpath "${COMMUX_DIR}/lib;${TORCH_LIBRARY_DIR}")
    if(SKBUILD)
      # See pnc above: sibling site-packages, reached relative to pinc/lib.
      set(_pnc_commux_rpath "$ORIGIN;$ORIGIN/../../commux/lib;$ORIGIN/../../torch/lib")
    endif()
    set_target_properties(pnetcdf PROPERTIES
      BUILD_RPATH "${COMMUX_DIR}/lib;${TORCH_LIBRARY_DIR}"
      INSTALL_RPATH "${_pnc_commux_rpath}")
  endif()
  target_link_libraries(pnetcdf PUBLIC pnc)
  set_target_properties(pnetcdf PROPERTIES
    OUTPUT_NAME pnetcdf
    VERSION 8.0.0
    SOVERSION 8)

  install(TARGETS pnetcdf
    ARCHIVE DESTINATION "${PNETCDF_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${PNETCDF_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
  install(TARGETS pnc
    ARCHIVE DESTINATION "${PNETCDF_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${PNETCDF_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf_comm.h"
    DESTINATION "${PNETCDF_INSTALL_INCLUDEDIR}")
  install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/include/pnetcdf_commux.h"
    DESTINATION "${PNETCDF_INSTALL_INCLUDEDIR}")
  if(PNETCDF_BUILD_PYTHON_COMMUX)
    find_package(Python3 QUIET COMPONENTS Development.Module)
    if(Python3_Development.Module_FOUND)
      add_library(pnetcdf_commux_py MODULE python/csrc/pnetcdf_commux.cpp)
      set(_pnc_py_install_rpath "$ORIGIN/lib;$ORIGIN;${COMMUX_DIR}/lib;${TORCH_LIBRARY_DIR}")
      if(SKBUILD)
        # Wheel layout: _pnetcdf_commux.so sits in pinc/, libpnc in pinc/lib,
        # and commux/torch are sibling site-packages one level up.
        set(_pnc_py_install_rpath "$ORIGIN/lib;$ORIGIN;$ORIGIN/../commux/lib;$ORIGIN/../torch/lib")
      endif()
      set_target_properties(pnetcdf_commux_py PROPERTIES
        PREFIX ""
        OUTPUT_NAME "_pnetcdf_commux"
        BUILD_RPATH "$<TARGET_FILE_DIR:pnc>;${COMMUX_DIR}/lib;${TORCH_LIBRARY_DIR}"
        INSTALL_RPATH "${_pnc_py_install_rpath}")
      target_compile_definitions(pnetcdf_commux_py PRIVATE _GLIBCXX_USE_CXX11_ABI=${TORCH_CXX11_ABI})
      target_include_directories(pnetcdf_commux_py
        PRIVATE
          "${CMAKE_CURRENT_SOURCE_DIR}/src/include"
          "${TORCH_INSTALL_PREFIX}/include"
          "${TORCH_INSTALL_PREFIX}/include/torch/csrc/api/include")
      target_link_libraries(pnetcdf_commux_py PRIVATE pnc Python3::Module)
      install(TARGETS pnetcdf_commux_py LIBRARY DESTINATION "${PNETCDF_INSTALL_PYDIR}")
    endif()
  endif()
  if(BUILD_TESTING)
    add_executable(commux_smoke test/commux_smoke.c)
    target_link_libraries(commux_smoke PRIVATE pnetcdf)
    add_test(NAME commux_smoke COMMAND commux_smoke)
    add_executable(commux_c_4proc test/commux_c_4proc.c)
    target_link_libraries(commux_c_4proc PRIVATE pnc)
    add_executable(ncmpix_cdf5_4proc test/ncmpix_cdf5_4proc.c)
    target_link_libraries(ncmpix_cdf5_4proc PRIVATE pnc netcdf)
    add_executable(ncmpix_record_4proc test/ncmpix_record_4proc.c)
    target_link_libraries(ncmpix_record_4proc PRIVATE pnc netcdf)

    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -c "import shutil; print(shutil.which('torchrun') or '')"
      OUTPUT_VARIABLE TORCHRUN_EXECUTABLE
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
    if(TORCHRUN_EXECUTABLE)
      add_test(NAME commux_c_4proc
        COMMAND "${TORCHRUN_EXECUTABLE}" --nproc-per-node=4 --no-python
                "$<TARGET_FILE:commux_c_4proc>")
      add_test(NAME ncmpix_cdf5_4proc
        COMMAND "${TORCHRUN_EXECUTABLE}" --nproc-per-node=4 --no-python
                "$<TARGET_FILE:ncmpix_cdf5_4proc>")
      add_test(NAME ncmpix_record_4proc
        COMMAND "${TORCHRUN_EXECUTABLE}" --nproc-per-node=4 --no-python
                "$<TARGET_FILE:ncmpix_record_4proc>")
      add_test(NAME commux_4proc
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test/commux_4proc.py")
      set_tests_properties(commux_4proc PROPERTIES SKIP_RETURN_CODE 77)
    endif()
  endif()
  return()
endif()

find_package(MPI REQUIRED COMPONENTS C)
find_program(M4_EXECUTABLE m4 REQUIRED)

set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_C)
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_DIRS})
check_include_file(stdlib.h HAVE_STDLIB_H)
check_include_file(sys/mount.h HAVE_SYS_MOUNT_H)
check_include_file(sys/param.h HAVE_SYS_PARAM_H)
check_include_file(sys/stat.h HAVE_SYS_STAT_H)
check_include_file(sys/statvfs.h HAVE_SYS_STATVFS_H)
check_include_file(sys/types.h HAVE_SYS_TYPES_H)
check_include_file(sys/vfs.h HAVE_SYS_VFS_H)
check_include_file(unistd.h HAVE_UNISTD_H)
check_function_exists(access HAVE_ACCESS)
check_function_exists(close HAVE_CLOSE)
check_function_exists(lstat HAVE_LSTAT)
check_function_exists(open HAVE_OPEN)
check_function_exists(stat HAVE_STAT)
check_function_exists(statfs HAVE_STATFS)
check_function_exists(statvfs HAVE_STATVFS)
check_function_exists(strerror HAVE_STRERROR)
check_function_exists(unlink HAVE_UNLINK)
test_big_endian(WORDS_BIGENDIAN)
if(WORDS_BIGENDIAN)
  set(FLOAT_WORDS_BIGENDIAN 1)
endif()

check_c_source_compiles("#include <stddef.h>\nint main(void){ptrdiff_t x=0; return (int)x;}" HAVE_PTRDIFF_T)
check_c_source_compiles("#include <sys/types.h>\nint main(void){uint x=0; return (int)x;}" HAVE_UINT)
check_c_source_compiles("#include <sys/types.h>\nint main(void){ushort x=0; return (int)x;}" HAVE_USHORT)
check_c_source_compiles("#include <sys/types.h>\nint main(void){uchar x=0; return (int)x;}" HAVE_UCHAR)
check_c_source_compiles("int main(void){longlong x=0; return (int)x;}" HAVE_LONGLONG)
check_c_source_compiles("int main(void){ulonglong x=0; return (int)x;}" HAVE_ULONGLONG)
check_c_source_compiles("int main(void){int64 x=0; return (int)x;}" HAVE_INT64)
check_c_source_compiles("int main(void){uint64 x=0; return (int)x;}" HAVE_UINT64)

set(CMAKE_EXTRA_INCLUDE_FILES mpi.h)
check_type_size("MPI_Aint" SIZEOF_MPI_AINT LANGUAGE C)
check_type_size("MPI_Count" SIZEOF_MPI_COUNT LANGUAGE C)
check_type_size("MPI_Offset" SIZEOF_MPI_OFFSET LANGUAGE C)
set(CMAKE_EXTRA_INCLUDE_FILES)
check_type_size("char" SIZEOF_CHAR LANGUAGE C)
check_type_size("signed char" SIZEOF_SIGNED_CHAR LANGUAGE C)
check_type_size("unsigned char" SIZEOF_UNSIGNED_CHAR LANGUAGE C)
check_type_size("short" SIZEOF_SHORT LANGUAGE C)
check_type_size("unsigned short" SIZEOF_UNSIGNED_SHORT LANGUAGE C)
check_type_size("unsigned short int" SIZEOF_UNSIGNED_SHORT_INT LANGUAGE C)
check_type_size("int" SIZEOF_INT LANGUAGE C)
check_type_size("unsigned int" SIZEOF_UNSIGNED_INT LANGUAGE C)
check_type_size("long" SIZEOF_LONG LANGUAGE C)
check_type_size("float" SIZEOF_FLOAT LANGUAGE C)
check_type_size("double" SIZEOF_DOUBLE LANGUAGE C)
check_type_size("long long" SIZEOF_LONG_LONG LANGUAGE C)
check_type_size("unsigned long long" SIZEOF_UNSIGNED_LONG_LONG LANGUAGE C)
check_type_size("off_t" SIZEOF_OFF_T LANGUAGE C)
check_type_size("ptrdiff_t" SIZEOF_PTRDIFF_T LANGUAGE C)
check_type_size("size_t" SIZEOF_SIZE_T LANGUAGE C)

if(SIZEOF_MPI_OFFSET EQUAL SIZEOF_LONG_LONG)
  set(MPI_OFFSET_FMT "%lld")
elseif(SIZEOF_MPI_OFFSET EQUAL SIZEOF_LONG)
  set(MPI_OFFSET_FMT "%ld")
else()
  set(MPI_OFFSET_FMT "%d")
endif()

foreach(sym
    MPI_COMBINER_DARRAY
    MPI_COMBINER_DUP
    MPI_COMBINER_HINDEXED_BLOCK
    MPI_COMBINER_HINDEXED_INTEGER
    MPI_COMBINER_HVECTOR_INTEGER
    MPI_COMBINER_INDEXED_BLOCK
    MPI_COMBINER_STRUCT_INTEGER
    MPI_COMBINER_SUBARRAY)
  check_c_source_compiles("#include <mpi.h>\nint main(void){int x=${sym}; return x == 0;}" HAVE_${sym})
endforeach()

check_symbol_exists(MPI_OFFSET mpi.h HAVE_DECL_MPI_OFFSET)
check_symbol_exists(MPI_STATUSES_IGNORE mpi.h HAVE_MPI_STATUSES_IGNORE)
check_symbol_exists(MPI_Type_size_c mpi.h HAVE_MPI_TYPE_SIZE_C)
check_symbol_exists(MPI_Type_size_x mpi.h HAVE_MPI_TYPE_SIZE_X)
check_symbol_exists(MPI_Get_count_c mpi.h HAVE_MPI_GET_COUNT_C)
if(HAVE_MPI_TYPE_SIZE_C AND HAVE_MPI_GET_COUNT_C)
  set(HAVE_MPI_LARGE_COUNT 1)
endif()

file(MAKE_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/src/include"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/include"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/common"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio"
  "${CMAKE_CURRENT_BINARY_DIR}/src/dispatchers")

configure_file(src/include/pnetcdf.h.in "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf.h" @ONLY)
file(READ "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf.h" _pnetcdf_h)
string(REPLACE "DIST_DATE" "${PNETCDF_RELEASE_DATE}" _pnetcdf_h "${_pnetcdf_h}")
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf.h" "${_pnetcdf_h}")
configure_file(src/include/pnetcdf_comm.h.in "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf_comm.h" @ONLY)
configure_file(cmake/config.h.in "${CMAKE_CURRENT_BINARY_DIR}/src/include/config.h" @ONLY)

function(pnetcdf_m4 out in)
  string(REPLACE ";" " " _m4_args "${ARGN}")
  add_custom_command(
    OUTPUT "${out}"
    COMMAND /bin/sh -c "\"${M4_EXECUTABLE}\" -I \"${CMAKE_CURRENT_SOURCE_DIR}/m4\" ${_m4_args} \"${in}\" > \"${out}\""
    DEPENDS "${in}"
    VERBATIM)
endfunction()

pnetcdf_m4("${CMAKE_CURRENT_BINARY_DIR}/src/drivers/include/ncx.h"
           "${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/include/ncx_h.m4")

set(_m4_defs)
if(PNETCDF_ENABLE_ERANGE_FILL)
  list(APPEND _m4_defs -DERANGE_FILL)
endif()
foreach(name ncx convert_swap)
  pnetcdf_m4("${CMAKE_CURRENT_BINARY_DIR}/src/drivers/common/${name}.c"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/common/${name}.m4"
             ${_m4_defs})
endforeach()
foreach(name ncmpio_getput ncmpio_i_getput ncmpio_varn ncmpio_i_varn ncmpio_attr)
  pnetcdf_m4("${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio/${name}.c"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/ncmpio/${name}.m4"
             ${_m4_defs})
endforeach()
foreach(name attr_getput var_getput)
  pnetcdf_m4("${CMAKE_CURRENT_BINARY_DIR}/src/dispatchers/${name}.c"
             "${CMAKE_CURRENT_SOURCE_DIR}/src/dispatchers/${name}.m4")
endforeach()

set(COMM_SOURCES src/comm/pnc_comm_mpi.c)
if(PNETCDF_COMM_BACKEND STREQUAL "commux")
  set(COMM_SOURCES src/comm/pnc_comm_commux.c)
endif()

set(COMMON_SOURCES
  src/drivers/common/utf8proc.c
  src/drivers/common/mem_alloc.c
  src/drivers/common/dtype_decode.c
  src/drivers/common/create_imaptype.c
  src/drivers/common/check_name.c
  src/drivers/common/pack_unpack.c
  src/drivers/common/utils.c
  src/drivers/common/error_posix2nc.c
  src/drivers/common/error_mpi2nc.c
  src/drivers/common/hash_map.c
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/common/ncx.c"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/common/convert_swap.c")

set(NCMPIO_SOURCES
  src/drivers/ncmpio/ncmpio_driver.c
  src/drivers/ncmpio/ncmpio_create.c
  src/drivers/ncmpio/ncmpio_open.c
  src/drivers/ncmpio/ncmpio_enddef.c
  src/drivers/ncmpio/ncmpio_close.c
  src/drivers/ncmpio/ncmpio_sync.c
  src/drivers/ncmpio/ncmpio_file_misc.c
  src/drivers/ncmpio/ncmpio_header_get.c
  src/drivers/ncmpio/ncmpio_header_put.c
  src/drivers/ncmpio/ncmpio_dim.c
  src/drivers/ncmpio/ncmpio_var.c
  src/drivers/ncmpio/ncmpio_bput.c
  src/drivers/ncmpio/ncmpio_wait.c
  src/drivers/ncmpio/ncmpio_fill.c
  src/drivers/ncmpio/ncmpio_util.c
  src/drivers/ncmpio/ncmpio_hash_func.c
  src/drivers/ncmpio/ncmpio_file_io.c
  src/drivers/ncmpio/ncmpio_intra_node.c
  src/drivers/ncmpio/ncmpio_fstype.c
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio/ncmpio_getput.c"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio/ncmpio_i_getput.c"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio/ncmpio_varn.c"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio/ncmpio_i_varn.c"
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio/ncmpio_attr.c")
if(PNETCDF_ENABLE_SUBFILING)
  list(APPEND NCMPIO_SOURCES src/drivers/ncmpio/ncmpio_subfile.c)
endif()

set(DISPATCHER_SOURCES
  src/dispatchers/lib_version.c
  src/dispatchers/file.c
  src/dispatchers/dimension.c
  src/dispatchers/variable.c
  src/dispatchers/attribute.c
  src/dispatchers/error_codes.c
  src/dispatchers/cdl_header_parser.c
  "${CMAKE_CURRENT_BINARY_DIR}/src/dispatchers/attr_getput.c"
  "${CMAKE_CURRENT_BINARY_DIR}/src/dispatchers/var_getput.c")

add_library(pnetcdf
  ${COMM_SOURCES}
  ${COMMON_SOURCES}
  ${NCMPIO_SOURCES}
  ${DISPATCHER_SOURCES}
  "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/include/ncx.h")

target_compile_definitions(pnetcdf PRIVATE HAVE_CONFIG_H=1)
target_include_directories(pnetcdf
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/src/include>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/include"
    "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/drivers/ncmpio"
    "${CMAKE_CURRENT_BINARY_DIR}/src/drivers/ncmpio")
target_link_libraries(pnetcdf PUBLIC MPI::MPI_C)
set_target_properties(pnetcdf PROPERTIES
  OUTPUT_NAME pnetcdf
  VERSION 8.0.0
  SOVERSION 8)

install(TARGETS pnetcdf
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf.h"
  "${CMAKE_CURRENT_BINARY_DIR}/src/include/pnetcdf_comm.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
