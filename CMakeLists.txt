cmake_minimum_required(VERSION 3.0.0)
project(downbrust VERSION 0.1.0)

include(CTest)
enable_testing()

set(SRC src)
include_directories(${SRC})
aux_source_directory(${PROJECT_SOURCE_DIR}/src DIR_SRCS)
add_executable(${PROJECT_NAME} ${DIR_SRCS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
