cmake_minimum_required(VERSION 2.8.12)

set(TARGET_NAME mimalloc)
set(EXTENSION_NAME ${TARGET_NAME}_extension)
project(${TARGET_NAME})

include_directories(include)
include_directories(mimalloc/include)

add_subdirectory(mimalloc)
set(MIMALLOC_EXTENSION_FILES mimalloc_extension.cpp ${MIMALLOC_OBJECT_FILES})

set(PARAMETERS "-warnings")
build_static_extension(mimalloc ${MIMALLOC_EXTENSION_FILES})
# we do not do build_loadable_extension here because mimalloc is static-only

install(
  TARGETS mimalloc_extension
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
