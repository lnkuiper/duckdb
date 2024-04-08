cmake_minimum_required(VERSION 2.8.12)

set(TARGET_NAME mimalloc)
set(EXTENSION_NAME ${TARGET_NAME}_extension)
project(${TARGET_NAME})

include_directories(include)
set(MIMALLOC_EXTENSION_FILES mimalloc_extension.cpp ${MIMALLOC_OBJECT_FILES})

set(PARAMETERS "-warnings")
build_static_extension(mimalloc ${MIMALLOC_EXTENSION_FILES})
# we do not do build_loadable_extension here because mimalloc is static-only

find_package(mimalloc CONFIG)

if(NOT ${mimalloc_FOUND})
  message(FATAL_ERROR "mimalloc not found, did you set up vcpkg correctly?")
endif()

# Static lib
target_link_libraries(${EXTENSION_NAME} mimalloc-static)
target_include_directories(${EXTENSION_NAME} PRIVATE ${MIMALLOC_DIR}/include)

install(
  TARGETS mimalloc_extension
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
