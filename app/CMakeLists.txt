# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.20.0)

find_package(Zephyr REQUIRED HINTS $ENV{ZEPHYR_BASE})
project(sample_are)

add_compile_options(-Wno-dev)

target_sources(app PRIVATE src/main.c)

target_compile_options(app PRIVATE -g)

target_link_libraries(app PUBLIC spf)

# Added to check for alignment error in filesz in zephyr.elf fragments. Do readelf -l zephyr.elf
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--align-segments=4")

include_directories(
	${ZEPHYR_ARE_MODULE_DIR}/fwk/spf/utils/cmn/inc
	${ZEPHYR_ARE_MODULE_DIR}/fwk/platform/posal/inc
	${ZEPHYR_ARE_MODULE_DIR}/fwk/platform/posal/inc/linux/
	${ZEPHYR_ARE_MODULE_DIR}/fwk/platform/posal/inc/linux/stringl
	${ZEPHYR_ARE_MODULE_DIR}/fwk/api/ar_utils
	${ZEPHYR_ARE_MODULE_DIR}/fwk/api/ar_utils/generic/
	${ZEPHYR_ARE_MODULE_DIR}/ar_osal/api
	${ZEPHYR_ARE_MODULE_DIR}/gpr/api
	${ZEPHYR_HAL_XTENSA_MODULE_DIR}/include
	${ZEPHYR_BASE}/include/zephyr/
	${ZEPHYR_BASE}/include/
	#${ZEPHYR_BASE}/include/zephyr/posix
	${CMAKE_BINARY_DIR}/zephyr/include/generated
  	${CMAKE_BINARY_DIR}/include/generated
	)

zephyr_include_directories()