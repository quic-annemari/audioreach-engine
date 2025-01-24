#[[
   @file version.cmake

   @brief

   @copyright
   Copyright (c) Qualcomm Innovation Center, Inc. All Rights Reserved.
   SPDX-License-Identifier: BSD-3-Clause

]]
function(get_version)

  message(STATUS "${CMAKE_CXX_COMPILER_VERSION}")

  set(GCC_VERSION ${CMAKE_CXX_COMPILER_VERSION})

  string(REPLACE "." ";" VERSION_LIST ${GCC_VERSION})

  if(${ARSPF_PLATFORM_QNX}) #TODO
	  list(GET VERSION_LIST 0 MAJOR)
	  list(GET VERSION_LIST 1 MINOR)
	  list(GET VERSION_LIST 2 MICRO)
  else()
	  list(GET VERSION_LIST 0 MAJOR)
	  list(GET VERSION_LIST 1 MINOR)
	  list(GET VERSION_LIST 2 MICRO)
  endif()


  math(EXPR MAJOR_VERSION "${MAJOR}" OUTPUT_FORMAT HEXADECIMAL)
  string(REPLACE "0x" "" GCC_MAJOR_VERSION ${MAJOR_VERSION})
  string(LENGTH ${GCC_MAJOR_VERSION} MAJOR_LEN)
  if(${MAJOR_LEN} EQUAL 1)
      string(CONCAT GCC_MAJOR_VERSION "0" ${GCC_MAJOR_VERSION}) #Pre fixing 0
  endif()

  math(EXPR MINOR_VERSION "${MINOR}" OUTPUT_FORMAT HEXADECIMAL)
  string(REPLACE "0x" "" GCC_MINOR_VERSION ${MINOR_VERSION})
  string(LENGTH ${GCC_MINOR_VERSION} MINOR_LEN)
  if(${MINOR_LEN} EQUAL 1)
      string(CONCAT GCC_MINOR_VERSION "0" ${GCC_MINOR_VERSION}) #Pre fixing 0
  endif()

  math(EXPR MICRO_VERSION "${MICRO}" OUTPUT_FORMAT HEXADECIMAL)
  string(REPLACE "0x" "" GCC_MICRO_VERSION ${MICRO_VERSION})
  string(LENGTH ${GCC_MICRO_VERSION} MICRO_LEN)
  if(${MICRO_LEN} EQUAL 1)
      string(CONCAT GCC_MICRO_VERSION "0" ${GCC_MICRO_VERSION}) #Pre fixing 0
  endif()
  STRING(JOIN "" CURRENT_GCC_VERSION "0x00" ${GCC_MAJOR_VERSION} ${GCC_MINOR_VERSION} ${GCC_MICRO_VERSION})
  set(CURRENT_GCC_VERSION ${CURRENT_GCC_VERSION} PARENT_SCOPE)
  message(STATUS "Converted GCC version: ${CURRENT_GCC_VERSION}")
endfunction()
get_version()
