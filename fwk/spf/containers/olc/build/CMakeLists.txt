#[[
   @file CMakeLists.txt

   @brief

   @copyright
   Copyright (c) Qualcomm Innovation Center, Inc. All Rights Reserved.
   SPDX-License-Identifier: BSD-3-Clause-Clear

]]
cmake_minimum_required(VERSION 3.10)

#Include directories
include_directories(
                    ../cmn/inc
                    ../core/src
                   )

#Add the sub directory
add_subdirectory(../core/build core)
