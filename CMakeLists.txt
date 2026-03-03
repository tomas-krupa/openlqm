# Copyright 2025, Noblis, Inc.
# https://fingerprint.nist.gov/openlqm

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at https://www.apache.org/licenses/LICENSE-2.0.

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# OpenLQM was developed by Noblis, Inc. under contract to the National
# Institute of Standards and Technology in 2024-2025, as a port of "LQMetric,"
# which was developed by Noblis, Inc. under contract to the Federal Bureau of
# Investigation's Criminal Justice Information Services Division in 2012-2014.

# CMake deprecated the OLD behavior of CMP0000 which allowed cmake_minimum_required to be included from another .cmake file with only a
# a warning. Addressing this by setting a minimum in all CMakeLists.txt as required, but still including the "real" call from a common
# initialization script, which will be greater than or equal to the fake minimum used here.
cmake_minimum_required(VERSION 3.30)
project("SUPERBUILD")
include(OpenLQM/openlqm_cmake_utils/init.cmake)

option(STATICALLY_LINK_LIBRARY "Link C++ runtime statically to the OpenLQM library" OFF)
option(STATICALLY_LINK_EXECUTABLE "Link C++ runtime statically to the OpenLQM executable" OFF)

option(PROFILE_RELEASE_BUILD OFF)

# Default install location
# if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#todo: cmake with the ninja generator can't handle spaces in paths on linux, so overriding this for now
if(NOT INSTALL_PREFIX_OVERRIDE)
	get_filename_component(output_directory_name "${CMAKE_BINARY_DIR}" NAME)
	set(INSTALL_PREFIX_OVERRIDE ${PROJECT_SOURCE_DIR}/out/install/${output_directory_name} CACHE PATH "CMake with the Ninja generator has an issue on linux with project paths containing spaces (no issue on Windows), resulting in a malformed CMAKE_INSTALL_PREFIX being set automatically. This override is a temporary workaround.")
endif()
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX_OVERRIDE} CACHE PATH "Install directory used by install(). OpenLQM defaults this to a project-local directory instead of the usual system directory" FORCE)

if(STATICALLY_LINK_LIBRARY)
	set(OpenCV_STATIC ON)
endif()

if(STATICALLY_LINK_LIBRARY)
	if(MSVC)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
		set(MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
		set(FORWARD_CMAKE_MSVC_RUNTIME_LIBRARY -DCMAKE_MSVC_RUNTIME_LIBRARY=${MSVC_RUNTIME_LIBRARY})
		set(FORWARD_MSVC_RUNTIME_LIBRARY -DMSVC_RUNTIME_LIBRARY=${MSVC_RUNTIME_LIBRARY})
	endif()
endif()

set(STATICALLY_LINK_OPENCV_CRT OFF)

if(MSVC)
	set(STATICALLY_LINK_OPENCV_CRT ${STATICALLY_LINK_LIBRARY})
endif()
set(OPENCV_STATIC_CRT_FLAG -DBUILD_WITH_STATIC_CRT=${STATICALLY_LINK_OPENCV_CRT})

set(OPENCV_REQUIRED_COMPONENTS core ml imgproc imgcodecs)
list(JOIN OPENCV_REQUIRED_COMPONENTS "," OPENCV_REQUIRED_COMPONENTS_DELIMITED)
foreach(component ${OPENCV_REQUIRED_COMPONENTS})
	list(APPEND OPENCV_REQUIRED_COMPONENT_TARGETS "opencv_${component}")
endforeach()
list(JOIN OPENCV_REQUIRED_COMPONENT_TARGETS "," OPENCV_REQUIRED_COMPONENT_TARGETS_DELIMITED)

add_subdirectory(OpenLQM)
