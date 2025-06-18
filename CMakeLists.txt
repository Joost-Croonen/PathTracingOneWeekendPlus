# CMakeList.txt : CMake project for TestSetup, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("PathTracingOneWeekendPlus")


find_package(OpenMP REQUIRED)

# Add source to this project's executable.
add_executable (PathTracingOneWeekendPlus   "main.cpp" "vec3.h" "color.h" "ray.h" "hittable.h" "sphere.h" "hittable_list.h" "rtweekend.h" "interval.h" "camera.h" "material.h" "aabb.h" "bvh.h" "texture.h" "rtw_stb_image.h" "perlin.h" "quad.h" "onb.h" "pdf.h")

target_link_libraries(PathTracingOneWeekendPlus PRIVATE OpenMP::OpenMP_CXX)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET PathTracingOneWeekendPlus PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
