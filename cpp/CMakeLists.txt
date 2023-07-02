cmake_minimum_required(VERSION 3.12...3.18)

project(qdft)

if(MSVC)
  # add_compile_options(/W3 /WX)
else()
  # add_compile_options(-Wall -Werror)
endif()

include("${CMAKE_CURRENT_LIST_DIR}/src/qdft/CMakeLists.txt")

include("${CMAKE_CURRENT_LIST_DIR}/lib/cpm.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/lib/python.cmake")

if(Python3_FOUND AND Python3_NumPy_FOUND)

  include("${CMAKE_CURRENT_LIST_DIR}/lib/matplotlibcpp.cmake")
  include("${CMAKE_CURRENT_LIST_DIR}/lib/numcpp.cmake")
  include("${CMAKE_CURRENT_LIST_DIR}/lib/pybind.cmake")

endif()

include("${CMAKE_CURRENT_LIST_DIR}/examples/CMakeLists.txt")
