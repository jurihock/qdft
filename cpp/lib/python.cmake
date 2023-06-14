add_library(python INTERFACE)

find_package(Python3 COMPONENTS Development NumPy)

if(Python3_FOUND)

  target_include_directories(python
    INTERFACE "${Python3_INCLUDE_DIRS}")

  target_link_directories(python
    INTERFACE "${Python3_LIBRARY_DIRS}")

  target_link_libraries(python
    INTERFACE "${Python3_LIBRARIES}")

endif()

if(Python3_NumPy_FOUND)

  target_include_directories(python
    INTERFACE "${Python3_NumPy_INCLUDE_DIRS}")

endif()
