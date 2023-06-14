# https://github.com/pybind/pybind11

CPMAddPackage(
  NAME pybind
  VERSION 2.10.3
  GITHUB_REPOSITORY pybind/pybind11
  DOWNLOAD_ONLY YES)

if(pybind_ADDED)

  add_library(pybind INTERFACE)

  target_include_directories(pybind
    INTERFACE "${pybind_SOURCE_DIR}/include")

  target_compile_features(pybind
    INTERFACE cxx_std_11)

endif()
