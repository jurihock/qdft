# https://github.com/dpilger26/NumCpp

CPMAddPackage(
  NAME numcpp
  VERSION 2.8.0
  GIT_TAG Version_2.8.0
  GITHUB_REPOSITORY dpilger26/NumCpp
  DOWNLOAD_ONLY YES)

if(numcpp_ADDED)

  add_library(numcpp INTERFACE)

  target_include_directories(numcpp
    INTERFACE "${numcpp_SOURCE_DIR}/include")

  target_compile_features(numcpp
    INTERFACE cxx_std_14)

  target_compile_definitions(numcpp
    INTERFACE -DNUMCPP_NO_USE_BOOST)

  target_compile_definitions(numcpp
    INTERFACE -DNUMCPP_INCLUDE_PYBIND_PYTHON_INTERFACE)

endif()
