# https://github.com/mackron/dr_libs

CPMAddPackage(
  NAME dr
  VERSION 0.13.7
  GIT_TAG dd762b861ecadf5ddd5fb03e9ca1db6707b54fbb
  GITHUB_REPOSITORY mackron/dr_libs
  DOWNLOAD_ONLY YES)

if(dr_ADDED)

  add_library(dr INTERFACE)

  target_include_directories(dr
    INTERFACE "${dr_SOURCE_DIR}")

  target_compile_definitions(dr
    INTERFACE -DDR_WAV_IMPLEMENTATION)

endif()
