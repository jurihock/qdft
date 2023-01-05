file(GLOB EXAMPLES "${CMAKE_CURRENT_LIST_DIR}/*.cmake")

foreach(EXAMPLE ${EXAMPLES})

  include("${EXAMPLE}")

endforeach()
