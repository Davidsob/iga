set(library iga)

file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

include_directories(
${CMAKE_SOURCE_DIR}
${TPL_DIR}/eigen-3.3.5
${CMAKE_SOURCE_DIR}/bsplines/
${CMAKE_SOURCE_DIR}/bsplines/utils
)

set(EXTERNAL_LIBS)

message(STATUS "Processing iga directory...")
message(STATUS "library named: ${library}")
message(STATUS "library headers: ${HEADERS}")
message(STATUS "library sources: ${SOURCES}")

add_library(${library} ${HEADERS} ${SOURCES})
target_link_libraries(${library} ${EXTERNAL_LIBS})

