set(library splines)

file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

include_directories(
${CMAKE_SOURCE_DIR}
)

set(EXTERNAL_LIBS)

message(STATUS "Processing splines library...")
message(STATUS "library named: ${library}")
message(STATUS "library headers: ${HEADERS}")
message(STATUS "library sources: ${SOURCES}")

add_library(${library} ${HEADERS} ${SOURCES})
target_link_libraries(${library} ${EXTERNAL_LIBS})
