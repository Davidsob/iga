file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cpp")

include_directories(
${PROJECT_SOURCE_DIR}
${CMAKE_SOURCE_DIR}
${TPL_DIR}/eigen-3.3.5
)

set(EXTERNAL_LIBS)

add_library(utils ${HEADERS} ${SOURCES})
target_link_libraries(utils ${EXTERNAL_LIBS})
