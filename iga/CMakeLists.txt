message(STATUS "==> Processing iga directory")

set(HEADERS ShapeFunctions.h)
set(SOURCES)

include_directories(
${CMAKE_SOURCE_DIR}
)

set(EXTERNAL_LIBRARIES)

# add_library(iga ${HEADERS} ${SOURCES})
# target_link_libraries(iga ${EXTERNAL_LIBRARIES})

message(STATUS "yo, ${CMAKE_SOURCE_DIR}")