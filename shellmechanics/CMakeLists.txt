set(library shellmechanics)

file(GLOB
	HEADERS
	"*.h"
	"boundary_conditions/*.h"
	"operators/*.h"
	"variables/*.h"
)

file(GLOB
	SOURCES	
	"*.cpp"
	"boundary_conditions/*.cpp"
	"operators/*.cpp"
	"variables/*.cpp"
)

include_directories(
${CMAKE_SOURCE_DIR}
${TPL_DIR}/eigen-3.3.5
#${CMAKE_SOURCE_DIR}/base
)

set(EXTERNAL_LIBS splines)

message(STATUS "Processing iga library...")
message(STATUS "library named: ${library}")
message(STATUS "library headers: ${HEADERS}")
message(STATUS "library sources: ${SOURCES}")

add_library(${library} ${HEADERS} ${SOURCES})
target_link_libraries(${library} ${EXTERNAL_LIBS})


