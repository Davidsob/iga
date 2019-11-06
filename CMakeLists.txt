project(bsplines)
cmake_minimum_required(VERSION 2.8)

# subdirectories
add_subdirectory(splines) # get down in the weeds and build/link other dirs
add_subdirectory(iga) # get down in the weeds and build/link other dirs

include_directories(
${CMAKE_SOURCE_DIR}
${TPL_DIR}/eigen-3.3.5
)

link_directories(
)

message(STATUS "${TPL_DIR}/eigen-3.3.5")

set(EXTERNAL_LIBS
#"/usr/lib64/libblas.so.3"
#"/usr/lib64/liblapack.so.3"
#utils
#grids
#materials
)




# set the executable
add_executable(Main.exe main.cxx)
target_link_libraries(Main.exe ${EXTERNAL_LIBS})

add_executable(Elbow.exe elbow.cxx)
target_link_libraries(Elbow.exe ${EXTERNAL_LIBS})
