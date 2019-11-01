project(bsplines)
cmake_minimum_required(VERSION 2.8)

# include libraries
include_directories(
)

# subdirectories
#add_subdirectory(math)

link_directories(
)


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
