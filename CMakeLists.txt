project(bsplines)
cmake_minimum_required(VERSION 2.8)

# subdirectories
add_subdirectory(splines) # get down in the weeds and build/link other dirs
add_subdirectory(iga) # get down in the weeds and build/link other dirs

include_directories(
${CMAKE_SOURCE_DIR}
${TPL_DIR}/eigen-3.3.5
${CMAKE_SOURCE_DIR}/iga
)

message(STATUS "Include Directories:")
message(STATUS "${CMAKE_SOURCE_DIR}")
message(STATUS "${TPL_DIR}/eigen-3.3.5")

link_directories()

set(EXTERNAL_LIBS iga splines)

# set the executable
add_executable(Main.exe main.cxx)
target_link_libraries(Main.exe ${EXTERNAL_LIBS})

add_executable(Elbow.exe elbow.cxx)
target_link_libraries(Elbow.exe ${EXTERNAL_LIBS})

add_executable(Heat.exe heat_conduction.cxx)
target_link_libraries(Heat.exe ${EXTERNAL_LIBS})

add_executable(Stress.exe solid_mechanics.cxx)
target_link_libraries(Stress.exe ${EXTERNAL_LIBS})

add_executable(PlaneStrain.exe plane_strain_ex.cxx)
target_link_libraries(PlaneStrain.exe ${EXTERNAL_LIBS})

add_executable(BBar.exe bbar_ex.cxx)
target_link_libraries(BBar.exe ${EXTERNAL_LIBS})