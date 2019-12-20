project(bsplines)
cmake_minimum_required(VERSION 2.8)

# subdirectories
add_subdirectory(splines) # get down in the weeds and build/link other dirs
add_subdirectory(iga) # get down in the weeds and build/link other dirs
add_subdirectory(utils) # get down in the weeds and build/link other dirs
add_subdirectory(thermal) # get down in the weeds and build/link other dirs
add_subdirectory(shellmechanics) # get down in the weeds and build/link other dirs
add_subdirectory(solidmechanics) # get down in the weeds and build/link other dirs

include_directories(
${CMAKE_SOURCE_DIR}
${TPL_DIR}/eigen-3.3.5
${CMAKE_SOURCE_DIR}/base
${CMAKE_SOURCE_DIR}/iga
${CMAKE_SOURCE_DIR}/utils
)

message(STATUS "Include Directories:")
message(STATUS "${CMAKE_SOURCE_DIR}")
message(STATUS "${TPL_DIR}/eigen-3.3.5")

link_directories()

set(EXTERNAL_LIBS iga splines utils)

# set the executable
# add_executable(Main.exe main.cxx)
# target_link_libraries(Main.exe ${EXTERNAL_LIBS})

# add_executable(Elbow.exe elbow.cxx)
# target_link_libraries(Elbow.exe ${EXTERNAL_LIBS})

# add_executable(Heat.exe heat_conduction.cxx)
# target_link_libraries(Heat.exe ${EXTERNAL_LIBS})

# add_executable(Stress.exe solid_mechanics.cxx)
# target_link_libraries(Stress.exe ${EXTERNAL_LIBS})

# add_executable(PlaneStrain.exe plane_strain_ex.cxx)
# target_link_libraries(PlaneStrain.exe ${EXTERNAL_LIBS})

# add_executable(BBar.exe bbar_ex.cxx)
# target_link_libraries(BBar.exe ${EXTERNAL_LIBS})

# add_executable(Integration.exe integration_ex.cxx)
# target_link_libraries(Integration.exe ${EXTERNAL_LIBS})

# add_executable(TimBeam.exe tim_beam_ex.cxx)
# target_link_libraries(TimBeam.exe ${EXTERNAL_LIBS})

# add_executable(AxialTimBeam.exe tim_beam_axial_ex.cxx)
# target_link_libraries(AxialTimBeam.exe ${EXTERNAL_LIBS})

# add_executable(KLShellEx.exe kl_shell_ex.cxx)
# target_link_libraries(KLShellEx.exe ${EXTERNAL_LIBS})

# add_executable(MRShellEx.exe mr_shell_ex.cxx)
# target_link_libraries(MRShellEx.exe ${EXTERNAL_LIBS})

# add_executable(SixDofShellEx.exe sixdof_shell_ex.cxx)
# target_link_libraries(SixDofShellEx.exe ${EXTERNAL_LIBS})

# add_executable(PinchedCylinder.exe pinched_shell_ex.cxx)
# target_link_libraries(PinchedCylinder.exe ${EXTERNAL_LIBS})

# add_executable(CurvedShellBeam.exe curved_shell_ex.cxx)
# target_link_libraries(CurvedShellBeam.exe ${EXTERNAL_LIBS})

# add_executable(Ht3D.exe ht3d_ex.cxx)
# target_link_libraries(Ht3D.exe ${EXTERNAL_LIBS})

# add_executable(Sm3D.exe sm3d_ex.cxx)
# target_link_libraries(Sm3D.exe ${EXTERNAL_LIBS})

# add_executable(NlSm3D.exe nlsm3d_ex.cxx)
# target_link_libraries(NlSm3D.exe ${EXTERNAL_LIBS})

add_executable(ShellEx1.exe shell_plane_ex.cxx)
target_link_libraries(ShellEx1.exe ${EXTERNAL_LIBS})

add_executable(ShellEx2.exe shell_pinched_cylinder_ex.cxx)
target_link_libraries(ShellEx2.exe ${EXTERNAL_LIBS})
