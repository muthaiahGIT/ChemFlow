cmake_minimum_required(VERSION 3.2)

project(chemflow)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -DNoRepository -ftemplate-depth-100 -fPIC -Xlinker --add-needed -Xlinker --no-as-needed")

add_subdirectory(src/ChemThermo)

include_directories(src /usr/include/eigen3 /opt/openfoam6/src/OpenFOAM/lnInclude /opt/openfoam6/src/thermophysicalModels/specie/lnInclude /opt/openfoam6/src/thermophysicalModels/reactionThermo/lnInclude /opt/openfoam6/src/thermophysicalModels/basic/lnInclude /opt/openfoam6/src/transportModels/compressible/lnInclude /opt/openfoam6/src/ODE/lnInclude /opt/openfoam6/src/thermophysicalModels/chemistryModel/lnInclude /opt/openfoam6/src/finiteVolume/lnInclude /opt/openfoam6/src/meshTools/lnInclude /opt/openfoam6/src/OpenFOAM/lnInclude /opt/openfoam6/src/OSspecific/POSIX/lnInclude)

link_directories(/opt/openfoam6/platforms/linux64GccDPInt32Opt/lib)

set(SOURCE_FILES src/main.cpp src/tdma.cpp)

add_executable(chemflow ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(chemflow ChemThermo chemistryModel finiteVolume meshTools OpenFOAM dl)
