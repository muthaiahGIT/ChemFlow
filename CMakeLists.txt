cmake_minimum_required(VERSION 3.2)

project(chemflow)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -O2 -DNoRepository -ftemplate-depth-100 -fPIC -Xlinker --add-needed -Xlinker --no-as-needed")
set(INC_EIGEN /usr/include)
set(INC_OF /opt/openfoam6/src)

add_subdirectory(src/ChemThermo)

include_directories(src ${INC_EIGEN}/eigen3 ${INC_OF}/OpenFOAM/lnInclude ${INC_OF}/thermophysicalModels/specie/lnInclude ${INC_OF}/thermophysicalModels/reactionThermo/lnInclude ${INC_OF}/thermophysicalModels/basic/lnInclude ${INC_OF}/transportModels/compressible/lnInclude ${INC_OF}/ODE/lnInclude ${INC_OF}/thermophysicalModels/chemistryModel/lnInclude ${INC_OF}/finiteVolume/lnInclude ${INC_OF}/meshTools/lnInclude ${INC_OF}/OpenFOAM/lnInclude /opt/openfoam6/src/OSspecific/POSIX/lnInclude)

link_directories(/opt/openfoam6/platforms/linux64GccDPInt32Opt/lib)

set(SOURCE_FILES src/main.cpp src/tdma.cpp)

add_executable(chemflow ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(chemflow ChemThermo chemistryModel finiteVolume meshTools OpenFOAM)
