aux_source_directory(. DIR_LIB_SRCS)

set(CMAKE_CXX_FLAGS "-m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -O2 -DNoRepository -ftemplate-depth-100 -fPIC")

include_directories(${INC_EIGEN}/eigen3 ${INC_OF}/OpenFOAM/lnInclude ${INC_OF}/thermophysicalModels/specie/lnInclude ${INC_OF}/thermophysicalModels/reactionThermo/lnInclude ${INC_OF}/thermophysicalModels/basic/lnInclude ${INC_OF}/transportModels/compressible/lnInclude ${INC_OF}/ODE/lnInclude ${INC_OF}/thermophysicalModels/chemistryModel/lnInclude ${INC_OF}/finiteVolume/lnInclude ${INC_OF}/meshTools/lnInclude ${INC_OF}/OpenFOAM/lnInclude /opt/openfoam6/src/OSspecific/POSIX/lnInclude)

add_library(ChemThermo ${DIR_LIB_SRCS})
