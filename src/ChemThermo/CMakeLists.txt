aux_source_directory(. DIR_LIB_SRCS)

set(CMAKE_CXX_FLAGS "-m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wnon-virtual-dtor -Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -DNoRepository -ftemplate-depth-100 -fPIC") 

include_directories(/usr/include/eigen3 /opt/openfoam6/src/thermophysicalModels/specie/lnInclude /opt/openfoam6/src/thermophysicalModels/reactionThermo/lnInclude /opt/openfoam6/src/thermophysicalModels/basic/lnInclude /opt/openfoam6/src/transportModels/compressible/lnInclude /opt/openfoam6/src/ODE/lnInclude /opt/openfoam6/src/thermophysicalModels/chemistryModel/lnInclude /opt/openfoam6/src/finiteVolume/lnInclude /opt/openfoam6/src/meshTools/lnInclude /opt/openfoam6/src/OpenFOAM/lnInclude /opt/openfoam6/src/OSspecific/POSIX/lnInclude) 

add_library(ChemThermo ${DIR_LIB_SRCS})
