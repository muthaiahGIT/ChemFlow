cc = g++

INC_EIGEN = /usr/include
INC_OF = /opt/openfoam6/src
LIB_OF = /opt/openfoam6/platforms/linux64GccDPInt32Opt/lib
LIBS = -lchemistryModel -lfiniteVolume -lmeshTools -lOpenFOAM -lpthread
SRC = ./src
INC_ALL = -I $(SRC) -I $(INC_EIGEN)/eigen3 -I $(INC_OF)/OpenFOAM/lnInclude \
-I $(INC_OF)/thermophysicalModels/specie/lnInclude -I $(INC_OF)/thermophysicalModels/reactionThermo/lnInclude \
-I $(INC_OF)/thermophysicalModels/basic/lnInclude -I $(INC_OF)/transportModels/compressible/lnInclude \
-I $(INC_OF)/ODE/lnInclude -I $(INC_OF)/thermophysicalModels/chemistryModel/lnInclude \
-I $(INC_OF)/finiteVolume/lnInclude -I $(INC_OF)/meshTools/lnInclude -I $(INC_OF)/OpenFOAM/lnInclude \
-I $(INC_OF)/OSspecific/POSIX/lnInclude
CFLAGS = -m64 -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -Wall -Wextra -Wnon-virtual-dtor \
-Wno-unused-parameter -Wno-invalid-offsetof -Wno-attributes -O2 -DNoRepository -ftemplate-depth-200 \
-fPIC -Xlinker --add-needed -Xlinker --no-as-needed


chemflow: $(SRC)/*.cpp $(SRC)/*.h
	$(cc) $(CFLAGS) -o chemflow $(SRC)/*.cpp $(INC_ALL) -L $(LIB_OF) $(LIBS)


clean:
	rm -f *.o chemflow 
