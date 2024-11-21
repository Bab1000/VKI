MPP_INSTALL_DIR =/Users/Bruninho/Codes/gsi_solvers/mutation++_git/install


CPP = g++
CPP_FLAGS = -g -O2 -w -I $(MPP_INSTALL_DIR)/include/mutation++
CPP_LIBS  = $(MPP_INSTALL_DIR)/lib/libmutation++.dylib

post_process : post_process.o StagLineSol.o
	$(CPP) -o $@ $(CPP_FLAGS)  $^ $(CPP_LIBS)

%.o : %.cpp
	$(CPP) -c $(CPP_FLAGS) $<

clean :
	rm *.o post_process
