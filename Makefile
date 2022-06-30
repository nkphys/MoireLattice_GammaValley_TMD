# Platform: linux

EXENAME  = MoireBands
### ------ Personal PC compilation ------------
CXX     = g++
CPPFLAGS = -std=c++11
LDFLAGS  = -llapack -lblas

### ------ Newton compilation ------------
### MUST USE: module load gcc/4.8.2
#CXX = icpc  ### Or use g++ (both works!)
#CPPFLAGS = -std=c++11
#LDFLAGS = /data/apps/lapack/3.5.0/lib/liblapack.a /data/apps/lapack/3.5.0/lib/libblas.a -lgfortran

# LDFLAGS  = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -L$(INTELROOT)/lib/intel64 -lpthread
# LDFLAGS  = -L$(MKLROOT)/lib/intel64  -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64  -lm  -liomp5  -qopenmp -lpthread
#LDFLAGS  =  /opt/intel/mkl/lib/libmkl_core.a  /opt/intel/mkl/lib/libmkl_intel_lp64.a  /opt/intel/mkl/lib/libmkl_intel_thread.a -ldl -lpthread -lm -I/opt/intel/mkl/include -qopenmp -liomp5
### --- turn on for production -----------
#CPPFLAGS += -Isrc  #### Look inside the src folder for header files
#CPPFLAGS += -Wall -Werror -Wextra #### Enable warnings and treat warnings as errors
#CPPFLAGS += -DNDEBUG #### This disables debugging
#CPPFLAGS += -O3 #### Optimization level here
#STRIP_COMMAND = strip #### "strips off" all the debugging lines of executable

### --- turn on for debugging -----------
CPPFLAGS += -Isrc
#CPPFLAGS += -Wall -Werror -Wextra #### This enables warnings with extra debugging
#CPPFLAGS += -g3 #### link gdb to file system of program
CPPFLAGS += -O3 #### Reduce compilation time and make debugging produce the expected results.
STRIP_COMMAND = true #### Keeps lines in the executable for debugging

$(EXENAME): clean main.o 
	$(CXX) $(CPPFLAGS) -o $(EXENAME)  main.o $(LDFLAGS) 
	$(STRIP_COMMAND) $(EXENAME)

all: $(EXENAME)

	 
clean:
	rm -f $(EXENAME) MoireBands *.o

######## End of Makefile ########
