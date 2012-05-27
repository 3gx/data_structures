CXX = g++
CC  = gcc
LD  = g++
F90 = gfortran

ifeq ($(OMP), 1)
	OMPFLAGS  = -fopenmp
endif

ifeq ($(OMPCXX), 1)
	OMPFLAGS += -D_GLIBCXX_PARALLEL
endif

OFLAGS = -O4 -g -Wall
OFLAGS += -funroll-all-loops

CXXFLAGS =  -fPIC $(OFLAGS) -Wstrict-aliasing=2 $(OMPFLAGS)

ifeq ($(SSE), 1)
	CXXFLAGS += -D__mySSE__
endif



LDFLAGS = 

SRCPATH = ./
SRC = octree.cpp
OBJ = $(SRC:%.cpp=%.o) 

PROG  = octree

RM = /bin/rm

all:	  $(PROG) $(PROG1) $(PROG2)


$(PROG): $(OBJ) 
	$(LD) $(LDFLAGS) $^ -o $@ $(OMPFLAGS)

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	/bin/rm -rf *.o $(PROG) 

clean_all:
	/bin/rm -rf *.o $(PROG) plummer.dat *~

$(OBJ): octree.h plummer.h mytimer.h boundary.h vector3.h







