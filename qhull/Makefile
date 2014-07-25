CXX=g++-mp-4.8
#CXX=clang++
#CXX=~/usr/local/llvm/bin-3.4/bin/clang++
LD=$(CXX)

CXXFLAGS=-O3 -std=c++11 -g -Wall
LDFLAGS=

OBJ=qhull.o

PROG=qhull

all: $(PROG)

$(PROG): $(OBJ)
	$(LD) $(LDFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ): qhull.h linsolve.h vector.h


