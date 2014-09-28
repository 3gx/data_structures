CXX=g++-mp-4.8
CXX=g++
#CXX=clang++
#CXX=~/usr/local/llvm/bin-3.4/bin/clang++
LD=$(CXX)

CXXFLAGS=-O2 -std=c++11 -g -Wall
LDFLAGS=

OBJ=qhull.o

PROG=qhull

all: $(PROG)

$(PROG): $(OBJ)
	$(LD) $(LDFLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	/bin/rm -rf $(OBJ) $(PROG)

$(OBJ): qhull.h linsolve.h vector.h


