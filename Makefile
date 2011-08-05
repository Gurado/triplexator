
# Makefile for linux 
# run with make 

# indicate location of BOOST include files
BOOST=
#BOOST=-I/usr/include/boost/ -DBOOST=1

# indicate location of BOOST libraries
BOOSTLIB=
#BOOSTLIB=/usr/lib/libz.dylib /usr/lib/libbz2.dylib /usr/lib/libboost_iostreams.so

# whether to use openMP
OPENMP=
# comment the next line if OpenMP is not supported by the compiler
OPENMP=-fopenmp -DSEQAN_ENABLE_PARALLELISM=1


# Compiler of choice (tested with gcc v4.4 on Mac osx and Linux x86)
CXX=g++

CXXFLAGS+=-Ilibrary -static
CXXFLAGS+=$(OPENMP) $(BOOST) 
CXXFLAGS+=-O3 -DNDEBUG -DSEQAN_ENABLE_TESTING=0 -DSEQAN_ENABLE_DEBUG=0
LDFLAGS+=$(OPENMP) $(BOOSTLIB) -static 

SRC=src

default: all

all: main test

main: triplexator.o
	$(CXX) $(LDFLAGS) -o bin/triplexator $(SRC)/triplexator.o

triplexator.o: $(SRC)/triplexator.cpp
	$(CXX) $(CXXFLAGS) -c -o $(SRC)/triplexator.o $(SRC)/triplexator.cpp

test: main
	./demos/smoketest_triplexator.sh ./bin/triplexator
	
clean:
	rm -f triplexator.o triplexator 
	rm -r -f ./demos/examples
	rm -r -f ./demos/tests

.PHONY: default all  clean

