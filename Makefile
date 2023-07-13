# Arch Linux does not tag GCC versions.
# Builds on Brad's Laptop
CXX := g++-10

DIR := src/include # Include directory
CXXFLAGS := -std=c++2a -Wall -fno-strict-aliasing -fno-omit-frame-pointer
CPPFLAGS := -DDIM=2 -I $(DIR)

ifdef PROD
CPPFLAGS := -DNDEBUG $(CPPFLAGS)
CXXFLAGS := $(CXXFLAGS) -O3
else ifdef EXP
CPPFLAGS := -DNDEBUG -DSTAT $(CPPFLAGS)
CXXFLAGS := -ggdb3 $(CXXFLAGS) -O3
else ifdef DEBUG
CXXFLAGS := -ggdb3 $(CXXFLAGS) -O0
else
CXXFLAGS := -ggdb3 $(CXXFLAGS)
endif

SRC := $(shell find . -path ./src/tests -prune -o \( -name '*.cpp' -a ! -name 'pencilPrinter.cpp' \) -print)
MAINS := ./src/main.cpp ./src/gen_tree.cpp
TREE_NODES := ./src/nirtree/node.o ./src/rtree/node.o ./src/rplustree/node.o ./src/rstartree/node.o ./src/quadtree/node.o ./src/revisedrstartree/node.o ./src/bulk_load.o
SRC := $(filter-out $(MAINS),$(SRC))
OBJ := $(SRC:.cpp=.o)
OBJ_TO_COPY := $(filter-out $(TREE_NODES),$(OBJ))


TESTSRC := $(shell find ./src/tests -name '*.cpp')
TESTOBJ := $(TESTSRC:.cpp=.o)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Disable building tests for now
#all: bin/main bin/gen_tree bin/tests
all: bin/main bin/gen_tree

src/main.o : src/main.cpp src/include/bench/randomPoints.h src/include/globals/globals.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

src/gen_tree.o : src/gen_tree.cpp src/bulk_load.o src/include/bench/randomPoints.h src/include/globals/globals.h
	cp src/bulk_load.o bin/bulk_load.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

bin/main: $(OBJ) src/main.o
	mkdir -p bin
	cp src/nirtree/node.o bin/nirtreenode.o
	cp src/rtree/node.o bin/rtreenode.o
	cp src/rplustree/node.o bin/rplustreenode.o
	cp src/rstartree/node.o bin/rstartreenode.o
	cp src/quadtree/node.o bin/quadtreenode.o
	cp src/revisedrstartree/node.o bin/revisedrstartreenode.o
	rm -rf test*.o
	cp $(OBJ_TO_COPY) bin/
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) bin/*.o src/main.o -o bin/main -I $(DIR)

bin/gen_tree: $(OBJ) bin/main src/gen_tree.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) bin/*.o src/gen_tree.o -o bin/gen_tree -I $(DIR)

bin/tests: $(TESTOBJ) bin/main
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) bin/*.o src/tests/*.o -o bin/tests

.PHONY: all clean prod

clean:
	rm -rf bin/*
	find . -name "*.o" -delete
	find . -name "*.d" -delete
