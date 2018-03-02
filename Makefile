CXX = g++

# On MacOS, recent versions of gcc are aliased as (compiler_name)-7
DARWIN_GCC := $(shell command -v g++-7 2> /dev/null)
ifdef DARWIN_GCC
  CXX=g++-7
endif

LIBS = -lsfml-window -lsfml-graphics -lsfml-system

CXXFLAGS = -std=c++17 -O3 -Wall -flto -fopenmp -march=native -mtune=native

.PHONY: default all clean

default: all
all: raytracer

clean:
	-rm -f raytracer

raytracer: raytracer.cpp
	$(CXX) $(CXXFLAGS) raytracer.cpp $(LIBS) -o $@

grind: raytracer
	valgrind --leak-check=yes --read-var-info=yes --track-origins=yes ./raytracer
