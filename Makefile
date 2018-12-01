CC = g++-7
OUT_FILE = a.out

DIRS := .
SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))

uname_S=$(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(uname_S),Linux)
        CC = g++
        OUT_FILE = tsp
endif

all:
		$(CC) $(SOURCES) -g -O2 -Wall -o $(OUT_FILE) -std=gnu++14

nearest_neighbor:
		$(CC) nearest_main.cpp k-opt.hpp k-opt.cpp  nearest_neighbor.cpp nearest_neighbor.hpp read.cpp read.hpp compute_distance.cpp compute_distance.hpp -g -O2 -Wall -o nn.out -std=gnu++14

clean:
		rm $(OUT_FILE)
