CC = g++-7
OUT_FILE = a.out

DIRS := . Blossom5 Blossom5/MinCost
SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))

uname_S=$(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(uname_S),Linux)
        CC = g++
        OUT_FILE = tsp
endif

all:
		$(CC) $(SOURCES) -g -O2 -Wall -o $(OUT_FILE) -std=gnu++14

clean:
		rm $(OUT_FILE)