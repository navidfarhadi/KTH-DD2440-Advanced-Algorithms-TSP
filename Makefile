CC = g++-7
OUT_FILE = a.out

uname_S=$(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(uname_S),Linux)
        CC = g++
        OUT_FILE = tsp
endif

all:
		$(CC) *.cpp -g -O2 -Wall -o $(OUT_FILE) -std=gnu++14

clean:
		rm $(OUT_FILE)