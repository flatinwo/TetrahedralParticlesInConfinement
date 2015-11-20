CXX=g++
CXXFLAGS=-O3 -W -Wall -Werror -g -std=c++0x
LDFLAGS=-g
SRC= $(wildcard *.cpp)
OBJ= $(SRC:.cpp=.o)

OUT=../libpatchycolloids.a

$(OUT): $(OBJ)
	    ar rcs $(OUT) $(OBJ) 

%.o: %.c
	    @$(CXX) -o $@ -c $< $(CXXFLAGS)

.PHONY: clean

clean:
	    @rm -rf *.o
