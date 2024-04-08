CXX = g++ ${CXXFLAGS}
CXXFLAGS = -std=c++11 -lpthread -I/usr/local/include/eigen3 -Wall
EXEC = sim.out

FILES = main.cpp

all:
	${CXX} ${FILES} -o ${EXEC}

.PHONY: clean
clean:
	rm -f *.o *~ ${EXEC}