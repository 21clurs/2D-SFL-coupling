CXX = g++ ${CXXFLAGS}
CXXFLAGS = -std=c++11 -lpthread -I/usr/local/include/eigen3 -Wall -O3
EXEC = sim.out

FILES = main.cpp liquidmesh.cpp rigidbody.cpp testinghelpers.cpp sim.cpp boundaryintegral.cpp mesh.cpp scenes.cpp simoptions.cpp

all:
	${CXX} ${FILES} -o ${EXEC}

.PHONY: clean
clean:
	rm -f *.o *~ ${EXEC}