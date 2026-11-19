CXX = g++
CXXFLAGS = -O2 -Wall

PTHREAD_FLAGS = -pthread
OPENMP_FLAGS = -fopenmp


PTHREAD_TARGET = src/pthreads_out
OPENMP_TARGET = src/openmp_out

all: $(PTHREAD_TARGET) $(OPENMP_TARGET)

$(PTHREAD_TARGET): src/pthreads.cpp
	$(CXX) $(CXXFLAGS) $(PTHREAD_FLAGS) src/pthreads.cpp -o $(PTHREAD_TARGET)

$(OPENMP_TARGET): src/openmp.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) src/openmp.cpp -o $(OPENMP_TARGET)
clean:
	rm -f $(PTHREAD_TARGET) $(OPENMP_TARGET_TARGET)

.PHONY: all clean
