CXX = g++
CXXFLAGS = -O2 -Wall -pthread

TARGET = src/out

all: $(TARGET)

$(TARGET): src/pthreads.cpp
	$(CXX) $(CXXFLAGS) src/pthreads.cpp -o $(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all clean
