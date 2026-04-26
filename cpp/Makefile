CXX    = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17

.PHONY: all demo test clean

all: demo

demo: demo.cpp universalfft.hpp
	$(CXX) $(CXXFLAGS) -o ufft_demo_cpp demo.cpp

test: demo
	mkdir -p ../tests/data
	./ufft_demo_cpp ../tests/data

clean:
	rm -f ufft_demo_cpp *.csv
