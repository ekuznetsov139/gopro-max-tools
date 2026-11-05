CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra

all: test_maxcam

test_maxcam: ../tests/test_maxcam.cpp maxcam.cpp maxcam.h
	$(CXX) $(CXXFLAGS) ../tests/test_maxcam.cpp maxcam.cpp -o $@

test: test_maxcam
	./test_maxcam ../calibration/model4.txt

clean:
	rm -f test_maxcam

.PHONY: all test clean
