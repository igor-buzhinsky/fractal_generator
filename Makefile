CXXFLAGS = -std=c++0x -O3 -Wall
LDLIBS = -lX11 -lpthread

all: fractal_generator

fractal_generator: src/fractal_generator.cpp src/CImg.h
	g++ $(CXXFLAGS) src/fractal_generator.cpp -o $@ $(LDLIBS)

tests/unit_tests: tests/unit_tests.cpp src/fractal_generator.cpp src/CImg.h
	g++ $(CXXFLAGS) tests/unit_tests.cpp -o $@ $(LDLIBS)

test: fractal_generator tests/unit_tests
	./tests/unit_tests
	./tests/cli_tests.sh

clean:
	rm -f fractal_generator tests/unit_tests

.PHONY: all test clean
