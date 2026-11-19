CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wshadow -Wconversion $(EXTRA_CXXFLAGS)
LDLIBS = -lpng -pthread

HEADERS = $(wildcard src/*.h) third_party/cimg/CImg.h
OBJECTS = src/fractal.o src/palette.o src/cli.o

all: fractal_generator

fractal_generator: $(OBJECTS) src/main.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

tests/unit_tests: $(OBJECTS) tests/unit_tests.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -Isrc -isystem third_party/cimg -c $< -o $@

test: fractal_generator tests/unit_tests
	./tests/unit_tests
	./tests/cli_tests.sh

clean:
	rm -f fractal_generator tests/unit_tests src/*.o tests/*.o

.PHONY: all test clean
