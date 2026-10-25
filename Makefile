CXXFLAGS = -std=c++17 -O3 -Wall
LDLIBS = -lX11 -pthread

HEADERS = $(wildcard src/*.h)
OBJECTS = src/fractal.o src/palette.o src/cli.o

all: fractal_generator

fractal_generator: $(OBJECTS) src/main.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

# Unit tests do not use CImg, so they do not need X11
tests/unit_tests: $(OBJECTS) tests/unit_tests.o
	$(CXX) $(CXXFLAGS) $^ -o $@ -pthread

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -Isrc -c $< -o $@

test: fractal_generator tests/unit_tests
	./tests/unit_tests
	./tests/cli_tests.sh

clean:
	rm -f fractal_generator tests/unit_tests src/*.o tests/*.o

.PHONY: all test clean
