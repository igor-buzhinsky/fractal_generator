CXXFLAGS = -std=c++20 -O3 -Wall -Wextra -Wpedantic -Wshadow -Wconversion $(EXTRA_CXXFLAGS)
LDLIBS = -lpng -pthread

HEADERS = $(wildcard src/*.h) third_party/cimg/CImg.h
OBJECTS = src/fractal.o src/palette.o src/cli.o src/image.o

all: fractal_generator

# Records the compiler and flags, so that changing them (e.g. EXTRA_CXXFLAGS for a
# sanitizer build) rebuilds everything instead of mixing old and new objects
.build_flags: FORCE
	@echo '$(CXX) $(CXXFLAGS) $(LDLIBS)' | cmp -s - $@ || echo '$(CXX) $(CXXFLAGS) $(LDLIBS)' > $@

fractal_generator: $(OBJECTS) src/main.o .build_flags
	$(CXX) $(CXXFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS)

tests/unit_tests: $(OBJECTS) tests/unit_tests.o .build_flags
	$(CXX) $(CXXFLAGS) $(filter %.o,$^) -o $@ $(LDLIBS)

%.o: %.cpp $(HEADERS) .build_flags
	$(CXX) $(CXXFLAGS) -Isrc -isystem third_party/cimg -c $< -o $@

test: fractal_generator tests/unit_tests
	./tests/unit_tests
	./tests/cli_tests.sh

clean:
	rm -f fractal_generator tests/unit_tests src/*.o tests/*.o .build_flags

.PHONY: all test clean FORCE
