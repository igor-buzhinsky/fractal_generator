all:
	g++ -std=c++0x -O3 -Wall -lX11 -lpthread src/*.cpp -o fractal_generator
