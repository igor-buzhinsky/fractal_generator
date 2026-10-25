all:
	g++ -std=c++0x -O3 -Wall src/*.cpp -o fractal_generator -lX11 -lpthread
