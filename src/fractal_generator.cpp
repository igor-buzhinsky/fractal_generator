//============================================================================
// Name        : fractal_generator.cpp
// Author      : Igor Buzhinsky
// Version     :
// Copyright   : 
//============================================================================

#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <pthread.h>

#include "CImg.h"

using namespace cimg_library;
using namespace std;

void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

double calc_coordinate(double min, double max, int i, int iter) {
	return min + (max - min) * ((double) i / (iter - 1));
}

typedef complex<double> dcomplex;
typedef dcomplex un_function(dcomplex);
string g_z1_name, g_z2_name, g_c_name, g_name;
un_function *g_z1, *g_z2, *g_c, *g;

dcomplex square(dcomplex z) {
	return z * z;
}
dcomplex invert(dcomplex z) {
	return dcomplex(1, 0) / z;
}
dcomplex core_function(dcomplex z, dcomplex c) {
	return g(g_z2(g_z1(z)) + g_c(c));
}

un_function *str_to_function(const string &str) {
	un_function *res;
	if (str == "sqr") {
		res = square;
	} else if (str == "component_abs") {
		res = [](dcomplex z) -> dcomplex { return dcomplex(abs(real(z)), abs(imag(z))); };
	} else if (str == "conj") {
		res = [](dcomplex z) -> dcomplex { return conj(z); };
	} else if (str == "cube") {
		res = [](dcomplex z) -> dcomplex { return z * z * z; };
	} else if (str == "inverse") {
		res = invert;
	} else if (str == "sin") {
		res = [](dcomplex z) -> dcomplex { return sin(z); };
	} else if (str == "cos") {
		res = [](dcomplex z) -> dcomplex { return cos(z); };
	} else if (str == "tan") {
		res = [](dcomplex z) -> dcomplex { return tan(z); };
	} else if (str == "cot") {
		res = [](dcomplex z) -> dcomplex { return cos(z) / sin(z); };
	} else if (str == "sinh") {
		res = [](dcomplex z) -> dcomplex { return sinh(z); };
	} else if (str == "cosh") {
		res = [](dcomplex z) -> dcomplex { return cosh(z); };
	} else if (str == "tanh") {
		res = [](dcomplex z) -> dcomplex { return tanh(z); };
	} else if (str == "coth") {
		res = [](dcomplex z) -> dcomplex { return cosh(z) / sinh(z); };
	} else if (str == "exp") {
		res = [](dcomplex z) -> dcomplex { return exp(z); };
	} else if (str == "log") {
		res = [](dcomplex z) -> dcomplex { return log(z); };
	} else {
		res = [](dcomplex z) -> dcomplex { return z; };
	}
	return res;
}

const double log_2 = log(2);
bool continuous_intensity = false;
bool julia_on = true;
dcomplex julia_c(0, 0);
double x_min, x_max, y_min, y_max, palette_skew = 1.0;
int points = 300000, palette_size = 1000, threads_num = 4;

double intensity(double x, double y, size_t palette_size) {
	dcomplex c = julia_on ? julia_c : dcomplex(x, y);
	dcomplex z0 = julia_on ? dcomplex(x, y) : dcomplex(0, 0);

	dcomplex v = core_function(z0, c);
	size_t i;

	if (continuous_intensity) {
		for (i = 1; i < (palette_size - 2) && abs(v) < (1 << 16); i++) {
			v = core_function(v, c);
		}

		if (i < palette_size - 2) {
			double log_zn = log(abs(v)) / 2;
			double nu = log(log_zn / log_2) / log_2;
			if (nu > 1) {
				nu = 1;
			} else if (nu < 0) {
				nu = 0;
			}
			return i + 1 - nu;
		}
	} else {
		for (i = 1; i < palette_size - 1; i++) {
			v = core_function(v, c);
			if (abs(v) > 2) {
				return i;
			}
		}
	}
	return 0;
}

struct color {
	unsigned char r, g, b;
	
	color(): r(0), g(0), b(0) {};
	
	color(unsigned char r, unsigned char g, unsigned char b): r(r), g(g), b(b) {};
	
	color(string str) {
		r = g = b = 0;
		if (str == "black") {
		} else if (str == "white") {
			r = g = b = 255;
		} else if (str == "red") {
			r = 255;
		} else if (str == "green") {
			g = 255;
		} else if (str == "blue") {
			b = 255;
		} else if (str == "yellow") {
			r = g = 255;
		} else if (str == "magenta") {
			b = r = 255;
		} else if (str == "cyan") {
			b = g = 255;
		} else {
		}
	}
};

color interpolate(color color1, color color2, double fraction) {
	return color {
		(unsigned char)(color1.r * (1 - fraction) + color2.r * fraction),
		(unsigned char)(color1.g * (1 - fraction) + color2.g * fraction),
		(unsigned char)(color1.b * (1 - fraction) + color2.b * fraction)
	};
}

void set_functions(const string &str_g, const string &str_g_z1, const string &str_g_z2,
		const string &str_g_c) {
	g_name = str_g;
	g_z1_name = str_g_z1;
	g_z2_name = str_g_z2;
	g_c_name = str_g_c;
	g = str_to_function(str_g);
	g_z1 = str_to_function(str_g_z1);
	g_z2 = str_to_function(str_g_z2);
	g_c = str_to_function(str_g_c);
}

void set_square(double x_min_, double x_max_, double y_min_, double y_max_) {
	x_min = x_min_;
	x_max = x_max_;
	y_min = y_min_;
	y_max = y_max_;
}

void mandelbrot() {
	set_functions("", "", "sqr", "");
	set_square(-1.6, 1.0, -1.05, 1.05);
	julia_on = false;
}

void tricorn() {
	set_functions("", "conj", "sqr", "");
	set_square(-1.6, 1.6, -1.6, 1.6);
	julia_on = false;
}

void burning_ship() {
	set_functions("", "component_abs", "sqr", "");
	set_square(-2.0, 1.0, -1.5, 0.5);
	julia_on = false;
}

void julia() {
	set_functions("", "", "sqr", "");
	set_square(-1.6, 1.6, -1.05, 1.05);
	julia_on = true;
	julia_c = dcomplex(-0.7, 0.27015);
}

struct palette {
	vector<color> colors;

	palette() {}
	
	void fill(size_t size, string color_sequence, color c_body, double skew) {
		colors.clear();
		colors.push_back(c_body);
		vector<string> str_colors;
		split(color_sequence, ':', str_colors);
		vector<color> reference_colors;
		for (size_t i = 0; i < str_colors.size(); i++) {
			reference_colors.push_back(color(str_colors[i]));
		}
		double interval_length = 1. / (reference_colors.size() - 1);
		for (size_t i = 1; i < size; i++) {
			double fraction = pow((double) i / (size - 1), skew);
			for (size_t j = 1; j < reference_colors.size(); j++) {
				double border = j * interval_length;
				if (fraction <= border) {
					color c1 = reference_colors[j - 1];
					color c2 = reference_colors[j];
					color c = interpolate(c1, c2, 1 - (border - fraction) / interval_length);
					colors.push_back(c);
					break;
				}
			}
		}
	}

	void apply(CImg<unsigned char> &img, size_t x_i, size_t y_i, double index) const {
		int int_index = (int) floor(index);
		color c1 = colors[int_index];
		color c2 = colors[int_index + 1];
		color c = interpolate(c1, c2, index - int_index);
		img(x_i, y_i, 0) = c.r;
		img(x_i, y_i, 1) = c.g;
		img(x_i, y_i, 2) = c.b;
	}
};

void zoom(double &x1, double &x2, double value) {
	double radius = (x2 - x1) / 2;
	double center = (x1 + x2) / 2;
	x1 = center - radius / value;
	x2 = center + radius / value;
}

void print_help() {
	cout << "Commands:" << endl;
	cout << "    (x_min | x_max | y_min | y_max) <double>" << endl;
	cout << "    (left | right | top | bottom | zoom) <double>" << endl;
	cout << "    palette_size <integer>" << endl;
	cout << "    palette_skew <double>" << endl;
	cout << "    points <integer>" << endl;
	cout << "    colors <color_1:color_2:...color_n>" << endl;
	cout << "    body_color <color>" << endl;
	cout << "    threads_num <integer>" << endl;
	cout << "    continuous_intensity <0 or 1>" << endl;
	cout << "    julia_on <0 or 1>" << endl;
	cout << "    julia_c <re:double> <im:double>" << endl;
	cout << "    (mandelbrot | julia | tricorn | burning_ship)" << endl;
	cout << "    (plot | help | exit)" << endl;
	cout << "    (g | g_z1 | g_z2 | g_c) <function_name>" << endl;
}

string color_sequence("black:red:yellow:green:blue:cyan:magenta:white");
string color_body_name = "black";
color color_body(color_body_name);
size_t x_iter, y_iter;
palette pal;
float **buffer;

void print_state() {
	printf("x in [%.17f, %.17f];\n", x_min, x_max);
	printf("y in [%.17f, %.17f];\n", y_min, y_max);
	printf("points = %d; palette_size = %d; palette_skew = %.2f; threads_num = %d;\n",
			points, palette_size, palette_skew, threads_num);
	printf("body_color = %s; continuous_intensity = %i; \n",
			color_body_name.c_str(), (int) continuous_intensity);
	printf("colors = %s;\n",
			color_sequence.c_str());
	printf("f_c(z) = %s(%s(%s(z)) + %s(c)); \n",
				g_name.c_str(), g_z2_name.c_str(), g_z1_name.c_str(), g_c_name.c_str());
	printf("julia_on = %i; julia_c = (%.6f, %.6f)\n",
			(int) julia_on, real(julia_c), imag(julia_c));
}

void *run(void *ptr) {
	int thread_index = *(int *)ptr;
	for (size_t x_i = thread_index; x_i < x_iter; x_i += threads_num) {
		buffer[x_i] = new float[y_iter];
		double x = calc_coordinate(x_min, x_max, x_i, x_iter);
		for (size_t y_i = 0; y_i < y_iter; y_i++) {
			double y = calc_coordinate(y_min, y_max, y_i, y_iter);
			float color = (float) intensity(x, y, palette_size);
			buffer[x_i][y_i] = color;
		}
	}
	return NULL;
}

void read_args(vector<string> &arguments, int number) {
	arguments.resize(number);
	for (int i = 0; i < number; i++) {
		cin >> arguments[i];
	}
}

void clean_buffer() {
	if (buffer != NULL) {
		for (size_t x_i = 0; x_i < x_iter; x_i++) {
			delete[] buffer[x_i];
		}
		delete[] buffer;
	}
}

void recalc_iter() {
	double x_range = x_max - x_min;
	double y_range = y_max - y_min;

	// x_iter * y_iter = calculations
	// x_iter / y_iter = x_range / y_range
	// x_iter = x_range / y_range * y_iter
	// x_range / y_range * y_iter^2 = calculations
	// y_iter = sqrt(calculations * y_range / x_range)
	y_iter = (size_t) sqrt(points * y_range / x_range);
	x_iter = (size_t) sqrt(points * x_range / y_range);
}

void recalc_fractal() {
	clean_buffer();
	recalc_iter();
	buffer = new float*[x_iter];
	vector<pthread_t> threads;
	int *thread_indices = new int[threads_num];
	for (int i = 0; i < threads_num; i++) {
		threads.push_back(0);
		thread_indices[i] = i;
		int code = pthread_create(&threads[i], NULL, run, (void *) &(thread_indices[i]));
		if (code) {
			cerr << "pthreads error!" << endl;
			exit(1);
		}
	}

	for (int i = 0; i < threads_num; i++) {
		pthread_join(threads[i], NULL);
	}
	delete[] thread_indices;
}

void prepare_image() {
	pal.fill(palette_size, color_sequence, color_body, palette_skew);
	CImg<unsigned char> img(x_iter, y_iter, 1, 3, 0);
	for (size_t x_i = 0; x_i < x_iter; x_i++) {
		for (size_t y_i = 0; y_i < y_iter; y_i++) {
			pal.apply(img, x_i, y_i, buffer[x_i][y_i]);
		}
	}
	img.save("out.png");
}

int main(int argc, char *argv[]) {
	print_help();
	julia();
	
	string command;
	vector<string> arguments;
	bool first_time = true;
	while (true) {
		arguments.clear();
		if (first_time) {
			command = "plot";
			first_time = false;
		} else {
			cout << "> ";
			cin >> command;
			int arg_number;
			if (command == "exit" || command == "help" || command == "plot" || command == "mandelbrot"
					 || command == "julia" || command == "tricorn" || command == "burning_ship") {
				arg_number = 0;
			} else if (command == "julia_c") {
				arg_number = 2;
			} else {
				arg_number = 1;
			}
			read_args(arguments, arg_number);
		}
		
		bool must_recalc_fractal = true;

		if (command == "exit") {
			clean_buffer();
			return 0;
		} else if (command == "help") {
			print_help();
			continue;
		} else if (command == "plot") {
			// do nothing
		} else if (command == "mandelbrot") {
			mandelbrot();
		} else if (command == "julia") {
			julia();
		} else if (command == "tricorn") {
			tricorn();
		} else if (command == "burning_ship") {
			burning_ship();
		} else if (command == "x_min") {
			double value = atof(arguments[0].c_str());
			if (value >= x_max) {
				cerr << "x_min must be less than x_max!" << endl;
				continue;
			}
			x_min = value;
		} else if (command == "x_max") {
			double value = atof(arguments[0].c_str());
			if (value <= x_min) {
				cerr << "x_max must be greater than x_min!" << endl;
				continue;
			}
			x_max = value;
		} else if (command == "y_min") {
			double value = atof(arguments[0].c_str());
			if (value >= y_max) {
				cerr << "y_min must be less than y_max!" << endl;
				continue;
			}
			y_min = value;
		} else if (command == "y_max") {
			double value = atof(arguments[0].c_str());
			if (value <= y_min) {
				cerr << "y_max must be greater than y_min!" << endl;
				continue;
			}
			y_max = value;
		} else if (command == "julia_c") {
			double re = atof(arguments[0].c_str());
			double im = atof(arguments[1].c_str());
			julia_c = dcomplex(re, im);
		} else if (command == "palette_size") {
			int value = atoi(arguments[0].c_str());
			if (value < 10) {
				cerr << "The value must be at least 10!" << endl;
				continue;
			}
			palette_size = value;
		} else if (command == "points") {
			int value = atoi(arguments[0].c_str());
			if (value < 100) {
				cerr << "The value must be at least 100!" << endl;
				continue;
			}
			points = value;
		} else if (command == "g") {
			set_functions(arguments[0], g_z1_name, g_z2_name, g_c_name);
		} else if (command == "g_z1") {
			set_functions(g_name, arguments[0], g_z2_name, g_c_name);
		} else if (command == "g_z2") {
			set_functions(g_name, g_z1_name, arguments[0], g_c_name);
		} else if (command == "g_c") {
			set_functions(g_name, g_z1_name, g_z2_name, arguments[0]);
		} else if (command == "colors") {
			color_sequence = arguments[0];
			// TODO check correctness
			must_recalc_fractal = false;
		} else if (command == "body_color") {
			color_body_name = arguments[0];
			color_body = color(color_body_name);
			must_recalc_fractal = false;
			// TODO check correctness
		} else if (command == "zoom") {
			double value = atof(arguments[0].c_str());
			if (value <= 0) {
				cerr << "The value must be positive!" << endl;
				continue;
			}
			zoom(x_min, x_max, value);
			zoom(y_min, y_max, value);
		} else if (command == "left") {
			double value = atof(arguments[0].c_str());
			if (value <= 0) {
				cerr << "The value must be positive!" << endl;
				continue;
			}
			x_max -= (x_max - x_min) * (1 - value);
		} else if (command == "right") {
			double value = atof(arguments[0].c_str());
			if (value <= 0) {
				cerr << "The value must be positive!" << endl;
				continue;
			}
			x_min += (x_max - x_min) * (1 - value);
		} else if (command == "top") {
			double value = atof(arguments[0].c_str());
			if (value <= 0) {
				cerr << "The value must be positive!" << endl;
				continue;
			}
			y_max -= (y_max - y_min) * (1 - value);
		} else if (command == "bottom") {
			double value = atof(arguments[0].c_str());
			if (value <= 0) {
				cerr << "The value must be positive!" << endl;
				continue;
			}
			y_min += (y_max - y_min) * (1 - value);
		} else if (command == "palette_skew") {
			double value = atof(arguments[0].c_str());
			if (value <= 0) {
				cerr << "The value must be positive!" << endl;
				continue;
			}
			palette_skew = value;
			must_recalc_fractal = false;
		} else if (command == "continuous_intensity") {
			int value = atoi(arguments[0].c_str());
			if (value != 0 && value != 1) {
				cerr << "Only valid values are 0 and 1!" << endl;
				continue;
			}
			continuous_intensity = (bool) value;
		} else if (command == "julia_on") {
			int value = atoi(arguments[0].c_str());
			if (value != 0 && value != 1) {
				cerr << "Only valid values are 0 and 1!" << endl;
				continue;
			}
			julia_on = (bool) value;
		} else if (command == "threads_num") {
			int value = atoi(arguments[0].c_str());
			if (value < 1 || value > 1024) {
				cerr << "The value must be between 1 and 1024!" << endl;
				continue;
			}
			threads_num = value;
			continue;
		} else {
			cerr << "Unknown command! Type 'help' to see the list of available ones." << endl;
			continue;
		}
		
		if (must_recalc_fractal) {
			recalc_fractal();
		}
		
		prepare_image();
		print_state();
	}
}
