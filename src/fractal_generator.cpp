//============================================================================
// Name        : fractal_generator.cpp
// Author      : Igor Buzhinsky
//============================================================================

#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <pthread.h>

#include "CImg.h"

using namespace cimg_library;
using namespace std;

void split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}

// Formats a number with the given number of significant digits, without trailing zeros
string format_number(double value, int precision = 15) {
	char result[64];
	snprintf(result, sizeof(result), "%.*g", precision, value);
	return result;
}

double calc_coordinate(double min, double max, int i, int count) {
	return min + (max - min) * ((double) i / (count - 1));
}

typedef complex<double> dcomplex;
typedef dcomplex complex_function(dcomplex);
string g_z1_name, g_z2_name, g_c_name, g_name;
complex_function *g_z1, *g_z2, *g_c, *g;

dcomplex square(dcomplex z) {
	return z * z;
}
dcomplex invert(dcomplex z) {
	return dcomplex(1, 0) / z;
}
dcomplex core_function(dcomplex z, dcomplex c) {
	return g(g_z2(g_z1(z)) + g_c(c));
}

complex_function *str_to_function(const string &str) {
	complex_function *res;
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

// A point escapes once |z| exceeds this radius
const double ESCAPE_RADIUS = 2;
// For smooth coloring, iterate further: the smoothing formula is only accurate for large |z|
const double SMOOTH_ESCAPE_RADIUS = 1 << 16;

bool continuous_intensity = false;
bool julia_on = true;
dcomplex julia_c(0, 0);
double x_min, x_max, y_min, y_max, palette_skew = 1.0;
// approximate number of pixels in the image (the 'points' command)
int pixel_count = 300000;
int palette_size = 1000, threads_num = 4;

// Returns the palette index for the point (x, y): 0 if the point does not escape
// (it is in the set), otherwise the number of iterations before it escaped, which is
// fractional for smooth coloring. The result never exceeds palette_size - 2, so that
// palette::apply() can always interpolate between colors at index and index + 1.
double iteration_count(double x, double y, size_t max_palette_size) {
	dcomplex c = julia_on ? julia_c : dcomplex(x, y);
	dcomplex z0 = julia_on ? dcomplex(x, y) : dcomplex(0, 0);

	dcomplex v = core_function(z0, c);
	size_t i;

	if (continuous_intensity) {
		for (i = 1; i < (max_palette_size - 2) && abs(v) < SMOOTH_ESCAPE_RADIUS; i++) {
			v = core_function(v, c);
		}

		if (i < max_palette_size - 2) {
			// normalized iteration count: subtract how far |z| overshot the bailout
			double nu = log2(log2(abs(v))) - 1;
			nu = min(max(nu, 0.0), 1.0);
			return i + 1 - nu;
		}
	} else {
		for (i = 1; i < max_palette_size - 1; i++) {
			v = core_function(v, c);
			if (abs(v) > ESCAPE_RADIUS) {
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

	color(const string &name) {
		r = g = b = 0;
		if (name == "white") {
			r = g = b = 255;
		} else if (name == "red") {
			r = 255;
		} else if (name == "green") {
			g = 255;
		} else if (name == "blue") {
			b = 255;
		} else if (name == "yellow") {
			r = g = 255;
		} else if (name == "magenta") {
			b = r = 255;
		} else if (name == "cyan") {
			b = g = 255;
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

void set_view(double new_x_min, double new_x_max, double new_y_min, double new_y_max) {
	x_min = new_x_min;
	x_max = new_x_max;
	y_min = new_y_min;
	y_max = new_y_max;
}

void mandelbrot() {
	set_functions("", "", "sqr", "");
	set_view(-1.6, 1.0, -1.05, 1.05);
	julia_on = false;
}

void tricorn() {
	set_functions("", "conj", "sqr", "");
	set_view(-1.6, 1.6, -1.6, 1.6);
	julia_on = false;
}

void burning_ship() {
	set_functions("", "component_abs", "sqr", "");
	set_view(-2.0, 1.0, -1.5, 0.5);
	julia_on = false;
}

void julia() {
	set_functions("", "", "sqr", "");
	set_view(-1.6, 1.6, -1.05, 1.05);
	julia_on = true;
	julia_c = dcomplex(-0.7, 0.27015);
}

struct palette {
	// colors[0] is the body color, the rest is the gradient for escaping points
	vector<color> colors;

	void fill(size_t size, const string &color_sequence, color c_body, double skew) {
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
				// the last border may be slightly below 1 due to rounding
				if (fraction <= border || j == reference_colors.size() - 1) {
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
	cout << R"(Presets:
    mandelbrot | julia | tricorn | burning_ship
                                     Load a fractal with its default view
View:
    x_min | x_max | y_min | y_max <double>
                                     Set a boundary of the visible area
    zoom <double>                    Zoom in around the center (< 1 zooms out)
    left | right | top | bottom <double>
                                     Keep only this fraction of the view on that side,
                                     e.g. 'left 0.5' keeps the left half
Image:
    points <integer>                 Total number of pixels in the output image
    palette_size <integer>           Number of colors, also the maximum iteration count
    palette_skew <double>            Exponent of color distribution (1 = uniform)
    colors <color_1:...:color_n>     Gradient for points outside the set; available colors:
                                     black, white, red, green, blue, yellow, magenta, cyan
    body_color <color>               Color of points inside the set
    continuous_intensity <0 or 1>    Smooth coloring instead of discrete bands
Formula f_c(z) = g(g_z2(g_z1(z)) + g_c(c)):
    g | g_z1 | g_z2 | g_c <function_name>
                                     Replace one of the functions in the formula;
                                     available: sqr, cube, inverse, conj, component_abs,
                                     sin, cos, tan, cot, sinh, cosh, tanh, coth, exp, log;
                                     any other name means identity
    julia_on <0 or 1>                1: Julia set (c is fixed, z_0 is the point),
                                     0: Mandelbrot-like set (z_0 = 0, c is the point)
    julia_c <re:double> <im:double>  Constant c used for Julia sets
Other:
    threads_num <integer>            Number of threads used for computation
    plot                             Redraw the image
    help                             Show this list
    exit                             Quit the program
Each new image is saved to out.png.
)";
}

string color_sequence("black:red:yellow:green:blue:cyan:magenta:white");
string color_body_name = "black";
color color_body(color_body_name);
size_t image_width, image_height;
palette pal;
float **buffer;

// Function application for display: an empty name stands for identity
string apply_name(const string &function_name, const string &argument) {
	return function_name.empty() ? argument : function_name + "(" + argument + ")";
}

void print_state() {
	string formula = apply_name(g_name,
			apply_name(g_z2_name, apply_name(g_z1_name, "z")) + " + " + apply_name(g_c_name, "c"));
	cout << endl;
	cout << "Formula:  z -> " << formula << endl;
	if (julia_on) {
		cout << "          Julia set: z starts at the point, c = " << format_number(real(julia_c))
				<< (imag(julia_c) < 0 ? " - " : " + ") << format_number(fabs(imag(julia_c)))
				<< "i (julia_c)" << endl;
	} else {
		cout << "          Mandelbrot-like set: z starts at 0, c is the point" << endl;
	}
	cout << "View:     x in [" << format_number(x_min) << ", " << format_number(x_max)
			<< "], y in [" << format_number(y_min) << ", " << format_number(y_max) << "]" << endl;
	cout << "          center (" << format_number((x_min + x_max) / 2) << ", "
			<< format_number((y_min + y_max) / 2) << "), width " << format_number(x_max - x_min, 6)
			<< ", height " << format_number(y_max - y_min, 6) << endl;
	cout << "Image:    " << image_width << " x " << image_height << " pixels (points = " << pixel_count
			<< "), saved to out.png" << endl;
	cout << "Colors:   " << color_sequence << " (colors), inside the set: " << color_body_name
			<< " (body_color)" << endl;
	cout << "          palette_size = " << palette_size << ", palette_skew = " << format_number(palette_skew)
			<< ", continuous_intensity = " << (continuous_intensity ? "on" : "off") << endl;
	cout << "Threads:  " << threads_num << endl;
}

void *run(void *ptr) {
	int thread_index = *(int *)ptr;
	for (size_t x_i = thread_index; x_i < image_width; x_i += threads_num) {
		buffer[x_i] = new float[image_height];
		double x = calc_coordinate(x_min, x_max, x_i, image_width);
		for (size_t y_i = 0; y_i < image_height; y_i++) {
			double y = calc_coordinate(y_min, y_max, y_i, image_height);
			buffer[x_i][y_i] = (float) iteration_count(x, y, palette_size);
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
		for (size_t x_i = 0; x_i < image_width; x_i++) {
			delete[] buffer[x_i];
		}
		delete[] buffer;
	}
}

// Chooses the image size so that it has about pixel_count pixels and the aspect ratio of the view
void compute_image_size() {
	double aspect_ratio = (x_max - x_min) / (y_max - y_min);
	image_height = (size_t) sqrt(pixel_count / aspect_ratio);
	image_width = (size_t) sqrt(pixel_count * aspect_ratio);
}

void recalc_fractal() {
	clean_buffer();
	compute_image_size();
	buffer = new float*[image_width];
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
	CImg<unsigned char> img(image_width, image_height, 1, 3, 0);
	for (size_t x_i = 0; x_i < image_width; x_i++) {
		for (size_t y_i = 0; y_i < image_height; y_i++) {
			pal.apply(img, x_i, y_i, buffer[x_i][y_i]);
		}
	}
	img.save("out.png");
}

// Unit tests include this file and provide their own main()
#ifndef FRACTAL_GENERATOR_NO_MAIN
int main() {
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
			if (!(cin >> command)) {
				// end of input, e.g. Ctrl+D or the end of a piped script
				command = "exit";
			}
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
			pixel_count = value;
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
			must_recalc_fractal = false;
		} else if (command == "body_color") {
			color_body_name = arguments[0];
			color_body = color(color_body_name);
			must_recalc_fractal = false;
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
#endif
