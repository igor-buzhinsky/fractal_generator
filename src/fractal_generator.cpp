//============================================================================
// Name        : fractal_generator.cpp
// Author      : Igor Buzhinsky
//============================================================================

#include <algorithm>
#include <climits>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <map>
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

// Parses the whole text as a finite number
bool parse_double(const string &text, double &value) {
	char *end;
	value = strtod(text.c_str(), &end);
	return !text.empty() && *end == '\0' && isfinite(value);
}

// Parses the whole text as an integer that fits into int
bool parse_int(const string &text, int &value) {
	char *end;
	long result = strtol(text.c_str(), &end, 10);
	value = (int) result;
	return !text.empty() && *end == '\0' && result >= INT_MIN && result <= INT_MAX;
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

// Functions available in the formula; "id" (identity) is stored as an empty name
const map<string, complex_function *> FUNCTIONS = {
	{"", [](dcomplex z) { return z; }},
	{"sqr", square},
	{"cube", [](dcomplex z) { return z * z * z; }},
	{"inverse", invert},
	{"conj", [](dcomplex z) { return conj(z); }},
	{"component_abs", [](dcomplex z) { return dcomplex(abs(real(z)), abs(imag(z))); }},
	{"sin", [](dcomplex z) { return sin(z); }},
	{"cos", [](dcomplex z) { return cos(z); }},
	{"tan", [](dcomplex z) { return tan(z); }},
	{"cot", [](dcomplex z) { return cos(z) / sin(z); }},
	{"sinh", [](dcomplex z) { return sinh(z); }},
	{"cosh", [](dcomplex z) { return cosh(z); }},
	{"tanh", [](dcomplex z) { return tanh(z); }},
	{"coth", [](dcomplex z) { return cosh(z) / sinh(z); }},
	{"exp", [](dcomplex z) { return exp(z); }},
	{"log", [](dcomplex z) { return log(z); }},
};

// Returns the stored name of a formula function, or false if there is no such function
bool normalize_function_name(const string &name, string &result) {
	result = name == "id" ? "" : name;
	return FUNCTIONS.count(result) > 0;
}

complex_function *str_to_function(const string &name) {
	return FUNCTIONS.at(name);
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
};

const map<string, color> COLORS = {
	{"black", color(0, 0, 0)},
	{"white", color(255, 255, 255)},
	{"red", color(255, 0, 0)},
	{"green", color(0, 255, 0)},
	{"blue", color(0, 0, 255)},
	{"yellow", color(255, 255, 0)},
	{"magenta", color(255, 0, 255)},
	{"cyan", color(0, 255, 255)},
};

bool is_color_name(const string &name) {
	return COLORS.count(name) > 0;
}

// Checks a sequence like "black:red:white": at least two known colors
bool is_color_sequence(const string &sequence) {
	vector<string> names;
	split(sequence, ':', names);
	return names.size() >= 2 && all_of(names.begin(), names.end(), is_color_name);
}

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

// Burning ship is traditionally drawn with the imaginary axis pointing down,
// conj(c) mirrors it to the usual orientation
void burning_ship() {
	set_functions("", "component_abs", "sqr", "conj");
	set_view(-2.0, 1.0, -0.5, 1.5);
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
			reference_colors.push_back(COLORS.at(str_colors[i]));
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
    colors <color_1:...:color_n>     Gradient (2+ colors) for points outside the set; colors:
                                     black, white, red, green, blue, yellow, magenta, cyan
    body_color <color>               Color of points inside the set
    continuous_intensity <0 or 1>    Smooth coloring instead of discrete bands
Formula f_c(z) = g(g_z2(g_z1(z)) + g_c(c)):
    g | g_z1 | g_z2 | g_c <function_name>
                                     Replace one of the functions in the formula;
                                     available: id (identity), sqr, cube, inverse, conj,
                                     component_abs, sin, cos, tan, cot, sinh, cosh, tanh,
                                     coth, exp, log
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
color color_body = COLORS.at(color_body_name);
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
			// image rows go from top to bottom, so the y axis points up as usual
			double y = calc_coordinate(y_max, y_min, y_i, image_height);
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

void print_error(const string &command, const string &message) {
	cerr << command << ": " << message << endl;
}

// The read_* functions parse a command argument and print an error if it is invalid
bool read_double(const string &command, const string &text, double &value) {
	if (!parse_double(text, value)) {
		print_error(command, "expected a number, got '" + text + "'");
		return false;
	}
	return true;
}

bool read_positive_double(const string &command, const string &text, double &value) {
	if (!read_double(command, text, value)) {
		return false;
	}
	if (value <= 0) {
		print_error(command, "the value must be positive");
		return false;
	}
	return true;
}

bool read_int(const string &command, const string &text, int min_value, int max_value, int &value) {
	if (!parse_int(text, value) || value < min_value || value > max_value) {
		print_error(command, "expected an integer from " + to_string(min_value) + " to "
				+ to_string(max_value) + ", got '" + text + "'");
		return false;
	}
	return true;
}

void clean_buffer() {
	if (buffer != NULL) {
		for (size_t x_i = 0; x_i < image_width; x_i++) {
			delete[] buffer[x_i];
		}
		delete[] buffer;
	}
}

// Chooses the image size so that it has about pixel_count pixels and the aspect ratio of the view.
// Each side is at least 2 pixels, so that coordinates can be computed for both edges,
// and at most MAX_IMAGE_SIDE, which image converters accept.
const double MAX_IMAGE_SIDE = 16384;

void compute_image_size() {
	const double min_side = 2, max_side = min(pixel_count / min_side, MAX_IMAGE_SIDE);
	double aspect_ratio = (x_max - x_min) / (y_max - y_min);
	image_height = (size_t) min(max(sqrt(pixel_count / aspect_ratio), min_side), max_side);
	image_width = (size_t) min(max(sqrt(pixel_count * aspect_ratio), min_side), max_side);
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
	try {
		img.save("out.png");
	} catch (const CImgException &e) {
		cerr << "Could not save out.png: " << e.what() << endl;
	}
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
		} else if (command == "x_min" || command == "x_max" || command == "y_min" || command == "y_max") {
			double value;
			if (!read_double(command, arguments[0], value)) {
				continue;
			}
			bool is_x = command[0] == 'x', is_min = command.substr(2) == "min";
			double &bound = is_x ? (is_min ? x_min : x_max) : (is_min ? y_min : y_max);
			double other = is_x ? (is_min ? x_max : x_min) : (is_min ? y_max : y_min);
			if (is_min ? value >= other : value <= other) {
				print_error(command, string("must be ") + (is_min ? "less than " : "greater than ")
						+ command.substr(0, 2) + (is_min ? "max" : "min"));
				continue;
			}
			bound = value;
		} else if (command == "julia_c") {
			double re, im;
			if (!read_double(command, arguments[0], re) || !read_double(command, arguments[1], im)) {
				continue;
			}
			julia_c = dcomplex(re, im);
		} else if (command == "palette_size") {
			int value;
			if (!read_int(command, arguments[0], 10, 10000000, value)) {
				continue;
			}
			palette_size = value;
		} else if (command == "points") {
			int value;
			if (!read_int(command, arguments[0], 100, 100000000, value)) {
				continue;
			}
			pixel_count = value;
		} else if (command == "g" || command == "g_z1" || command == "g_z2" || command == "g_c") {
			string name;
			if (!normalize_function_name(arguments[0], name)) {
				print_error(command, "unknown function '" + arguments[0] + "'");
				continue;
			}
			set_functions(command == "g" ? name : g_name, command == "g_z1" ? name : g_z1_name,
					command == "g_z2" ? name : g_z2_name, command == "g_c" ? name : g_c_name);
		} else if (command == "colors") {
			if (!is_color_sequence(arguments[0])) {
				print_error(command, "expected at least two known colors separated by ':', got '"
						+ arguments[0] + "'");
				continue;
			}
			color_sequence = arguments[0];
			must_recalc_fractal = false;
		} else if (command == "body_color") {
			if (!is_color_name(arguments[0])) {
				print_error(command, "unknown color '" + arguments[0] + "'");
				continue;
			}
			color_body_name = arguments[0];
			color_body = COLORS.at(color_body_name);
			must_recalc_fractal = false;
		} else if (command == "zoom" || command == "left" || command == "right" || command == "top"
				|| command == "bottom") {
			double value;
			if (!read_positive_double(command, arguments[0], value)) {
				continue;
			}
			if (command == "zoom") {
				zoom(x_min, x_max, value);
				zoom(y_min, y_max, value);
			} else if (command == "left") {
				x_max -= (x_max - x_min) * (1 - value);
			} else if (command == "right") {
				x_min += (x_max - x_min) * (1 - value);
			} else if (command == "top") {
				y_min += (y_max - y_min) * (1 - value);
			} else {
				y_max -= (y_max - y_min) * (1 - value);
			}
		} else if (command == "palette_skew") {
			double value;
			if (!read_positive_double(command, arguments[0], value)) {
				continue;
			}
			palette_skew = value;
			must_recalc_fractal = false;
		} else if (command == "continuous_intensity" || command == "julia_on") {
			int value;
			if (!read_int(command, arguments[0], 0, 1, value)) {
				continue;
			}
			(command == "julia_on" ? julia_on : continuous_intensity) = value;
		} else if (command == "threads_num") {
			int value;
			if (!read_int(command, arguments[0], 1, 1024, value)) {
				continue;
			}
			threads_num = value;
			continue;
		} else {
			print_error(command, "unknown command, type 'help' to see the list of available ones");
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
