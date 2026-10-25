#include "fractal.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>

using namespace std;

namespace {

const map<string, complex_function *> FUNCTIONS = {
	{"", [](dcomplex z) { return z; }},
	{"sqr", [](dcomplex z) { return z * z; }},
	{"cube", [](dcomplex z) { return z * z * z; }},
	{"inverse", [](dcomplex z) { return dcomplex(1, 0) / z; }},
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

named_function function_by_name(const string &name) {
	named_function result;
	find_function(name, result);
	return result;
}

// Function application for display: an empty name stands for identity
string apply_name(const string &function_name, const string &argument) {
	return function_name.empty() ? argument : function_name + "(" + argument + ")";
}

struct preset {
	formula f;
	view area;
	bool julia_on;
};

const map<string, preset> PRESETS = {
	{"mandelbrot", {formula("", "", "sqr", ""), {-1.6, 1.0, -1.05, 1.05}, false}},
	{"julia", {formula("", "", "sqr", ""), {-1.6, 1.6, -1.05, 1.05}, true}},
	{"tricorn", {formula("", "conj", "sqr", ""), {-1.6, 1.6, -1.6, 1.6}, false}},
	// Burning ship is traditionally drawn with the imaginary axis pointing down,
	// conj(c) mirrors it to the usual orientation
	{"burning_ship", {formula("", "component_abs", "sqr", "conj"), {-2.0, 1.0, -0.5, 1.5}, false}},
};

const dcomplex DEFAULT_JULIA_C(-0.7, 0.27015);

// A point escapes once |z| exceeds this radius
const double ESCAPE_RADIUS = 2;
// For smooth coloring, iterate further: the smoothing formula is only accurate for large |z|
const double SMOOTH_ESCAPE_RADIUS = 1 << 16;

// Coordinate of pixel i out of count, so that pixel 0 is at first and the last one at last
double pixel_coordinate(double first, double last, size_t i, size_t count) {
	return first + (last - first) * ((double) i / (count - 1));
}

}

vector<string> function_names() {
	vector<string> names;
	for (const auto &entry : FUNCTIONS) {
		names.push_back(entry.first.empty() ? "id" : entry.first);
	}
	return names;
}

bool find_function(const string &name, named_function &result) {
	string stored_name = name == "id" ? "" : name;
	auto entry = FUNCTIONS.find(stored_name);
	if (entry == FUNCTIONS.end()) {
		return false;
	}
	result = {stored_name, entry->second};
	return true;
}

formula::formula(const string &g, const string &g_z1, const string &g_z2, const string &g_c)
	: g(function_by_name(g)), g_z1(function_by_name(g_z1)), g_z2(function_by_name(g_z2)),
	g_c(function_by_name(g_c)) {}

string formula::to_string() const {
	return apply_name(g.name, apply_name(g_z2.name, apply_name(g_z1.name, "z")) + " + "
			+ apply_name(g_c.name, "c"));
}

void view::zoom(double factor) {
	double x_center = (x_min + x_max) / 2, x_radius = width() / 2 / factor;
	double y_center = (y_min + y_max) / 2, y_radius = height() / 2 / factor;
	*this = {x_center - x_radius, x_center + x_radius, y_center - y_radius, y_center + y_radius};
}

vector<string> preset_names() {
	vector<string> names;
	for (const auto &entry : PRESETS) {
		names.push_back(entry.first);
	}
	return names;
}

bool apply_preset(const string &name, fractal_settings &settings) {
	auto entry = PRESETS.find(name);
	if (entry == PRESETS.end()) {
		return false;
	}
	settings.f = entry->second.f;
	settings.area = entry->second.area;
	settings.julia_on = entry->second.julia_on;
	if (settings.julia_on) {
		settings.julia_c = DEFAULT_JULIA_C;
	}
	return true;
}

double iteration_count(const fractal_settings &settings, dcomplex point) {
	const formula &f = settings.f;
	size_t palette_size = settings.palette_size;
	dcomplex c = settings.julia_on ? settings.julia_c : point;
	dcomplex v = f.apply(settings.julia_on ? point : dcomplex(0, 0), c);
	size_t i;

	if (settings.continuous_intensity) {
		for (i = 1; i < palette_size - 2 && abs(v) < SMOOTH_ESCAPE_RADIUS; i++) {
			v = f.apply(v, c);
		}

		if (i < palette_size - 2) {
			// normalized iteration count: subtract how far |z| overshot the bailout
			double nu = log2(log2(abs(v))) - 1;
			nu = min(max(nu, 0.0), 1.0);
			return i + 1 - nu;
		}
	} else {
		for (i = 1; i < palette_size - 1; i++) {
			v = f.apply(v, c);
			if (abs(v) > ESCAPE_RADIUS) {
				return i;
			}
		}
	}
	return 0;
}

image_size compute_image_size(const view &area, int pixel_count) {
	const double min_side = 2;
	const double max_side = min(pixel_count / min_side, (double) MAX_IMAGE_SIDE);
	double aspect_ratio = area.width() / area.height();
	return {
		(size_t) min(max(sqrt(pixel_count * aspect_ratio), min_side), max_side),
		(size_t) min(max(sqrt(pixel_count / aspect_ratio), min_side), max_side)
	};
}

iteration_grid compute_fractal(const fractal_settings &settings, int threads_num) {
	image_size size = compute_image_size(settings.area, settings.pixel_count);
	iteration_grid grid;
	grid.width = size.width;
	grid.height = size.height;
	grid.values.resize(grid.width * grid.height);
	const view &area = settings.area;

	// thread i computes rows i, i + threads_num, ...
	auto compute_rows = [&](int thread_index) {
		for (size_t y_i = thread_index; y_i < grid.height; y_i += threads_num) {
			// image rows go from top to bottom, so the y axis points up as usual
			double y = pixel_coordinate(area.y_max, area.y_min, y_i, grid.height);
			for (size_t x_i = 0; x_i < grid.width; x_i++) {
				double x = pixel_coordinate(area.x_min, area.x_max, x_i, grid.width);
				grid.values[y_i * grid.width + x_i] = (float) iteration_count(settings, dcomplex(x, y));
			}
		}
	};

	vector<thread> threads;
	for (int i = 0; i < threads_num; i++) {
		threads.emplace_back(compute_rows, i);
	}
	for (thread &t : threads) {
		t.join();
	}
	return grid;
}
