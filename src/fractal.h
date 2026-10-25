// Fractal definition and computation: formula, view, presets and per-pixel iteration counts.

#pragma once

#include <complex>
#include <string>
#include <vector>

typedef std::complex<double> dcomplex;
typedef dcomplex complex_function(dcomplex);

// A function that can be used in the formula
struct named_function {
	std::string name; // empty for the identity function
	complex_function *apply;
};

// Names of all available functions, "id" being the identity function
std::vector<std::string> function_names();

// Looks up a function by name; returns false if there is no such function
bool find_function(const std::string &name, named_function &result);

// f_c(z) = g(g_z2(g_z1(z)) + g_c(c))
struct formula {
	named_function g, g_z1, g_z2, g_c;

	// Each name must be valid for find_function(), "" is the identity function
	formula(const std::string &g, const std::string &g_z1, const std::string &g_z2,
			const std::string &g_c);

	dcomplex apply(dcomplex z, dcomplex c) const {
		return g.apply(g_z2.apply(g_z1.apply(z)) + g_c.apply(c));
	}

	// Human-readable form, e.g. "sqr(z) + c"
	std::string to_string() const;
};

// Visible area of the complex plane
struct view {
	double x_min, x_max, y_min, y_max;

	double width() const { return x_max - x_min; }
	double height() const { return y_max - y_min; }

	// Zooms in around the center by the given factor (< 1 zooms out)
	void zoom(double factor);
};

struct fractal_settings {
	formula f = formula("", "", "sqr", "");
	// true: Julia set, z starts at the point and c = julia_c;
	// false: Mandelbrot-like set, z starts at 0 and c is the point
	bool julia_on = false;
	dcomplex julia_c = dcomplex(0, 0);
	view area = {-2, 2, -2, 2};
	// approximate number of pixels in the image (the 'points' command)
	int pixel_count = 300000;
	// maximum iteration count, equal to the number of colors in the palette
	int palette_size = 1000;
	// fractional iteration counts for smooth coloring
	bool continuous_intensity = false;
};

// Names of all presets
std::vector<std::string> preset_names();

// Sets the formula, view and mode of a preset; returns false if there is no such preset
bool apply_preset(const std::string &name, fractal_settings &settings);

// Returns the palette index for a point: 0 if the point does not escape (it is in the set),
// otherwise the number of iterations before it escaped, which is fractional for smooth
// coloring. The result never exceeds palette_size - 2, so that palette::at() can always
// interpolate between colors at index and index + 1.
double iteration_count(const fractal_settings &settings, dcomplex point);

// Each image side is at most this long, to keep extremely wide or tall images reasonable
const size_t MAX_IMAGE_SIDE = 16384;

struct image_size {
	size_t width, height;
};

// Chooses the image size so that it has about pixel_count pixels and the aspect ratio
// of the view. Each side is at least 2 pixels, so that coordinates can be computed for
// both edges, and at most MAX_IMAGE_SIDE.
image_size compute_image_size(const view &area, int pixel_count);

// Iteration counts of all pixels
struct iteration_grid {
	size_t width = 0, height = 0;
	std::vector<float> values; // row by row, top row first

	float at(size_t x_i, size_t y_i) const { return values[y_i * width + x_i]; }
};

// Computes iteration counts for the whole image using the given number of threads
iteration_grid compute_fractal(const fractal_settings &settings, int threads_num);
