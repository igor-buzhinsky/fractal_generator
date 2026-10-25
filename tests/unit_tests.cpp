//============================================================================
// Unit tests for the computation and command code. Run with: make test
//============================================================================

#include <cmath>
#include <iostream>
#include <sstream>

#include "cli.h"
#include "fractal.h"
#include "palette.h"

using namespace std;

int g_checks = 0, g_failures = 0;

#define CHECK(condition) do { \
	g_checks++; \
	if (!(condition)) { \
		g_failures++; \
		cerr << __FILE__ << ":" << __LINE__ << ": check failed: " #condition << endl; \
	} \
} while (0)

bool near(double a, double b, double eps = 1e-9) {
	return fabs(a - b) < eps;
}

bool near(dcomplex a, dcomplex b, double eps = 1e-9) {
	return abs(a - b) < eps;
}

bool same_color(color a, color b) {
	return a.r == b.r && a.g == b.g && a.b == b.b;
}

color named_color(const string &name) {
	color result;
	CHECK(find_color(name, result));
	return result;
}

fractal_settings preset_settings(const string &name) {
	fractal_settings settings;
	CHECK(apply_preset(name, settings));
	return settings;
}

void test_functions() {
	dcomplex z(1, 2);
	named_function f;
	CHECK(find_function("sqr", f) && f.name == "sqr" && near(f.apply(z), z * z));
	CHECK(find_function("cube", f) && near(f.apply(z), z * z * z));
	CHECK(find_function("conj", f) && near(f.apply(z), dcomplex(1, -2)));
	CHECK(find_function("component_abs", f) && near(f.apply(dcomplex(-1, -2)), dcomplex(1, 2)));
	CHECK(find_function("inverse", f) && near(f.apply(z), dcomplex(1, 0) / z));
	CHECK(find_function("exp", f) && near(f.apply(z), exp(z)));
	CHECK(find_function("cot", f) && near(f.apply(z), cos(z) / sin(z)));
	// identity is stored as an empty name, so that it is not shown in the formula
	CHECK(find_function("id", f) && f.name == "" && near(f.apply(z), z));
	CHECK(!find_function("no_such_function", f));

	// every listed name can be found
	for (const string &name : function_names()) {
		CHECK(find_function(name, f));
	}
}

void test_formula() {
	formula mandelbrot("", "", "sqr", "");
	CHECK(near(mandelbrot.apply(dcomplex(1, 2), dcomplex(3, 4)), dcomplex(1, 2) * dcomplex(1, 2) + dcomplex(3, 4)));
	CHECK(mandelbrot.to_string() == "sqr(z) + c");

	formula tricorn("", "conj", "sqr", "");
	CHECK(near(tricorn.apply(dcomplex(1, 2), dcomplex(0, 0)), dcomplex(1, -2) * dcomplex(1, -2)));
	CHECK(tricorn.to_string() == "sqr(conj(z)) + c");

	// the outer function g is applied last
	formula outer("exp", "", "", "sin");
	CHECK(near(outer.apply(dcomplex(1, 0), dcomplex(1, 0)), exp(dcomplex(1, 0) + sin(dcomplex(1, 0)))));
	CHECK(outer.to_string() == "exp(z + sin(c))");
}

void test_iteration_count() {
	fractal_settings settings = preset_settings("mandelbrot");
	settings.palette_size = 100;
	// 0 is inside the Mandelbrot set
	CHECK(iteration_count(settings, dcomplex(0, 0)) == 0);
	// -1 is inside as well (period-2 cycle)
	CHECK(iteration_count(settings, dcomplex(-1, 0)) == 0);
	// far outside points escape immediately
	CHECK(iteration_count(settings, dcomplex(2, 2)) == 1);
	// points closer to the set take more iterations to escape
	CHECK(iteration_count(settings, dcomplex(0.3, 0)) > iteration_count(settings, dcomplex(1, 0)));
	// the result never exceeds palette_size - 2
	for (double x = -2; x <= 1; x += 0.01) {
		CHECK(iteration_count(settings, dcomplex(x, 0.01)) <= settings.palette_size - 2);
	}

	settings.continuous_intensity = true;
	CHECK(iteration_count(settings, dcomplex(0, 0)) == 0);
	double smooth = iteration_count(settings, dcomplex(0.3, 0));
	CHECK(smooth > 0 && smooth <= settings.palette_size - 2);

	// for c = 0, the filled Julia set is the unit disk
	settings = preset_settings("julia");
	settings.julia_c = dcomplex(0, 0);
	CHECK(iteration_count(settings, dcomplex(0.5, 0.5)) == 0);
	CHECK(iteration_count(settings, dcomplex(1.5, 0)) > 0);
}

void test_presets() {
	fractal_settings settings = preset_settings("mandelbrot");
	CHECK(!settings.julia_on && settings.f.to_string() == "sqr(z) + c");
	settings = preset_settings("tricorn");
	CHECK(!settings.julia_on && settings.f.to_string() == "sqr(conj(z)) + c");
	settings = preset_settings("burning_ship");
	CHECK(!settings.julia_on && settings.f.to_string() == "sqr(component_abs(z)) + conj(c)");
	settings = preset_settings("julia");
	CHECK(settings.julia_on && near(settings.julia_c, dcomplex(-0.7, 0.27015)));
	CHECK(!apply_preset("no_such_preset", settings));

	CHECK(preset_names().size() == 4);
	for (const string &name : preset_names()) {
		settings = preset_settings(name);
		CHECK(settings.area.x_min < settings.area.x_max && settings.area.y_min < settings.area.y_max);
	}
}

void test_colors() {
	CHECK(same_color(named_color("black"), color{0, 0, 0}));
	CHECK(same_color(named_color("white"), color{255, 255, 255}));
	CHECK(same_color(named_color("yellow"), color{255, 255, 0}));
	CHECK(same_color(named_color("cyan"), color{0, 255, 255}));
	color unused;
	CHECK(!find_color("unknown", unused));
	CHECK(color_names().size() == 8);

	CHECK(same_color(interpolate(color{0, 0, 0}, color{200, 100, 50}, 0), color{0, 0, 0}));
	CHECK(same_color(interpolate(color{0, 0, 0}, color{200, 100, 50}, 1), color{200, 100, 50}));
	CHECK(same_color(interpolate(color{0, 0, 0}, color{200, 100, 50}, 0.5), color{100, 50, 25}));

	CHECK(is_color_sequence("black:white"));
	CHECK(is_color_sequence("black:red:yellow:green:blue:cyan:magenta:white"));
	CHECK(!is_color_sequence("black"));
	CHECK(!is_color_sequence("black:purple"));
	CHECK(!is_color_sequence("black::white"));
	CHECK(!is_color_sequence(""));
}

void test_palette() {
	palette p(100, "black:red:white", named_color("blue"), 1.0);
	CHECK(p.colors.size() == 100);
	// index 0 is the body color, then the gradient goes from the first to the last color
	CHECK(same_color(p.colors[0], named_color("blue")));
	CHECK(p.colors[1].r < 10 && p.colors[1].g == 0 && p.colors[1].b == 0);
	CHECK(same_color(p.colors[99], named_color("white")));

	// fractional indices interpolate between neighboring colors
	CHECK(same_color(p.at(0), named_color("blue")));
	CHECK(same_color(p.at(97.5), interpolate(p.colors[97], p.colors[98], 0.5)));

	// the size must not depend on the number of colors in the gradient
	string sequence = "black";
	for (int n = 2; n <= 64; n++) {
		sequence += ":white";
		CHECK(palette(50, sequence, named_color("black"), 1.0).colors.size() == 50);
	}

	// skew does not change the endpoints
	palette skewed(100, "black:red:white", named_color("blue"), 3.0);
	CHECK(skewed.colors.size() == 100);
	CHECK(same_color(skewed.colors[99], named_color("white")));
}

void test_view() {
	view area = {-1, 3, 0, 2};
	CHECK(near(area.width(), 4) && near(area.height(), 2));
	area.zoom(2);
	CHECK(near(area.x_min, 0) && near(area.x_max, 2) && near(area.y_min, 0.5) && near(area.y_max, 1.5));
	area.zoom(0.5);
	CHECK(near(area.x_min, -1) && near(area.x_max, 3) && near(area.y_min, 0) && near(area.y_max, 2));
}

void test_image_size() {
	image_size size = compute_image_size({-2, 2, -1, 1}, 20000);
	CHECK(size.width == 200 && size.height == 100);
	// extremely wide view: the height must still be at least 2 pixels
	size = compute_image_size({0, 1, 0, 1e-9}, 10000);
	CHECK(size.width == 5000 && size.height == 2);
	// extremely tall view
	size = compute_image_size({0, 1e-9, 0, 1}, 10000);
	CHECK(size.width == 2 && size.height == 5000);
	// a side is never longer than image converters accept
	size = compute_image_size({0, 1, 0, 1e-9}, 300000);
	CHECK(size.width == MAX_IMAGE_SIDE && size.height == 2);
}

void test_compute_fractal() {
	fractal_settings settings = preset_settings("mandelbrot");
	settings.area = {-2, 1, -1.5, 1.5};
	settings.pixel_count = 10000;
	settings.palette_size = 100;

	iteration_grid single_threaded = compute_fractal(settings, 1);
	CHECK(single_threaded.width == 100 && single_threaded.height == 100);
	CHECK(single_threaded.values.size() == 100 * 100);
	// the center row passes through the main cardioid, the image corner is outside
	CHECK(single_threaded.at(single_threaded.width * 2 / 3, single_threaded.height / 2) == 0);
	CHECK(single_threaded.at(0, 0) > 0);

	// multithreading must not change the result, even with more threads than rows
	CHECK(compute_fractal(settings, 7).values == single_threaded.values);
	settings.pixel_count = 100;
	CHECK(compute_fractal(settings, 64).values == compute_fractal(settings, 1).values);
}

void test_y_axis_orientation() {
	// upper half of the Mandelbrot set: it touches the bottom edge of the image, not the top
	fractal_settings settings = preset_settings("mandelbrot");
	settings.area = {-2, 1, 0, 1.5};
	settings.pixel_count = 20000;
	settings.palette_size = 100;
	iteration_grid grid = compute_fractal(settings, 2);
	size_t x_zero = grid.width * 2 / 3;
	CHECK(grid.at(x_zero, grid.height - 1) == 0);
	CHECK(grid.at(x_zero, 0) > 0);
}

void test_parse_numbers() {
	double d;
	CHECK(parse_double("1.5", d) && d == 1.5);
	CHECK(parse_double("-2e-3", d) && d == -2e-3);
	CHECK(!parse_double("", d));
	CHECK(!parse_double("abc", d));
	CHECK(!parse_double("1.5x", d));
	CHECK(!parse_double("nan", d));
	CHECK(!parse_double("inf", d));

	int i;
	CHECK(parse_int("42", i) && i == 42);
	CHECK(parse_int("-7", i) && i == -7);
	CHECK(!parse_int("", i));
	CHECK(!parse_int("4.2", i));
	CHECK(!parse_int("abc", i));
	CHECK(!parse_int("99999999999999999999", i));

	CHECK(format_number(1.6) == "1.6");
	CHECK(format_number(-0.00105) == "-0.00105");
	CHECK(format_number(2.123456789, 3) == "2.12");
}

// Runs a command line that must succeed and returns its outcome
outcome run(app_state &state, const string &line) {
	stringstream errors;
	outcome result = run_command_line(state, line, errors);
	CHECK(errors.str().empty());
	return result;
}

// Checks that the command is rejected with the given error and does not change the state
void check_rejected(const string &line, const string &expected_error) {
	app_state state;
	apply_preset("mandelbrot", state.fractal);
	stringstream before, after, errors;
	print_state(before, state, {0, 0});
	CHECK(run_command_line(state, line, errors) == outcome::nothing);
	print_state(after, state, {0, 0});
	CHECK(before.str() == after.str());
	g_checks++;
	if (errors.str().find(expected_error) == string::npos) {
		g_failures++;
		cerr << "'" << line << "': expected error '" << expected_error << "', got '" << errors.str() << "'" << endl;
	}
}

void test_commands() {
	app_state state;
	CHECK(run(state, "mandelbrot") == outcome::recompute);
	CHECK(!state.fractal.julia_on);
	CHECK(run(state, "") == outcome::nothing);

	CHECK(run(state, "x_min -2") == outcome::recompute);
	CHECK(run(state, "x_max 2") == outcome::recompute);
	CHECK(run(state, "  y_min   -1  ") == outcome::recompute);
	CHECK(run(state, "y_max 1") == outcome::recompute);
	CHECK(near(state.fractal.area.x_min, -2) && near(state.fractal.area.x_max, 2));
	CHECK(near(state.fractal.area.y_min, -1) && near(state.fractal.area.y_max, 1));

	CHECK(run(state, "zoom 2") == outcome::recompute);
	CHECK(near(state.fractal.area.x_min, -1) && near(state.fractal.area.y_max, 0.5));
	CHECK(run(state, "left 0.5") == outcome::recompute);
	CHECK(near(state.fractal.area.x_min, -1) && near(state.fractal.area.x_max, 0));
	CHECK(run(state, "top 0.5") == outcome::recompute);
	CHECK(near(state.fractal.area.y_min, 0) && near(state.fractal.area.y_max, 0.5));
	CHECK(run(state, "right 0.5") == outcome::recompute);
	CHECK(near(state.fractal.area.x_min, -0.5) && near(state.fractal.area.x_max, 0));
	CHECK(run(state, "bottom 0.5") == outcome::recompute);
	CHECK(near(state.fractal.area.y_min, 0) && near(state.fractal.area.y_max, 0.25));

	CHECK(run(state, "g_z1 conj") == outcome::recompute);
	CHECK(run(state, "g sin") == outcome::recompute);
	CHECK(state.fractal.f.to_string() == "sin(sqr(conj(z)) + c)");
	CHECK(run(state, "g id") == outcome::recompute);
	CHECK(state.fractal.f.to_string() == "sqr(conj(z)) + c");

	CHECK(run(state, "julia_on 1") == outcome::recompute);
	CHECK(run(state, "julia_c 0.3 -0.5") == outcome::recompute);
	CHECK(state.fractal.julia_on && near(state.fractal.julia_c, dcomplex(0.3, -0.5)));

	// color changes do not require computing the fractal again
	CHECK(run(state, "colors black:white") == outcome::recolor);
	CHECK(run(state, "body_color red") == outcome::recolor);
	CHECK(run(state, "palette_skew 2") == outcome::recolor);
	CHECK(state.colors.colors == "black:white" && state.colors.body_color == "red" && state.colors.skew == 2);

	CHECK(run(state, "points 1000") == outcome::recompute && state.fractal.pixel_count == 1000);
	CHECK(run(state, "palette_size 50") == outcome::recompute && state.fractal.palette_size == 50);
	CHECK(run(state, "continuous_intensity 1") == outcome::recompute && state.fractal.continuous_intensity);
	CHECK(run(state, "threads_num 3") == outcome::nothing && state.threads_num == 3);
	CHECK(run(state, "plot") == outcome::recompute);
	CHECK(run(state, "exit") == outcome::quit);

	check_rejected("no_such_command", "no_such_command: unknown command");
	check_rejected("zoom", "zoom: expected 1 argument(s): <double>");
	check_rejected("zoom 2 3", "zoom: expected 1 argument(s): <double>");
	check_rejected("julia_c 1", "julia_c: expected 2 argument(s): <re:double> <im:double>");
	check_rejected("zoom abc", "zoom: expected a number, got 'abc'");
	check_rejected("zoom 0", "zoom: the value must be positive");
	check_rejected("left -1", "left: the value must be positive");
	check_rejected("x_min 100", "x_min: must be less than x_max");
	check_rejected("y_max -100", "y_max: must be greater than y_min");
	check_rejected("points 10", "points: expected an integer from 100 to 100000000, got '10'");
	check_rejected("points 1.5", "points: expected an integer from 100 to 100000000, got '1.5'");
	check_rejected("palette_size 1", "palette_size: expected an integer from 10 to 10000000, got '1'");
	check_rejected("continuous_intensity 2", "continuous_intensity: expected an integer from 0 to 1, got '2'");
	check_rejected("threads_num 0", "threads_num: expected an integer from 1 to 1024, got '0'");
	check_rejected("julia_c 1 x", "julia_c: expected a number, got 'x'");
	check_rejected("colors black", "colors: expected at least two known colors separated by ':', got 'black'");
	check_rejected("body_color purple", "body_color: unknown color 'purple'");
	check_rejected("g_z1 no_such_function", "g_z1: unknown function 'no_such_function'");
}

void test_help() {
	stringstream help;
	print_help(help);
	// every command is listed
	for (const command &c : commands()) {
		for (const string &name : c.names) {
			CHECK(help.str().find(name) != string::npos);
			CHECK(find_command(name) == &c);
		}
	}
	CHECK(find_command("no_such_command") == nullptr);
	// no line is too long
	string line;
	while (getline(help, line)) {
		CHECK(line.size() <= 93);
	}
}

int main() {
	test_functions();
	test_formula();
	test_iteration_count();
	test_presets();
	test_colors();
	test_palette();
	test_view();
	test_image_size();
	test_compute_fractal();
	test_y_axis_orientation();
	test_parse_numbers();
	test_commands();
	test_help();

	cout << "Unit tests: " << g_checks - g_failures << "/" << g_checks << " checks passed" << endl;
	return g_failures == 0 ? 0 : 1;
}
