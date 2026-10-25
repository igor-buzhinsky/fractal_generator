//============================================================================
// Unit tests for fractal_generator.cpp. Run with: make test
//============================================================================

#define FRACTAL_GENERATOR_NO_MAIN
#include "../src/fractal_generator.cpp"

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

void test_split() {
	vector<string> parts;
	split("black:red:white", ':', parts);
	CHECK(parts.size() == 3);
	CHECK(parts[0] == "black" && parts[1] == "red" && parts[2] == "white");

	parts.clear();
	split("single", ':', parts);
	CHECK(parts.size() == 1 && parts[0] == "single");
}

void test_calc_coordinate() {
	CHECK(near(calc_coordinate(-2, 2, 0, 5), -2));
	CHECK(near(calc_coordinate(-2, 2, 2, 5), 0));
	CHECK(near(calc_coordinate(-2, 2, 4, 5), 2));
}

void test_functions() {
	dcomplex z(1, 2);
	CHECK(near(str_to_function("sqr")(z), z * z));
	CHECK(near(str_to_function("cube")(z), z * z * z));
	CHECK(near(str_to_function("conj")(z), dcomplex(1, -2)));
	CHECK(near(str_to_function("component_abs")(dcomplex(-1, -2)), dcomplex(1, 2)));
	CHECK(near(str_to_function("inverse")(z), dcomplex(1, 0) / z));
	CHECK(near(str_to_function("exp")(z), exp(z)));
	CHECK(near(str_to_function("cot")(z), cos(z) / sin(z)));
	// identity is stored as an empty name
	CHECK(near(str_to_function("")(z), z));

	string name;
	CHECK(normalize_function_name("id", name) && name == "");
	CHECK(normalize_function_name("sin", name) && name == "sin");
	CHECK(!normalize_function_name("no_such_function", name));
}

void test_core_function() {
	set_functions("", "", "sqr", "");
	CHECK(near(core_function(dcomplex(1, 2), dcomplex(3, 4)), dcomplex(1, 2) * dcomplex(1, 2) + dcomplex(3, 4)));

	// tricorn: conj(z)^2 + c
	set_functions("", "conj", "sqr", "");
	CHECK(near(core_function(dcomplex(1, 2), dcomplex(0, 0)), dcomplex(1, -2) * dcomplex(1, -2)));

	// outer function g is applied last
	set_functions("exp", "", "", "");
	CHECK(near(core_function(dcomplex(1, 0), dcomplex(1, 0)), exp(dcomplex(2, 0))));
}

void test_intensity() {
	mandelbrot();
	continuous_intensity = false;
	// 0 is inside the Mandelbrot set
	CHECK(iteration_count(0, 0, 100) == 0);
	// -1 is inside as well (period-2 cycle)
	CHECK(iteration_count(-1, 0, 100) == 0);
	// far outside points escape immediately
	CHECK(iteration_count(2, 2, 100) == 1);
	// points closer to the set take more iterations to escape
	CHECK(iteration_count(0.3, 0, 100) > iteration_count(1, 0, 100));

	continuous_intensity = true;
	CHECK(iteration_count(0, 0, 100) == 0);
	double smooth = iteration_count(0.3, 0, 100);
	CHECK(smooth > 0 && smooth < 100);
	continuous_intensity = false;

	julia();
	CHECK(julia_on);
	// for c = 0, the filled Julia set is the unit disk
	julia_c = dcomplex(0, 0);
	CHECK(iteration_count(0.5, 0.5, 100) == 0);
	CHECK(iteration_count(1.5, 0, 100) > 0);
}

void test_presets() {
	mandelbrot();
	CHECK(!julia_on && g_z2_name == "sqr" && g_z1_name == "");
	tricorn();
	CHECK(!julia_on && g_z1_name == "conj");
	burning_ship();
	CHECK(!julia_on && g_z1_name == "component_abs" && g_c_name == "conj");
	julia();
	CHECK(julia_on && near(julia_c, dcomplex(-0.7, 0.27015)));
	for (auto preset : {mandelbrot, tricorn, burning_ship, julia}) {
		preset();
		CHECK(x_min < x_max && y_min < y_max);
	}
}

void test_colors() {
	CHECK(same_color(COLORS.at("black"), color(0, 0, 0)));
	CHECK(same_color(COLORS.at("white"), color(255, 255, 255)));
	CHECK(same_color(COLORS.at("yellow"), color(255, 255, 0)));
	CHECK(same_color(COLORS.at("cyan"), color(0, 255, 255)));

	CHECK(is_color_name("magenta"));
	CHECK(!is_color_name("unknown"));
	CHECK(is_color_sequence("black:white"));
	CHECK(is_color_sequence("black:red:yellow:green:blue:cyan:magenta:white"));
	CHECK(!is_color_sequence("black"));
	CHECK(!is_color_sequence("black:purple"));
	CHECK(!is_color_sequence("black::white"));
	CHECK(!is_color_sequence(""));

	CHECK(same_color(interpolate(color(0, 0, 0), color(200, 100, 50), 0), color(0, 0, 0)));
	CHECK(same_color(interpolate(color(0, 0, 0), color(200, 100, 50), 1), color(200, 100, 50)));
	CHECK(same_color(interpolate(color(0, 0, 0), color(200, 100, 50), 0.5), color(100, 50, 25)));
}

void test_palette() {
	palette p;
	p.fill(100, "black:red:white", COLORS.at("blue"), 1.0);
	CHECK(p.colors.size() == 100);
	// index 0 is the body color, then the gradient goes from the first to the last color
	CHECK(same_color(p.colors[0], COLORS.at("blue")));
	CHECK(p.colors[1].r < 10 && p.colors[1].g == 0 && p.colors[1].b == 0);
	CHECK(same_color(p.colors[99], COLORS.at("white")));

	// the size must not depend on the number of colors in the gradient
	string sequence = "black";
	for (int n = 2; n <= 64; n++) {
		sequence += ":white";
		p.fill(50, sequence, COLORS.at("black"), 1.0);
		CHECK(p.colors.size() == 50);
	}

	// skew does not change the endpoints
	p.fill(100, "black:red:white", COLORS.at("blue"), 3.0);
	CHECK(p.colors.size() == 100);
	CHECK(same_color(p.colors[99], COLORS.at("white")));
}

void test_zoom() {
	double x1 = -1, x2 = 3;
	zoom(x1, x2, 2);
	CHECK(near(x1, 0) && near(x2, 2));
	zoom(x1, x2, 0.5);
	CHECK(near(x1, -1) && near(x2, 3));
}

void test_recalc_iter() {
	set_view(-2, 2, -1, 1);
	pixel_count = 20000;
	compute_image_size();
	CHECK(image_width == 200 && image_height == 100);
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
}

void test_image_size_limits() {
	pixel_count = 10000;
	// extremely wide view: the height must still be at least 2 pixels
	set_view(0, 1, 0, 1e-9);
	compute_image_size();
	CHECK(image_height == 2 && image_width == 5000);
	// extremely tall view
	set_view(0, 1e-9, 0, 1);
	compute_image_size();
	CHECK(image_width == 2 && image_height == 5000);
	// a side is never longer than image converters accept
	pixel_count = 300000;
	set_view(0, 1, 0, 1e-9);
	compute_image_size();
	CHECK(image_width == MAX_IMAGE_SIDE && image_height == 2);
}

void test_apply_name() {
	CHECK(apply_name("", "z") == "z");
	CHECK(apply_name("sqr", "z") == "sqr(z)");
	CHECK(apply_name("sin", apply_name("conj", "z")) == "sin(conj(z))");
}

void test_recalc_fractal() {
	mandelbrot();
	set_view(-2, 1, -1.5, 1.5);
	pixel_count = 10000;
	palette_size = 100;

	threads_num = 1;
	recalc_fractal();
	CHECK(image_width == 100 && image_height == 100);
	// the center column passes through the main cardioid, the image corner is outside
	CHECK(buffer[image_width * 2 / 3][image_height / 2] == 0);
	CHECK(buffer[0][0] > 0);

	vector<float> single_threaded;
	for (size_t x_i = 0; x_i < image_width; x_i++) {
		single_threaded.insert(single_threaded.end(), buffer[x_i], buffer[x_i] + image_height);
	}

	// multithreading must not change the result
	threads_num = 7;
	recalc_fractal();
	bool same = true;
	for (size_t x_i = 0; x_i < image_width; x_i++) {
		for (size_t y_i = 0; y_i < image_height; y_i++) {
			same &= buffer[x_i][y_i] == single_threaded[x_i * image_height + y_i];
		}
	}
	CHECK(same);
	clean_buffer();
	buffer = NULL;
}

void test_y_axis_orientation() {
	// upper half of the Mandelbrot set: it touches the bottom edge of the image, not the top
	mandelbrot();
	set_view(-2, 1, 0, 1.5);
	pixel_count = 20000;
	palette_size = 100;
	threads_num = 2;
	recalc_fractal();
	size_t x_zero = image_width * 2 / 3;
	CHECK(buffer[x_zero][image_height - 1] == 0);
	CHECK(buffer[x_zero][0] > 0);
	clean_buffer();
	buffer = NULL;
}

int main() {
	test_split();
	test_calc_coordinate();
	test_functions();
	test_core_function();
	test_intensity();
	test_presets();
	test_colors();
	test_palette();
	test_zoom();
	test_recalc_iter();
	test_parse_numbers();
	test_image_size_limits();
	test_apply_name();
	test_recalc_fractal();
	test_y_axis_orientation();

	cout << "Unit tests: " << g_checks - g_failures << "/" << g_checks << " checks passed" << endl;
	return g_failures == 0 ? 0 : 1;
}
