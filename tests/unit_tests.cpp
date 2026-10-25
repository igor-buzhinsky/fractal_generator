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
	// empty and unknown names mean identity
	CHECK(near(str_to_function("")(z), z));
	CHECK(near(str_to_function("no_such_function")(z), z));
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
	CHECK(intensity(0, 0, 100) == 0);
	// -1 is inside as well (period-2 cycle)
	CHECK(intensity(-1, 0, 100) == 0);
	// far outside points escape immediately
	CHECK(intensity(2, 2, 100) == 1);
	// points closer to the set take more iterations to escape
	CHECK(intensity(0.3, 0, 100) > intensity(1, 0, 100));

	continuous_intensity = true;
	CHECK(intensity(0, 0, 100) == 0);
	double smooth = intensity(0.3, 0, 100);
	CHECK(smooth > 0 && smooth < 100);
	continuous_intensity = false;

	julia();
	CHECK(julia_on);
	// for c = 0, the filled Julia set is the unit disk
	julia_c = dcomplex(0, 0);
	CHECK(intensity(0.5, 0.5, 100) == 0);
	CHECK(intensity(1.5, 0, 100) > 0);
}

void test_presets() {
	mandelbrot();
	CHECK(!julia_on && g_z2_name == "sqr" && g_z1_name == "");
	tricorn();
	CHECK(!julia_on && g_z1_name == "conj");
	burning_ship();
	CHECK(!julia_on && g_z1_name == "component_abs");
	julia();
	CHECK(julia_on && near(julia_c, dcomplex(-0.7, 0.27015)));
	for (auto preset : {mandelbrot, tricorn, burning_ship, julia}) {
		preset();
		CHECK(x_min < x_max && y_min < y_max);
	}
}

void test_colors() {
	CHECK(same_color(color("black"), color(0, 0, 0)));
	CHECK(same_color(color("white"), color(255, 255, 255)));
	CHECK(same_color(color("yellow"), color(255, 255, 0)));
	CHECK(same_color(color("cyan"), color(0, 255, 255)));
	CHECK(same_color(color("unknown"), color(0, 0, 0)));

	CHECK(same_color(interpolate(color(0, 0, 0), color(200, 100, 50), 0), color(0, 0, 0)));
	CHECK(same_color(interpolate(color(0, 0, 0), color(200, 100, 50), 1), color(200, 100, 50)));
	CHECK(same_color(interpolate(color(0, 0, 0), color(200, 100, 50), 0.5), color(100, 50, 25)));
}

void test_palette() {
	palette p;
	p.fill(100, "black:red:white", color("blue"), 1.0);
	CHECK(p.colors.size() == 100);
	// index 0 is the body color, then the gradient goes from the first to the last color
	CHECK(same_color(p.colors[0], color("blue")));
	CHECK(p.colors[1].r < 10 && p.colors[1].g == 0 && p.colors[1].b == 0);
	CHECK(same_color(p.colors[99], color("white")));

	// the size must not depend on the number of colors in the gradient
	string sequence = "black";
	for (int n = 2; n <= 64; n++) {
		sequence += ":white";
		p.fill(50, sequence, color("black"), 1.0);
		CHECK(p.colors.size() == 50);
	}

	// skew does not change the endpoints
	p.fill(100, "black:red:white", color("blue"), 3.0);
	CHECK(p.colors.size() == 100);
	CHECK(same_color(p.colors[99], color("white")));
}

void test_zoom() {
	double x1 = -1, x2 = 3;
	zoom(x1, x2, 2);
	CHECK(near(x1, 0) && near(x2, 2));
	zoom(x1, x2, 0.5);
	CHECK(near(x1, -1) && near(x2, 3));
}

void test_recalc_iter() {
	set_square(-2, 2, -1, 1);
	points = 20000;
	recalc_iter();
	CHECK(x_iter == 200 && y_iter == 100);
}

void test_apply_name() {
	CHECK(apply_name("", "z") == "z");
	CHECK(apply_name("sqr", "z") == "sqr(z)");
	CHECK(apply_name("sin", apply_name("conj", "z")) == "sin(conj(z))");
}

void test_recalc_fractal() {
	mandelbrot();
	set_square(-2, 1, -1.5, 1.5);
	points = 10000;
	palette_size = 100;

	threads_num = 1;
	recalc_fractal();
	CHECK(x_iter == 100 && y_iter == 100);
	// the center column passes through the main cardioid, the image corner is outside
	CHECK(buffer[x_iter * 2 / 3][y_iter / 2] == 0);
	CHECK(buffer[0][0] > 0);

	vector<float> single_threaded;
	for (size_t x_i = 0; x_i < x_iter; x_i++) {
		single_threaded.insert(single_threaded.end(), buffer[x_i], buffer[x_i] + y_iter);
	}

	// multithreading must not change the result
	threads_num = 7;
	recalc_fractal();
	bool same = true;
	for (size_t x_i = 0; x_i < x_iter; x_i++) {
		for (size_t y_i = 0; y_i < y_iter; y_i++) {
			same &= buffer[x_i][y_i] == single_threaded[x_i * y_iter + y_i];
		}
	}
	CHECK(same);
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
	test_apply_name();
	test_recalc_fractal();

	cout << "Unit tests: " << g_checks - g_failures << "/" << g_checks << " checks passed" << endl;
	return g_failures == 0 ? 0 : 1;
}
