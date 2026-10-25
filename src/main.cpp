//============================================================================
// Name        : fractal_generator
// Author      : Igor Buzhinsky
//============================================================================

#include <iostream>
#include <string>

#include "cli.h"
#include "fractal.h"
#include "palette.h"

// The program only saves images and never opens a window, so CImg does not need X11
#define cimg_display 0
// Save PNG files with libpng instead of calling an external converter
#define cimg_use_png
#include "CImg.h"

using namespace cimg_library;
using namespace std;

void save_image(const iteration_grid &grid, const app_state &state, const string &path) {
	color body_color;
	find_color(state.colors.body_color, body_color);
	palette pal(state.fractal.palette_size, state.colors.colors, body_color, state.colors.skew);
	CImg<unsigned char> img(grid.width, grid.height, 1, 3, 0);
	for (size_t y_i = 0; y_i < grid.height; y_i++) {
		for (size_t x_i = 0; x_i < grid.width; x_i++) {
			color c = pal.at(grid.at(x_i, y_i));
			img(x_i, y_i, 0) = c.r;
			img(x_i, y_i, 1) = c.g;
			img(x_i, y_i, 2) = c.b;
		}
	}
	try {
		img.save(path.c_str());
	} catch (const CImgException &e) {
		cerr << "Could not save " << path << ": " << e.what() << endl;
	}
}

int main() {
	print_help(cout);
	app_state state;
	apply_preset("julia", state.fractal);

	iteration_grid grid;
	// draw the initial image before reading the first command
	outcome next = outcome::recompute;
	while (next != outcome::quit) {
		if (next == outcome::recompute) {
			grid = compute_fractal(state.fractal, state.threads_num);
		}
		if (next == outcome::recompute || next == outcome::recolor) {
			save_image(grid, state, "out.png");
			print_state(cout, state, {grid.width, grid.height});
		}

		cout << "> ";
		string line;
		// end of input, e.g. Ctrl+D or the end of a piped script, quits as well
		next = getline(cin, line) ? run_command_line(state, line, cerr) : outcome::quit;
	}
	return 0;
}
