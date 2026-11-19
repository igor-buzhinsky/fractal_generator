//============================================================================
// Name        : fractal_generator
// Author      : Igor Buzhinsky
//============================================================================

#include <iostream>
#include <string>

#include "cli.h"
#include "fractal.h"
#include "image.h"

using namespace std;

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
			save_image(grid, state.colors, state.fractal.palette_size, "out.png", cerr);
			print_state(cout, state, {grid.width, grid.height});
		}

		cout << "> ";
		string line;
		// end of input, e.g. Ctrl+D or the end of a piped script, quits as well
		next = getline(cin, line) ? run_command_line(state, line, cerr) : outcome::quit;
	}
	return 0;
}
