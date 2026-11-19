#include "session.h"

#include "image.h"

using namespace std;

void run_session(app_state &state, istream &in, ostream &out, ostream &errors,
		const string &image_path, const compute_function &compute) {
	iteration_grid grid;
	// draw the initial image before reading the first command
	outcome next = outcome::recompute;
	while (next != outcome::quit) {
		if (next == outcome::recompute) {
			grid = compute(state.fractal, state.threads_num);
		}
		if (next == outcome::recompute || next == outcome::recolor) {
			save_image(grid, state.colors, state.fractal.palette_size, image_path, errors);
			print_state(out, state, {grid.width, grid.height});
		}

		out << "> ";
		string line;
		// end of input, e.g. Ctrl+D or the end of a piped script, quits as well
		next = getline(in, line) ? run_command_line(state, line, errors) : outcome::quit;
	}
}
