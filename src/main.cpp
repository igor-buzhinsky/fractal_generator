//============================================================================
// Name        : fractal_generator
// Author      : Igor Buzhinsky
//============================================================================

#include <iostream>

#include "cli.h"
#include "fractal.h"
#include "session.h"

using namespace std;

int main() {
	print_help(cout);
	app_state state;
	apply_preset("julia", state.fractal);
	run_session(state, cin, cout, cerr, "out.png");
	return 0;
}
