// The interactive session: reads commands, computes the fractal and saves the image.

#pragma once

#include <functional>
#include <iostream>
#include <string>

#include "cli.h"
#include "fractal.h"

// Computes the iteration counts with the given number of threads, like compute_fractal()
using compute_function = std::function<iteration_grid(const fractal_settings &, int threads_num)>;

// Draws the initial image, then runs commands read line by line until 'exit' or the end
// of input. After each command that changes the image, saves it to image_path and prints
// the state. The image is computed with compute, which tests can replace.
void run_session(app_state &state, std::istream &in, std::ostream &out, std::ostream &errors,
		const std::string &image_path, const compute_function &compute = compute_fractal);
