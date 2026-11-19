// Coloring iteration counts and saving the result as a PNG image.

#pragma once

#include <ostream>
#include <string>

#include "cli.h"
#include "fractal.h"

// Colors the iteration counts with the palette given by the coloring and saves them
// to a PNG file. Returns false and prints the reason to errors if saving fails.
bool save_image(const iteration_grid &grid, const coloring &colors, int palette_size,
		const std::string &path, std::ostream &errors);
