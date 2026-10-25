// Colors and palettes that map iteration counts to colors.

#pragma once

#include <string>
#include <vector>

struct color {
	unsigned char r = 0, g = 0, b = 0;
};

// Linear interpolation: fraction 0 gives color1, fraction 1 gives color2
color interpolate(color color1, color color2, double fraction);

// Names of all available colors
std::vector<std::string> color_names();

// Looks up a color by name; returns false if there is no such color
bool find_color(const std::string &name, color &result);

// Checks a sequence like "black:red:white": at least two known colors
bool is_color_sequence(const std::string &sequence);

struct palette {
	// colors[0] is the body color, the rest is the gradient for escaping points
	std::vector<color> colors;

	// color_sequence must be valid for is_color_sequence(). The gradient is spread over
	// the palette with the given skew: fraction f of the palette gets gradient position f^skew.
	palette(size_t size, const std::string &color_sequence, color body_color, double skew);

	// Color for a possibly fractional index from iteration_count()
	color at(double index) const;
};
