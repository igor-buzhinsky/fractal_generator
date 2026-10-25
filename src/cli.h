// Interactive commands: the command table, argument parsing, help and state output.

#pragma once

#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "fractal.h"

struct coloring {
	std::string colors = "black:red:yellow:green:blue:cyan:magenta:white";
	std::string body_color = "black";
	double skew = 1.0;
};

struct app_state {
	fractal_settings fractal;
	coloring colors;
	int threads_num = 4;
};

// What has to be done after a command
enum class outcome {
	nothing,
	recolor,   // redraw the image from the existing iteration counts
	recompute, // compute the fractal again and redraw the image
	quit
};

struct command {
	std::string section;            // help section the command is listed in
	std::vector<std::string> names; // several names share the description and the handler
	std::vector<std::string> arguments; // for help, e.g. "<double>"
	std::string description;        // '\n' starts a new paragraph
	// Throws std::invalid_argument if an argument is invalid; the state is then unchanged
	std::function<outcome(app_state &state, const std::string &name,
			const std::vector<std::string> &args)> run;
};

const std::vector<command> &commands();

// Returns nullptr if there is no such command
const command *find_command(const std::string &name);

// Splits an input line into a command name and its arguments, and runs the command.
// Prints errors to the given stream.
outcome run_command_line(app_state &state, const std::string &line, std::ostream &errors);

void print_help(std::ostream &out);

void print_state(std::ostream &out, const app_state &state, const image_size &size);

// Parses the whole text as a finite number
bool parse_double(const std::string &text, double &value);

// Parses the whole text as an integer that fits into int
bool parse_int(const std::string &text, int &value);

// Formats a number with the given number of significant digits, without trailing zeros
std::string format_number(double value, int precision = 15);
