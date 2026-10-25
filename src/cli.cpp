#include "cli.h"

#include <climits>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>

#include "palette.h"

using namespace std;

namespace {

double read_double(const string &text) {
	double value;
	if (!parse_double(text, value)) {
		throw invalid_argument("expected a number, got '" + text + "'");
	}
	return value;
}

double read_positive_double(const string &text) {
	double value = read_double(text);
	if (value <= 0) {
		throw invalid_argument("the value must be positive");
	}
	return value;
}

int read_int(const string &text, int min_value, int max_value) {
	int value;
	if (!parse_int(text, value) || value < min_value || value > max_value) {
		throw invalid_argument("expected an integer from " + to_string(min_value) + " to "
				+ to_string(max_value) + ", got '" + text + "'");
	}
	return value;
}

string join(const vector<string> &parts, const string &separator) {
	string result;
	for (size_t i = 0; i < parts.size(); i++) {
		result += (i == 0 ? "" : separator) + parts[i];
	}
	return result;
}

// Splits text into lines of at most the given width; '\n' always starts a new line
vector<string> wrap(const string &text, size_t width) {
	vector<string> lines;
	stringstream paragraphs(text);
	string paragraph;
	while (getline(paragraphs, paragraph)) {
		stringstream words(paragraph);
		string word, line;
		while (words >> word) {
			if (!line.empty() && line.size() + 1 + word.size() > width) {
				lines.push_back(line);
				line.clear();
			}
			line += (line.empty() ? "" : " ") + word;
		}
		lines.push_back(line);
	}
	return lines;
}

const string FORMULA_SECTION = "Formula f_c(z) = g(g_z2(g_z1(z)) + g_c(c))";

vector<command> make_commands() {
	return {
		{"Presets", preset_names(), {}, "Load a fractal with its default view",
			[](app_state &state, const string &name, const vector<string> &) {
				apply_preset(name, state.fractal);
				return outcome::recompute;
			}},
		{"View", {"x_min", "x_max", "y_min", "y_max"}, {"<double>"}, "Set a boundary of the visible area",
			[](app_state &state, const string &name, const vector<string> &args) {
				double value = read_double(args[0]);
				view &area = state.fractal.area;
				bool is_x = name[0] == 'x', is_min = name.substr(2) == "min";
				double &bound = is_x ? (is_min ? area.x_min : area.x_max) : (is_min ? area.y_min : area.y_max);
				double other = is_x ? (is_min ? area.x_max : area.x_min) : (is_min ? area.y_max : area.y_min);
				if (is_min ? value >= other : value <= other) {
					throw invalid_argument(string("must be ") + (is_min ? "less than " : "greater than ")
							+ name.substr(0, 2) + (is_min ? "max" : "min"));
				}
				bound = value;
				return outcome::recompute;
			}},
		{"View", {"zoom"}, {"<double>"}, "Zoom in around the center (< 1 zooms out)",
			[](app_state &state, const string &, const vector<string> &args) {
				state.fractal.area.zoom(read_positive_double(args[0]));
				return outcome::recompute;
			}},
		{"View", {"left", "right", "top", "bottom"}, {"<double>"},
			"Keep only this fraction of the view on that side, e.g. 'left 0.5' keeps the left half",
			[](app_state &state, const string &name, const vector<string> &args) {
				double removed = 1 - read_positive_double(args[0]);
				view &area = state.fractal.area;
				if (name == "left") {
					area.x_max -= area.width() * removed;
				} else if (name == "right") {
					area.x_min += area.width() * removed;
				} else if (name == "top") {
					area.y_min += area.height() * removed;
				} else {
					area.y_max -= area.height() * removed;
				}
				return outcome::recompute;
			}},
		{"Image", {"points"}, {"<integer>"}, "Total number of pixels in the output image",
			[](app_state &state, const string &, const vector<string> &args) {
				state.fractal.pixel_count = read_int(args[0], 100, 100000000);
				return outcome::recompute;
			}},
		{"Image", {"palette_size"}, {"<integer>"}, "Number of colors, also the maximum iteration count",
			[](app_state &state, const string &, const vector<string> &args) {
				state.fractal.palette_size = read_int(args[0], 10, 10000000);
				return outcome::recompute;
			}},
		{"Image", {"palette_skew"}, {"<double>"}, "Exponent of color distribution (1 = uniform)",
			[](app_state &state, const string &, const vector<string> &args) {
				state.colors.skew = read_positive_double(args[0]);
				return outcome::recolor;
			}},
		{"Image", {"colors"}, {"<color_1:...:color_n>"},
			"Gradient (2+ colors) for points outside the set; colors: " + join(color_names(), ", "),
			[](app_state &state, const string &, const vector<string> &args) {
				if (!is_color_sequence(args[0])) {
					throw invalid_argument("expected at least two known colors separated by ':', got '"
							+ args[0] + "'");
				}
				state.colors.colors = args[0];
				return outcome::recolor;
			}},
		{"Image", {"body_color"}, {"<color>"}, "Color of points inside the set",
			[](app_state &state, const string &, const vector<string> &args) {
				color unused;
				if (!find_color(args[0], unused)) {
					throw invalid_argument("unknown color '" + args[0] + "'");
				}
				state.colors.body_color = args[0];
				return outcome::recolor;
			}},
		{"Image", {"continuous_intensity"}, {"<0 or 1>"}, "Smooth coloring instead of discrete bands",
			[](app_state &state, const string &, const vector<string> &args) {
				state.fractal.continuous_intensity = read_int(args[0], 0, 1);
				return outcome::recompute;
			}},
		{FORMULA_SECTION, {"g", "g_z1", "g_z2", "g_c"}, {"<function_name>"},
			"Replace one of the functions in the formula; available: " + join(function_names(), ", "),
			[](app_state &state, const string &name, const vector<string> &args) {
				named_function function;
				if (!find_function(args[0], function)) {
					throw invalid_argument("unknown function '" + args[0] + "'");
				}
				formula &f = state.fractal.f;
				(name == "g" ? f.g : name == "g_z1" ? f.g_z1 : name == "g_z2" ? f.g_z2 : f.g_c) = function;
				return outcome::recompute;
			}},
		{FORMULA_SECTION, {"julia_on"}, {"<0 or 1>"},
			"1: Julia set (c is fixed, z_0 is the point),\n0: Mandelbrot-like set (z_0 = 0, c is the point)",
			[](app_state &state, const string &, const vector<string> &args) {
				state.fractal.julia_on = read_int(args[0], 0, 1);
				return outcome::recompute;
			}},
		{FORMULA_SECTION, {"julia_c"}, {"<re:double>", "<im:double>"}, "Constant c used for Julia sets",
			[](app_state &state, const string &, const vector<string> &args) {
				double re = read_double(args[0]), im = read_double(args[1]);
				state.fractal.julia_c = dcomplex(re, im);
				return outcome::recompute;
			}},
		{"Other", {"threads_num"}, {"<integer>"}, "Number of threads used for computation",
			[](app_state &state, const string &, const vector<string> &args) {
				state.threads_num = read_int(args[0], 1, 1024);
				return outcome::nothing;
			}},
		{"Other", {"plot"}, {}, "Redraw the image",
			[](app_state &, const string &, const vector<string> &) {
				return outcome::recompute;
			}},
		{"Other", {"help"}, {}, "Show this list",
			[](app_state &, const string &, const vector<string> &) {
				print_help(cout);
				return outcome::nothing;
			}},
		{"Other", {"exit"}, {}, "Quit the program",
			[](app_state &, const string &, const vector<string> &) {
				return outcome::quit;
			}},
	};
}

}

const vector<command> &commands() {
	static const vector<command> all_commands = make_commands();
	return all_commands;
}

const command *find_command(const string &name) {
	for (const command &c : commands()) {
		for (const string &command_name : c.names) {
			if (command_name == name) {
				return &c;
			}
		}
	}
	return nullptr;
}

outcome run_command_line(app_state &state, const string &line, ostream &errors) {
	stringstream words(line);
	string name, word;
	vector<string> args;
	if (!(words >> name)) {
		return outcome::nothing;
	}
	while (words >> word) {
		args.push_back(word);
	}

	const command *c = find_command(name);
	if (c == nullptr) {
		errors << name << ": unknown command, type 'help' to see the list of available ones" << endl;
		return outcome::nothing;
	}
	if (args.size() != c->arguments.size()) {
		errors << name << ": expected " << c->arguments.size() << " argument(s): "
				<< join(c->arguments, " ") << endl;
		return outcome::nothing;
	}
	try {
		return c->run(state, name, args);
	} catch (const invalid_argument &e) {
		errors << name << ": " << e.what() << endl;
		return outcome::nothing;
	}
}

void print_help(ostream &out) {
	const size_t description_column = 37, max_width = 93;
	string section;
	for (const command &c : commands()) {
		if (c.section != section) {
			section = c.section;
			out << section << ":" << endl;
		}
		string usage = "    " + join(c.names, " | ") + (c.arguments.empty() ? "" : " ")
				+ join(c.arguments, " ");
		vector<string> lines = wrap(c.description, max_width - description_column);
		if (usage.size() < description_column) {
			out << usage << string(description_column - usage.size(), ' ') << lines[0] << endl;
			lines.erase(lines.begin());
		} else {
			out << usage << endl;
		}
		for (const string &line : lines) {
			out << string(description_column, ' ') << line << endl;
		}
	}
	out << "Each new image is saved to out.png." << endl;
}

void print_state(ostream &out, const app_state &state, const image_size &size) {
	const fractal_settings &s = state.fractal;
	const view &area = s.area;
	out << endl;
	out << "Formula:  z -> " << s.f.to_string() << endl;
	if (s.julia_on) {
		out << "          Julia set: z starts at the point, c = " << format_number(real(s.julia_c))
				<< (imag(s.julia_c) < 0 ? " - " : " + ") << format_number(fabs(imag(s.julia_c)))
				<< "i (julia_c)" << endl;
	} else {
		out << "          Mandelbrot-like set: z starts at 0, c is the point" << endl;
	}
	out << "View:     x in [" << format_number(area.x_min) << ", " << format_number(area.x_max)
			<< "], y in [" << format_number(area.y_min) << ", " << format_number(area.y_max) << "]" << endl;
	out << "          center (" << format_number((area.x_min + area.x_max) / 2) << ", "
			<< format_number((area.y_min + area.y_max) / 2) << "), width " << format_number(area.width(), 6)
			<< ", height " << format_number(area.height(), 6) << endl;
	out << "Image:    " << size.width << " x " << size.height << " pixels (points = " << s.pixel_count
			<< "), saved to out.png" << endl;
	out << "Colors:   " << state.colors.colors << " (colors), inside the set: " << state.colors.body_color
			<< " (body_color)" << endl;
	out << "          palette_size = " << s.palette_size << ", palette_skew = "
			<< format_number(state.colors.skew) << ", continuous_intensity = "
			<< (s.continuous_intensity ? "on" : "off") << endl;
	out << "Threads:  " << state.threads_num << endl;
}

bool parse_double(const string &text, double &value) {
	char *end;
	value = strtod(text.c_str(), &end);
	return !text.empty() && *end == '\0' && isfinite(value);
}

bool parse_int(const string &text, int &value) {
	char *end;
	long result = strtol(text.c_str(), &end, 10);
	value = (int) result;
	return !text.empty() && *end == '\0' && result >= INT_MIN && result <= INT_MAX;
}

string format_number(double value, int precision) {
	char result[64];
	snprintf(result, sizeof(result), "%.*g", precision, value);
	return result;
}
