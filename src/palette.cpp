#include "palette.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

using namespace std;

namespace {

const map<string, color> COLORS = {
	{"black", {0, 0, 0}},
	{"white", {255, 255, 255}},
	{"red", {255, 0, 0}},
	{"green", {0, 255, 0}},
	{"blue", {0, 0, 255}},
	{"yellow", {255, 255, 0}},
	{"magenta", {255, 0, 255}},
	{"cyan", {0, 255, 255}},
};

vector<string> split(const string &s, char delimiter) {
	vector<string> parts;
	stringstream ss(s);
	string part;
	while (getline(ss, part, delimiter)) {
		parts.push_back(part);
	}
	return parts;
}

}

color interpolate(color color1, color color2, double fraction) {
	return {
		(unsigned char) (color1.r * (1 - fraction) + color2.r * fraction),
		(unsigned char) (color1.g * (1 - fraction) + color2.g * fraction),
		(unsigned char) (color1.b * (1 - fraction) + color2.b * fraction)
	};
}

vector<string> color_names() {
	vector<string> names;
	for (const auto &entry : COLORS) {
		names.push_back(entry.first);
	}
	return names;
}

bool find_color(const string &name, color &result) {
	auto entry = COLORS.find(name);
	if (entry == COLORS.end()) {
		return false;
	}
	result = entry->second;
	return true;
}

bool is_color_sequence(const string &sequence) {
	vector<string> names = split(sequence, ':');
	return names.size() >= 2 && all_of(names.begin(), names.end(),
			[](const string &name) { return COLORS.count(name) > 0; });
}

palette::palette(size_t size, const string &color_sequence, color body_color, double skew) {
	colors.push_back(body_color);
	vector<color> gradient;
	for (const string &name : split(color_sequence, ':')) {
		gradient.push_back(COLORS.at(name));
	}
	double interval_length = 1. / (gradient.size() - 1);
	for (size_t i = 1; i < size; i++) {
		double fraction = pow((double) i / (size - 1), skew);
		for (size_t j = 1; j < gradient.size(); j++) {
			double border = j * interval_length;
			// the last border may be slightly below 1 due to rounding
			if (fraction <= border || j == gradient.size() - 1) {
				colors.push_back(interpolate(gradient[j - 1], gradient[j],
						1 - (border - fraction) / interval_length));
				break;
			}
		}
	}
}

color palette::at(double index) const {
	size_t int_index = (size_t) floor(index);
	return interpolate(colors[int_index], colors[int_index + 1], index - int_index);
}
