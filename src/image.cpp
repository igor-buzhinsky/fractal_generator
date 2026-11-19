#include "image.h"

#include "palette.h"

// The program only saves images and never opens a window, so CImg does not need X11
#define cimg_display 0
// Save PNG files with libpng instead of calling an external converter
#define cimg_use_png
#include <CImg.h>

using namespace cimg_library;
using namespace std;

bool save_image(const iteration_grid &grid, const coloring &colors, int palette_size,
		const string &path, ostream &errors) {
	color body_color;
	find_color(colors.body_color, body_color);
	palette pal((size_t) palette_size, colors.colors, body_color, colors.skew);
	// image sides are at most MAX_IMAGE_SIDE, so they fit into CImg's unsigned int
	const auto width = (unsigned) grid.width, height = (unsigned) grid.height;
	CImg<unsigned char> img(width, height, 1, 3, 0);
	for (unsigned y_i = 0; y_i < height; y_i++) {
		for (unsigned x_i = 0; x_i < width; x_i++) {
			color c = pal.at(grid.at(x_i, y_i));
			img(x_i, y_i, 0) = c.r;
			img(x_i, y_i, 1) = c.g;
			img(x_i, y_i, 2) = c.b;
		}
	}
	try {
		// report errors only through the exception, without CImg printing them as well
		cimg::exception_mode(0);
		img.save(path.c_str());
		return true;
	} catch (const CImgException &e) {
		errors << "Could not save " << path << ": " << e.what() << endl;
		return false;
	}
}
