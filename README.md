# Fractal generator

This is a simple image generator for Mandelbrot and Julia sets and some derived fractals. It has an interactive command-line interface (e.g. for zooming) and generates its output as an image in the same directory.

Build on Linux:

> make

The build requires a C++17 compiler, and saving PNG images requires ImageMagick. For example, on Ubuntu:

> sudo apt install build-essential imagemagick

Run the tests:

> make test

This runs unit tests for the computation code (`tests/unit_tests.cpp`) and end-to-end tests that feed commands to the program and check its output and the generated image (`tests/cli_tests.sh`).

Run:

> ./fractal_generator

The tool will print its list of commands.

Example to plot the Mandelbrot set and then zoom the top left corner:

```
mandelbrot
top 0.5
left 0.5
```
