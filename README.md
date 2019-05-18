# Fractal generator

This is a simple image generator for Mandelbrot and Julia sets and some derived fractals. It has an interactive command-line interface (e.g. for zooming) and generates its output as an image in the same directory.

Build on Linux:

> make

Run:

> ./fractal_generator

The tool will print its list of commands.

Example to plot the Mandelbrot set and then zoom the top left corner:

> mandelbrot
> top 0.5
> left 0.5
