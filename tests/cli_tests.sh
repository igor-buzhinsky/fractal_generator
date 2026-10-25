#!/usr/bin/env bash
# End-to-end tests: feed commands to the program and check its output and out.png.
# Run with: make test

BINARY="$(cd "$(dirname "$0")/.." && pwd)/fractal_generator"
WORK_DIR="$(mktemp -d)"
trap 'rm -rf "$WORK_DIR"' EXIT
cd "$WORK_DIR" || exit 1

checks=0
failures=0

fail() {
	echo "FAIL: $1"
	failures=$((failures + 1))
}

# Runs the program with the given commands (one per line) and stores its output
run() {
	rm -f out.png
	output="$(printf '%s\n' "$@" | timeout 60 "$BINARY" 2>&1)"
	status=$?
}

expect_success() {
	checks=$((checks + 1))
	[ "$status" -eq 0 ] || fail "$1: exit code $status"
}

expect_output() {
	checks=$((checks + 1))
	grep -qF -- "$2" <<< "$output" || fail "$1: output does not contain '$2'"
}

# Checks that out.png is a PNG image of the given size
expect_png() {
	checks=$((checks + 1))
	if [ ! -f out.png ]; then
		fail "$1: out.png not created"
		return
	fi
	local header width height
	header="$(od -An -tx1 -N8 out.png | tr -d ' \n')"
	width=$(( 16#$(od -An -tx1 -j16 -N4 out.png | tr -d ' \n') ))
	height=$(( 16#$(od -An -tx1 -j20 -N4 out.png | tr -d ' \n') ))
	[ "$header" = "89504e470d0a1a0a" ] || fail "$1: out.png is not a PNG file"
	[ "$width" -eq "$2" ] && [ "$height" -eq "$3" ] \
		|| fail "$1: out.png is ${width}x${height}, expected $2x$3"
}

run "points 10000" "exit"
expect_success "default Julia set"
expect_output "default Julia set" "Presets:"
expect_output "default Julia set" "Julia set: z starts at the point"
expect_png "default Julia set" 123 81

for preset in mandelbrot tricorn burning_ship julia; do
	run "points 10000" "$preset" "exit"
	expect_success "$preset"
	checks=$((checks + 1))
	[ -f out.png ] || fail "$preset: out.png not created"
done

run "points 10000" "mandelbrot" "exit"
expect_output "mandelbrot formula" "z -> sqr(z) + c"
expect_output "mandelbrot mode" "Mandelbrot-like set"

run "mandelbrot" "x_min -2" "x_max 2" "y_min -1" "y_max 1" "points 20000" "exit"
expect_output "view bounds" "x in [-2, 2], y in [-1, 1]"
expect_png "view bounds" 200 100

run "mandelbrot" "x_min -2" "x_max 2" "y_min -2" "y_max 2" "zoom 4" "exit"
expect_output "zoom" "x in [-0.5, 0.5], y in [-0.5, 0.5]"

run "mandelbrot" "x_min -2" "x_max 2" "y_min -2" "y_max 2" "left 0.25" "top 0.5" "exit"
expect_output "left and top" "x in [-2, -1], y in [-2, 0]"

run "mandelbrot" "x_min -2" "x_max 2" "y_min -2" "y_max 2" "right 0.25" "bottom 0.5" "exit"
expect_output "right and bottom" "x in [1, 2], y in [0, 2]"

run "points 10000" "g_z1 conj" "g sin" "exit"
expect_output "custom formula" "z -> sin(sqr(conj(z)) + c)"

run "points 10000" "julia_c 0.3 -0.5" "exit"
expect_output "julia_c" "c = 0.3 - 0.5i"

run "points 10000" "colors black:white" "body_color red" "palette_skew 2" "continuous_intensity 1" "exit"
expect_output "colors" "black:white (colors), inside the set: red (body_color)"
expect_output "colors" "palette_skew = 2, continuous_intensity = on"

run "points 10000" "threads_num 1" "plot" "exit"
expect_output "threads_num" "Threads:  1"

run "x_min 100" "points 10" "palette_size 1" "zoom 0" "continuous_intensity 2" "threads_num 0" "exit"
expect_success "invalid values"
expect_output "invalid values" "x_min must be less than x_max!"
expect_output "invalid values" "The value must be at least 100!"
expect_output "invalid values" "The value must be at least 10!"
expect_output "invalid values" "The value must be positive!"
expect_output "invalid values" "Only valid values are 0 and 1!"
expect_output "invalid values" "The value must be between 1 and 1024!"

run "no_such_command" "exit"
expect_output "unknown command" "Unknown command!"

run "help" "exit"
expect_success "help"
[ "$(grep -c "Presets:" <<< "$output")" -eq 2 ] || fail "help: help text not printed again"
checks=$((checks + 1))

# The program must also quit when the input ends without 'exit'
run "points 10000"
expect_success "end of input"

echo "CLI tests: $((checks - failures))/$checks checks passed"
[ "$failures" -eq 0 ]
