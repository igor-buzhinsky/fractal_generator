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

# Runs the program with the given commands (one per line) and stores its output.
# Every run must exit successfully, even with invalid commands, and without sanitizer
# reports: leaks and data races are only reported when the program exits.
run() {
	rm -f out.png
	# generous timeout: sanitizer builds are much slower, it only guards against hangs
	output="$(printf '%s\n' "$@" | timeout 600 "$BINARY" 2>&1)"
	local status=$?
	# the state printed last, after the last command that changed the image
	last_state="$(awk '/^Formula:/ { block = "" } { block = block $0 "\n" } END { printf "%s", block }' \
			<<< "$output")"
	checks=$((checks + 1))
	[ "$status" -eq 0 ] || fail "'$*': exit code $status"
	checks=$((checks + 1))
	! grep -qE 'Sanitizer|runtime error:' <<< "$output" || fail "'$*': sanitizer report in output"
}

expect_output() {
	checks=$((checks + 1))
	grep -qF -- "$2" <<< "$output" || fail "$1: output does not contain '$2'"
}

# Checks the state printed last, so that a state printed before the commands cannot match
expect_state() {
	checks=$((checks + 1))
	grep -qF -- "$2" <<< "$last_state" || fail "$1: last printed state does not contain '$2'"
}

# Checks that out.png has the size given in the state printed last
expect_png_of_last_state() {
	local size
	size="$(sed -nE 's/^Image: +([0-9]+) x ([0-9]+) pixels.*/\1 \2/p' <<< "$last_state")"
	# shellcheck disable=SC2086 # split into width and height
	expect_png "$1" ${size:-0 0}
}

# Checks that out.png is a PNG image of the given size
expect_png() {
	checks=$((checks + 1))
	# the PNG signature and the image size in the IHDR chunk take 24 bytes
	if [ ! -f out.png ] || [ "$(wc -c < out.png)" -lt 24 ]; then
		fail "$1: out.png not created or truncated"
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
expect_output "default Julia set" "Presets:"
expect_output "default Julia set" "Julia set: z starts at the point"
expect_png "default Julia set" 123 81

# each preset has its own view, so the image size shows that the image was computed again
run "points 10000" "mandelbrot" "exit"
expect_state "mandelbrot" "z -> sqr(z) + c"
expect_state "mandelbrot" "Mandelbrot-like set"
expect_state "mandelbrot" "x in [-1.6, 1], y in [-1.05, 1.05]"
expect_png "mandelbrot" 111 89
run "points 10000" "tricorn" "exit"
expect_state "tricorn" "z -> sqr(conj(z)) + c"
expect_png "tricorn" 100 100
run "points 10000" "burning_ship" "exit"
expect_state "burning_ship" "z -> sqr(component_abs(z)) + conj(c)"
expect_png "burning_ship" 122 81
run "points 10000" "mandelbrot" "julia" "exit"
expect_state "julia" "Julia set: z starts at the point"
expect_png_of_last_state "julia"

run "mandelbrot" "x_min -2" "x_max 2" "y_min -1" "y_max 1" "points 20000" "exit"
expect_state "view bounds" "x in [-2, 2], y in [-1, 1]"
expect_png "view bounds" 200 100

run "mandelbrot" "x_min -2" "x_max 2" "y_min -2" "y_max 2" "zoom 4" "exit"
expect_state "zoom" "x in [-0.5, 0.5], y in [-0.5, 0.5]"

run "mandelbrot" "x_min -2" "x_max 2" "y_min -2" "y_max 2" "left 0.25" "top 0.5" "exit"
expect_state "left and top" "x in [-2, -1], y in [0, 2]"

run "mandelbrot" "x_min -2" "x_max 2" "y_min -2" "y_max 2" "right 0.25" "bottom 0.5" "exit"
expect_state "right and bottom" "x in [1, 2], y in [-2, 0]"

run "points 10000" "g_z1 conj" "g sin" "exit"
expect_state "custom formula" "z -> sin(sqr(conj(z)) + c)"

run "points 10000" "g_z1 conj" "g_z1 id" "exit"
expect_state "identity function" "z -> sqr(z) + c"

run "points 10000" "julia_c 0.3 -0.5" "exit"
expect_state "julia_c" "c = 0.3 - 0.5i"

# color changes come last: they must save the image and print the state without computing again
run "points 10000" "continuous_intensity 1" "exit"
expect_state "smooth coloring" "continuous_intensity = on"
plain_png="$(cksum < out.png)"
run "points 10000" "continuous_intensity 1" "colors black:white" "body_color red" "palette_skew 2" "exit"
expect_state "colors" "black:white (colors), inside the set: red (body_color)"
expect_state "colors" "palette_skew = 2, continuous_intensity = on"
expect_png_of_last_state "colors"
checks=$((checks + 1))
[ "$(cksum < out.png)" != "$plain_png" ] || fail "colors: out.png did not change"

run "points 10000" "threads_num 1" "plot" "exit"
expect_state "threads_num" "Threads:  1"

# 'plot' prints the state after the rejected commands
run "x_min 100" "y_max -5" "points 10" "palette_size 1" "zoom 0" "continuous_intensity 2" \
	"threads_num 0" "plot" "exit"
expect_output "invalid values" "x_min: must be less than x_max"
expect_output "invalid values" "y_max: must be greater than y_min"
expect_output "invalid values" "points: expected an integer from 100 to 100000000, got '10'"
expect_output "invalid values" "palette_size: expected an integer from 10 to 10000000, got '1'"
expect_output "invalid values" "zoom: the value must be positive"
expect_output "invalid values" "continuous_intensity: expected an integer from 0 to 1, got '2'"
expect_output "invalid values" "threads_num: expected an integer from 1 to 1024, got '0'"
# rejected values must not change the state
expect_state "invalid values" "x in [-1.6, 1.6], y in [-1.05, 1.05]"
expect_state "invalid values" "pixels (points = 300000)"
expect_state "invalid values" "palette_size = 1000, palette_skew = 1, continuous_intensity = off"
expect_state "invalid values" "Threads:  4"

run "zoom abc" "points 1.5" "julia_c 1 x" "left 0.5x" "exit"
expect_output "not numbers" "zoom: expected a number, got 'abc'"
expect_output "not numbers" "points: expected an integer from 100 to 100000000, got '1.5'"
expect_output "not numbers" "julia_c: expected a number, got 'x'"
expect_output "not numbers" "left: expected a number, got '0.5x'"

run "points 10000" "colors black" "colors black:purple" "body_color purple" "g_z1 no_such_function" \
	"plot" "exit"
expect_output "invalid names" "colors: expected at least two known colors separated by ':', got 'black'"
expect_output "invalid names" "colors: expected at least two known colors separated by ':', got 'black:purple'"
expect_output "invalid names" "body_color: unknown color 'purple'"
expect_output "invalid names" "g_z1: unknown function 'no_such_function'"
expect_state "invalid names" "black:red:yellow:green:blue:cyan:magenta:white (colors), inside the set: black (body_color)"
expect_state "invalid names" "z -> sqr(z) + c"

# a very wide view still produces a valid image
run "x_min -2" "x_max 2" "y_min 0" "y_max 0.000000001" "points 10000" "exit"
expect_png "very wide view" 5000 2

run "no_such_command" "exit"
expect_output "unknown command" "no_such_command: unknown command"

run "help" "exit"
[ "$(grep -c "Presets:" <<< "$output")" -eq 2 ] || fail "help: help text not printed again"
checks=$((checks + 1))

# The program must also quit when the input ends without 'exit'
run "points 10000"

echo "CLI tests: $((checks - failures))/$checks checks passed"
[ "$failures" -eq 0 ]
