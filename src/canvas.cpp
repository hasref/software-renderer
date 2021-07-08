#include "canvas.h"

#include <fmt/core.h>

#include <fstream>

#include "array2d.h"
#include "intdefines.h"
#include "vec.h"

Canvas::Canvas(usize width, usize height)
    : canvas_width(width), canvas_height(height), screen(height, width) {}

const auto& Canvas::get_canvas() const { return screen; }

usize Canvas::get_width() const { return canvas_width; }
usize Canvas::get_height() const { return canvas_height; }

void Canvas::put_pixel(i64 x, i64 y, const Color& color) {
  i64 width_start = -static_cast<i64>(canvas_width / 2);
  i64 width_end = static_cast<i64>(canvas_width / 2);

  i64 height_start = -static_cast<i64>(canvas_height / 2);
  i64 height_end = static_cast<i64>(canvas_height / 2);

  if (!((x >= width_start) && (x < width_end)) ||
      !(y >= height_start && y < height_end)) {
    fmt::print("Invalid write to pixel location {}, {}\n", x, y);
    // with the PpmWriter, we can still write the value to file
    // although the image may be a little garbled in this case.
  }

  // this conversion is fine - want it to happen
  usize screen_x = x + (canvas_width / 2);

  // reverse y range for screen coordinates
  usize screen_y = (canvas_height - 1) - (y + (canvas_height / 2));

  screen(screen_y, screen_x) = color;
}

void Canvas::write_to_file(const std::string& filename) {
  std::ofstream file(filename);
  if (!file) {
    fmt::print(stderr, "Failed to open file {}", filename);
    return;
  }
  file << "P3\n";
  file << canvas_width << " " << canvas_height << "\n";
  file << 255 << "\n";

  for (usize row = 0; row < canvas_height; row++) {
    for (usize col = 0; col < canvas_width; col++) {
      Color color = screen(row, col);
      file << +color.red << " " << +color.green << " " << +color.blue << "\n";
    }
  }
}
