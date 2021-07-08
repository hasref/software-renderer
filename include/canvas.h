#ifndef CANVAS_H_
#define CANVAS_H_

#include <fmt/core.h>

#include <fstream>

#include "array2d.h"
#include "intdefines.h"
#include "vec.h"

class Canvas {
 private:
  usize canvas_width{};
  usize canvas_height{};
  Array2d<Color> screen{800, 800};

 public:
  Canvas() = default;

  Canvas(usize width, usize height);

  const auto& get_canvas() const;

  usize get_width() const;
  usize get_height() const;

  /*
   * Writes pixel color at location x,y of the canvas. Note that x belongs
   * to
   * [-canvas_width / 2, canvas_width / 2), while y belongs to
   * [-canvas_height / 2, canvas_height / 2).
   *
   * @param x: horizontal canvas pixel in [-canvas_width / 2, canvas_width
   * / 2)
   * @param y: vertical canvas pixel in [-canvas_height / 2, canvas_height
   * / 2)
   */
  void put_pixel(i64 x, i64 y, const Color& color);

  void write_to_file(const std::string& filename = "output.ppm");
};

#endif
