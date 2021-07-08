#include <fmt/core.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "array2d.h"
#include "canvas.h"
#include "helpers.h"
#include "options.h"
#include "ray.h"
#include "sphere.h"
#include "vec.h"
#include "vec_utils.h"

// make sure canvas_width / canvas_height is divisible by 2

Color trace_ray(const Ray& ray, const std::vector<Sphere>& scene_spheres,
                const SceneOptions& options) {
  /*
   * Traces ray from min_z to max_z, computes the intersection with each
   * sphere and returns the color of the sphere closest to the min_z.
   */
  double t_min = options.t_min;
  double t_max = options.t_max;
  double closest_point = std::numeric_limits<double>::infinity();
  std::optional<Sphere> closest_sphere{};

  for (auto&& sphere : scene_spheres) {
    const auto& [intersect_optional_one, intersect_optional_two] =
        sphere.intersects_with(ray);

    if (intersect_optional_one) {
      double point_one = intersect_optional_one.value();
      if (in_range_inclusive(point_one, t_min, t_max) &&
          point_one < closest_point) {
        closest_point = point_one;
        closest_sphere = sphere;
      }
    }

    if (intersect_optional_two) {
      double point_two = intersect_optional_two.value();
      if (in_range_inclusive(point_two, t_min, t_max) &&
          point_two < closest_point) {
        closest_point = point_two;
        closest_sphere = sphere;
      }
    }
  }

  if (!closest_sphere) {
    return options.background_color;
  } else {
    return closest_sphere.value().color;
  }
}

/*
 * canvas_x and canvas_y are pixel coordinates in canvas space.
 * The canvas is organized so that the origin (0,0) point is in
 * the middle and not the top left.
 */
Vec3 canvas_to_viewport(i64 canvas_x, i64 canvas_y,
                        const SceneOptions& scene_options) {
  double viewport_x = (double)canvas_x * (scene_options.viewport_width /
                                          (double)scene_options.canvas_width);

  double viewport_y = (double)canvas_y * (scene_options.viewport_height /
                                          (double)scene_options.canvas_height);

  return {viewport_x, viewport_y, scene_options.dist_to_proj_plane};
}

int main() {
  SceneOptions options;
  options.origin = {0.0, 0.0, 0.0};

  Sphere red(1, Point3{0.0, -1, 3}, Color{255, 0, 0});
  Sphere blue(1, Point3{2, 0.0, 4.0}, Color{0, 0, 255});
  Sphere green(1, Point3{-2, 0.0, 4.0}, Color{0, 255, 0});

  std::vector<Sphere> spheres{red, green, blue};

  Color white(255, 255, 255);
  options.background_color = white;
  options.canvas_width = 400;
  options.canvas_height = 400;

  i64 width_start = -static_cast<i64>(options.canvas_width / 2);
  i64 width_end = static_cast<i64>(options.canvas_width / 2);

  i64 height_start = -static_cast<i64>(options.canvas_height / 2);
  i64 height_end = static_cast<i64>(options.canvas_height / 2);

  Canvas canvas(options.canvas_width, options.canvas_height);
  Ray ray(options.origin, Vec3{0, 0, 0});

  assert((options.canvas_height % 2 == 0 && options.canvas_width % 2 == 0) &&
         "canvas width and height should be even");

  for (i64 x = width_start; x < width_end; x++) {
    for (i64 y = height_start; y < height_end; y++) {
      Vec3 ray_direction = canvas_to_viewport(x, y, options);
      ray.set_ray_direction(ray_direction);
      Color pixel_color = trace_ray(ray, spheres, options);
      canvas.put_pixel(x, y, pixel_color);
    }
  }

  canvas.write_to_file();
}
