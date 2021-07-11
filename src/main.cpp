#include <fmt/core.h>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "array2d.h"
#include "canvas.h"
#include "helpers.h"
#include "intdefines.h"
#include "options.h"
#include "ray.h"
#include "sphere.h"
#include "vec.h"
#include "vec_utils.h"

// TODO: be more consistent with what constructs need encapsulation
// and which do not.
//
// TODO: consider refactoring the lighting and color calculations
// into a scene class?

/*
 * For now, these classes only implement diffuse lighting.
 *
 */
class Light {
 public:
  virtual double compute_intensity_at(const Point3& point, const Vec3& normal,
                                      const Vec3& view, i64 specular) const = 0;
  virtual ~Light() {}
};

/*
 * Ambient light is a cheap way of
 * not considering secondary light rays.
 */
class AmbientLight : public Light {
 private:
  double intensity{};

 public:
  explicit AmbientLight(double intensity) : intensity(intensity) {}

  // this function is a victim of not fitting into the pattern.
  double compute_intensity_at(const Point3& point, const Vec3& normal,
                              const Vec3& view, i64 specular) const override {
    return intensity;
  }
};

class PointLight : public Light {
 private:
  double intensity{};
  Point3 position{};

 public:
  PointLight(double intensity, const Point3& position)
      : intensity(intensity), position(position) {}

  double compute_intensity_at(const Point3& point, const Vec3& normal,
                              const Vec3& view, i64 specular) const override {
    // direction_vector == L (from book)
    Vec3 direction_vector = position - point;
    double n_dot_l = dot(direction_vector, normal);
    double retval = 0.0;
    if (n_dot_l < 0) {
      return 0;
    } else {
      retval +=
          intensity * n_dot_l / (length(normal) * length(direction_vector));
    }

    if (specular != -1) {
      Vec3 reflected_ray = 2.0 * normal * n_dot_l - direction_vector;
      double r_dot_v = dot(reflected_ray, view);
      if (r_dot_v > 0) {
        retval += intensity *
                  std::pow(r_dot_v / (length(reflected_ray) * length(view)),
                           specular);
      }
    }

    return retval;
  }
};

/* A directional light has the same direction
 * at all points in the scene. This represents
 * a light source which is a very large distance
 * away (e.g. like the sun where basically get
 * parallel light vectors).
 */
class DirectionalLight : public Light {
 private:
  double intensity;
  Vec3 direction{};

 public:
  DirectionalLight(double intensity, const Vec3& direction)
      : intensity(intensity), direction(direction) {}


  double compute_intensity_at(const Point3& point, const Vec3& normal,
                              const Vec3& view, i64 specular) const override {
    double n_dot_l = dot(direction, normal);
    double retval = 0.0;
    if (n_dot_l < 0) {
      return 0;
    } else {
      retval += intensity * n_dot_l / (length(normal) * length(direction));
    }

    if (specular != -1) {
      Vec3 reflected_ray = 2.0 * normal * n_dot_l - direction;
      double r_dot_v = dot(reflected_ray, view);
      if (r_dot_v > 0) {
        retval += intensity *
                  std::pow(r_dot_v / (length(reflected_ray) * length(view)),
                           specular);
      }
    }

    return retval;
  }
};

/*
 * For each point in the scene, computes the total intensity
 * of the light.
 * @param point: point to compute lighting at
 * @param normal: surface normal at "point"
 * @param view: vector from point "point" to the camera.
 * @param specular: specular coefficient at point (depends on material)
 * @param lights: vector of all lights in the scene
 */
double compute_lighting(const Vec3& point, const Vec3& normal, const Vec3& view,
                        i64 specular,
                        const std::vector<std::unique_ptr<Light>>& lights) {
  double intensity = 0.0;
  for (auto&& light : lights) {
    intensity += light->compute_intensity_at(point, normal, view, specular);
  }

  return intensity;
}

// make sure canvas_width / canvas_height is divisible by 2

Color trace_ray(const Ray& ray, const std::vector<Sphere>& scene_spheres,
                const std::vector<std::unique_ptr<Light>>& lights,
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
    //
    Point3 point_on_surface =
        ray.get_origin() + closest_point * ray.get_direction();
    Vec3 surface_normal = point_on_surface - closest_sphere.value().center;
    surface_normal = surface_normal / length(surface_normal);

    return closest_sphere.value().color *
           compute_lighting(point_on_surface, surface_normal,
                            -ray.get_direction(),
                            closest_sphere.value().specular, lights);
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

  Sphere red(1, 500, Point3{0.0, -1, 3}, Color{255, 0, 0});
  Sphere blue(1, 500, Point3{2, 0.0, 4.0}, Color{0, 0, 255});
  Sphere green(1, 10, Point3{-2, 0.0, 4.0}, Color{0, 255, 0});
  Sphere large_yellow(5000, 1000, Point3{0, -5001, 0}, Color{255, 255, 0});

  std::vector<Sphere> spheres{red, green, blue, large_yellow};

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

  // TODO: these need not be dynamically allocated, we could just store
  // references to them in the vector and create them e.g.
  // in a Scene class or just directly in main.
  std::vector<std::unique_ptr<Light>> lights{};
  lights.reserve(3);

  lights.push_back(std::make_unique<AmbientLight>(0.2));
  lights.push_back(std::make_unique<PointLight>(0.6, Point3{2, 1, 0}));
  lights.push_back(std::make_unique<DirectionalLight>(0.2, Vec3{1, 4, 4}));

  assert((options.canvas_height % 2 == 0 && options.canvas_width % 2 == 0) &&
         "canvas width and height should be even");

  for (i64 x = width_start; x < width_end; x++) {
    for (i64 y = height_start; y < height_end; y++) {
      Vec3 ray_direction = canvas_to_viewport(x, y, options);
      ray.set_ray_direction(ray_direction);
      Color pixel_color = trace_ray(ray, spheres, lights, options);
      canvas.put_pixel(x, y, pixel_color);
    }
  }

  canvas.write_to_file();
}
