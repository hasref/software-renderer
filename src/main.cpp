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

// make sure canvas_width / canvas_height is divisible by 2

struct Vec2 {
  double x = 0;
  double y = 0;

  Vec2() = default;
  Vec2(double x, double y) : x(x), y(y) {}
};

using Point2 = Vec2;
using Tuple2 = Vec2;

struct Vec3 {
  double x = 0;
  double y = 0;
  double z = 0;

  Vec3() = default;
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
};

using Point3 = Vec3;
using u8 = std::uint8_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using usize = std::size_t;

struct Color {
  u8 red = 0;
  u8 green = 0;
  u8 blue = 0;

  Color() = default;
  Color(u8 red, u8 green, u8 blue) : red(red), green(green), blue(blue) {}
};

struct SceneOptions {
  Point3 origin = {0.0, 0.0, 0.0};
  double viewport_width = 1;
  double viewport_height = 1;
  double dist_to_proj_plane = 1.0;
  usize canvas_width = 1920;   // num pixels
  usize canvas_height = 1080;  // num pixels
  double t_min = 1.0;
  double t_max = std::numeric_limits<double>::infinity();
  Color background_color{0, 0, 0};
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) {
  return Vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vec3& operator+=(Vec3& a, const Vec3& b) {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

inline Vec3 operator-(const Vec3& a) { return Vec3(-a.x, -a.y, -a.z); }

inline Vec3 operator-(const Vec3& a, const Vec3& b) {
  return operator+(a, operator-(b));
}

inline Vec3& operator-=(Vec3& a, const Vec3& b) {
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}

inline Vec3 operator*(double scale, const Vec3& a) {
  return Vec3(scale * a.x, scale * a.y, scale * a.z);
}

inline Vec3 operator*(const Vec3& a, double scale) {
  return operator*(scale, a);
}

double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

class Ray {
 private:
  Point3 origin{};
  Vec3 direction{};

 public:
  Ray() = default;
  Ray(const Point3& origin, const Vec3& direction)
      : origin(origin), direction(direction) {}

  Point3 point_at(double t) const {
    // we generally want t to be greater than d - the distance from camera to
    // the viewport. However, this is not a requirement here.
    return origin + t * direction;
  }

  const Point3& get_origin() const { return origin; }
  const Vec3& get_direction() const { return direction; }
  void set_ray_direction(const Vec3& direction) { this->direction = direction; }
};

// TODO: It may be worth creating a
// base class SceneObject that requires that
// each child implements intersects_with.
// The main issue here would be dealing with
// the return value - should this be a vector
// of all the intersection points between
// say a ray and some arbitrary polygon?

struct Sphere {
  double radius = 0;
  Point3 center{};
  Color color{};

  Sphere() = default;
  Sphere(double radius, const Point3& center, const Color& color)
      : radius(radius), center(center), color(color) {}

  using OptDoublePair = std::pair<std::optional<double>, std::optional<double>>;
  OptDoublePair intersects_with(const Ray& ray) const {
    /*
     * the intersection of a ray O(rigin) + t * D(irection)
     * with a sphere of radius r centered at C is given by solving
     * the follwing equation:
     * t^2 * dot(D, D) + 2t * ( dot(CO, D) ) + dot(CO, CO) - r^2 = 0
     * where CO is the vector from the center of the sphere to the origin
     */

    Vec3 D = ray.get_direction();
    Vec3 CO = ray.get_origin() - center;

    double a = dot(D, D);
    double b = 2 * dot(CO, D);
    double c = dot(CO, CO) - radius * radius;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
      return std::pair(std::nullopt, std::nullopt);
    }

    double t1 = (-b + std::sqrt(discriminant)) / (2 * a);
    double t2 = (-b - std::sqrt(discriminant)) / (2 * a);

    return std::pair(t1, t2);
  }
};

/*
 * TODO: May be better to use partial order or something like that here.
 */
template <typename T>
bool in_range_inclusive(const T& val, const T& a, const T& b) {
  if constexpr (std::is_arithmetic_v<T>) {
    return (val >= a) && (val <= b);
  } else {
    static_assert(std::is_arithmetic_v<T>);
  }
}

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

template <typename T>
class Array2d {
 private:
  const usize num_rows{};
  const usize num_cols{};
  std::vector<T> arr{};

 public:
  Array2d(usize num_rows, usize num_cols)
      : num_rows(num_rows), num_cols(num_cols), arr(num_rows * num_cols, T{}) {}

  explicit Array2d(const T& fill_value)
      : arr(num_rows * num_cols, fill_value) {}

  inline T& operator()(usize row, usize col) {
    return arr[row * num_cols + col];
  }

  inline const T& operator()(usize row, usize col) const {
    return arr[row * num_cols + col];
  }

  usize get_numrows() const { return num_rows; }
  usize get_numcols() const { return num_cols; }

  // all of the iterator methods!
  // we may need some sort of a slice later.
  auto begin() { return arr.begin(); }
  auto end() { return arr.end(); }
  auto cbegin() const { return arr.cbegin(); }
  auto cend() const { return arr.cend(); }
  auto rbegin() { return arr.rbegin(); }
  auto rend() { return arr.rend(); }
  auto crbegin() const { return arr.crbegin(); }
  auto crend() const { return arr.crend(); }
};

class PpmWriter {
 private:
  std::ofstream outfile{};
  usize count_written{};

 public:
  PpmWriter(usize num_rows, std::size_t num_cols,
            const std::string& filename = "./output.ppm")
      : outfile(filename), count_written(0) {
    if (!outfile) {
      fmt::print(stderr, "Failed to open output file: {}\n", filename);
      // TODO: can I redirect the output to cout on failure?
      return;
    }
    outfile << "P3\n";
    // first width, then height
    outfile << num_cols << " " << num_rows << "\n";
    outfile << 255 << "\n";
  }

  PpmWriter(const PpmWriter&) = delete;
  PpmWriter& operator=(const PpmWriter&) = delete;

  ~PpmWriter() { outfile.close(); }

  void write(const Color& value) {
    outfile << +value.red << " " << +value.green << " " << +value.blue << "\n";
    count_written++;
  }

  usize num_written() const { return count_written; }
};

class Canvas {
 private:
  usize canvas_width{};
  usize canvas_height{};
  Array2d<Color> screen{800, 800};

 public:
  Canvas() = default;

  Canvas(usize width, usize height)
      : canvas_width(width), canvas_height(height), screen(height, width) {}

  const auto& get_canvas() const { return screen; }

  usize get_width() const { return canvas_width; }
  usize get_height() const { return canvas_height; }

  /*
   * Writes pixel color at location x,y of the canvas. Note that x belongs
   * to
   * [-canvas_width / 2, canvas_width / 2), while y belongs to
   * [-canvas_height / 2, canvas_height / 2).
   *
   * TODO: Note that right now, this method does not actually do what
   * it says it does - this is because the PpmWriter used here outputs
   * values in sequence so if x and y arrive out of order, they will be
   * written out of order. We can of course fix this by using a buffer.
   *
   * @param x: horizontal canvas pixel in [-canvas_width / 2, canvas_width
   * / 2)
   * @param y: vertical canvas pixel in [-canvas_height / 2, canvas_height
   * / 2)
   */
  void put_pixel(i64 x, i64 y, const Color& color) {
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

  void write_to_file(const std::string& filename = "output.ppm") {
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
};

Vec3 canvas_to_viewport(i64 canvas_x, i64 canvas_y,
                        const SceneOptions& scene_options) {
  /*
   * canvas_x and canvas_y are pixel coordinates in canvas space.
   * The canvas is organized so that the origin (0,0) point is in
   * the middle and not the top left.
   */
  // these implicit conversion are a little bit of an issue.
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
