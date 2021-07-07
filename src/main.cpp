#include <fmt/core.h>

#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

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

struct Color {
  double r = 0;
  double g = 0;
  double b = 0;

  Color() = default;
  Color(double r, double g, double b) : r(r), g(g), b(b) {}
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
};

// TODO: It may be worth creating a
// base class SceneObject that requires that
// each child implements intersects_with.
// The main issue here would be dealing with
// the return value - should this be a vector
// of all the intersection points between
// say a ray and some arbitrary polygon?
using OptDoublePair = std::pair<std::optional<double>, std::optional<double>>;

struct Sphere {
  double radius = 0;
  Point3 center{};

  Sphere() = default;
  Sphere(double radius, const Point3& center)
      : radius(radius), center(center) {}

  OptDoublePair intersects_with(const Ray& ray) {
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

int main() {
  fmt::print("Hello, World\n");
  return 0;
}
