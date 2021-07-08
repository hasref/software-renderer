#ifndef SPHERE_H_
#define SPHERE_H_

#include <cmath>
#include <optional>

#include "ray.h"
#include "vec.h"

struct Sphere {
  double radius = 0;
  Point3 center{};
  Color color{};

  Sphere() = default;
  Sphere(double radius, const Point3& center, const Color& color);

  using OptDoublePair = std::pair<std::optional<double>, std::optional<double>>;
  OptDoublePair intersects_with(const Ray& ray) const;
};

#endif
