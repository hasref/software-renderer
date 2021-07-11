#ifndef SPHERE_H_
#define SPHERE_H_

#include <cmath>
#include <optional>

#include "intdefines.h"
#include "ray.h"
#include "vec.h"

struct Sphere {
  double radius = 0;
  i64 specular = -1;
  double reflectiveness = 0.0;
  Point3 center{};
  Color color{};

  Sphere() = default;
  Sphere(double radius, const Point3& center, const Color& color);

  /*
   * @param radius
   * @param specular: specular coefficient (-1 means do not apply specular
   * reflection)
   * @param center: center in view coordinate system
   * @param color: surface color
   *
   *
   */
  Sphere(double radius, i64 specular, const Point3& center, const Color& color);
  Sphere(double radius,
         i64 specular,
         double reflectiveness,
         const Point3& center,
         const Color& color);

  using OptDoublePair = std::pair<std::optional<double>, std::optional<double>>;
  OptDoublePair intersects_with(const Ray& ray) const;

  bool is_reflective() const;
};

#endif
