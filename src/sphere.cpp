#include "sphere.h"

#include "ray.h"
#include "vec.h"
#include "vec_utils.h"

Sphere::Sphere(double radius, const Point3& center, const Color& color)
    : radius(radius), center(center), color(color) {}

Sphere::Sphere(double radius,
               i64 specular,
               const Point3& center,
               const Color& color)
    : radius(radius), specular(specular), center(center), color(color) {}

Sphere::Sphere(double radius,
               i64 specular,
               double reflectiveness,  // between 0 and 1
               const Point3& center,
               const Color& color)
    : radius(radius),
      specular(specular),
      reflectiveness(reflectiveness),
      center(center),
      color(color) {}

std::pair<std::optional<double>, std::optional<double>> Sphere::intersects_with(
    const Ray& ray) const {
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

bool Sphere::is_reflective() const { return reflectiveness > 0.001; }
