#ifndef RAY_H_
#define RAY_H_

#include "vec.h"
#include "vec_utils.h"

class Ray {
 private:
  Point3 origin{};
  Vec3 direction{};

 public:
  Ray() = default;
  Ray(const Point3& origin, const Vec3& direction);

  Point3 point_at(double t) const;

  const Point3& get_origin() const;
  const Vec3& get_direction() const;
  void set_ray_direction(const Vec3& direction);
};
#endif
