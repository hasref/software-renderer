
#include "ray.h"

#include "vec.h"
#include "vec_utils.h"

Ray::Ray(const Point3& origin, const Vec3& direction)
    : origin(origin), direction(direction) {}

Point3 Ray::point_at(double t) const {
  // we generally want t to be greater than d - the distance from camera to
  // the viewport. However, this is not a requirement here.
  return origin + t * direction;
}

const Point3& Ray::get_origin() const { return origin; }

const Vec3& Ray::get_direction() const { return direction; }

void Ray::set_ray_direction(const Vec3& direction) {
  this->direction = direction;
}
