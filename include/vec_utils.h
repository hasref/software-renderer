#ifndef VEC_UTILS_H_
#define VEC_UTILS_H_

#include <algorithm>
#include <cassert>
#include <cmath>

#include "vec.h"

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

inline Vec3 operator/(const Vec3& a, double scalar) {
  assert((scalar != 0.0) && "Cannot divide Vec3 by 0");
  return Vec3{a.x / scalar, a.y / scalar, a.z / scalar};
}

inline double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double length(const Vec3& a) { return std::sqrt(dot(a, a)); }

inline Color operator*(const Color& color, double scalar) {
  // maybe it would be better to just use c-style casts here?
  usize red = static_cast<usize>(static_cast<double>(color.red) * scalar);
  usize green = static_cast<usize>(static_cast<double>(color.green) * scalar);
  usize blue = static_cast<usize>(static_cast<double>(color.blue) * scalar);

  Color out_color{};
  // need explicit casting due to how std::clamp is defined
  out_color.red = static_cast<u8>(std::clamp(red, usize{0}, usize{255}));
  out_color.green = static_cast<u8>(std::clamp(green, usize{0}, usize{255}));
  out_color.blue = static_cast<u8>(std::clamp(blue, usize{0}, usize{255}));
  return out_color;
}

inline Color operator*(double scalar, const Color& color) {
  return operator*(color, scalar);
}

inline Color& operator*=(Color& color, double scalar) {
  color = operator*(color, scalar);
  return color;
}

#endif
