#ifndef VEC_UTILS_H_
#define VEC_UTILS_H_

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

inline double dot(const Vec3& a, const Vec3& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

#endif
