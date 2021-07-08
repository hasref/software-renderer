#ifndef VEC_H_
#define VEC_H_

#include "intdefines.h"

struct Vec2 {
  double x = 0;
  double y = 0;

  Vec2() = default;
  Vec2(double x, double y) : x(x), y(y) {}
};

struct Vec3 {
  double x = 0;
  double y = 0;
  double z = 0;

  Vec3() = default;
  Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct Color {
  u8 red = 0;
  u8 green = 0;
  u8 blue = 0;

  Color() = default;
  Color(u8 red, u8 green, u8 blue) : red(red), green(green), blue(blue) {}
};

using Point2 = Vec2;
using Tuple2 = Vec2;
using Point3 = Vec3;

#endif
