#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <limits>

#include "vec.h"

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

#endif
