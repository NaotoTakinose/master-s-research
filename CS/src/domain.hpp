#pragma once

#include "common.hpp"
#include "math.hpp"

class Domain {
 private:
 public:
  double          xMin, xMax;
  double          yMin, yMax;
  double          xLength, yLength;
  Eigen::Vector2d vertexes[4];  // 左下→右下→右上→左上
  Line            sides[4];     // 下→右→上→左

  void generate(
      const double xMin, const double xMax,
      const double yMin, const double yMax) {
    this->xMin = xMin;
    this->xMax = xMax;
    this->yMin = yMin;
    this->yMax = yMax;

    this->xLength = xMax - xMin;
    this->yLength = yMax - yMin;

    vertexes[0] << xMin, yMin;
    vertexes[1] << xMax, yMin;
    vertexes[2] << xMax, yMax;
    vertexes[3] << xMin, yMax;

    rep(i, 0, 4) {
      sides[i] = Line(vertexes[i], vertexes[(i + 1) % 4]);
    }
  };
};