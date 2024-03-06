#pragma once

#include "common.hpp"

class Segment {
 private:
 public:
  Eigen::Vector2d p1, p2;

  Segment() {
    p1.Zero();
    p2.Zero();
  };
  Segment(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2) {
    this->p1 = p1;
    this->p2 = p2;
  };

  double getDis(const Eigen::Vector2d& p) const {
    Eigen::Vector2d base  = p2 - p1;
    double          cross = base.x() * (p - p1).y() - base.y() * (p - p1).x();
    return abs(cross) / base.norm();
  }

  Eigen::Vector2d project(const Eigen::Vector2d& p) const {
    Eigen::Vector2d base = this->p2 - this->p1;
    double          r    = (p - this->p1).dot(base) / base.squaredNorm();
    return this->p1 + base * r;
  }

  Eigen::Vector2d reflect(const Eigen::Vector2d& p) const {
    return p + (this->project(p) - p) * 2.0;
  }
};
typedef Segment Line;