#pragma once
#include "common.hpp"

enum class ParticleType {
  Ghost,
  Fluid,
  Wall
};

class Neighbor {
 private:
 public:
  int id;

  Neighbor(const int& id) {
    this->id = id;
  }
};

class Particle {
 private:
 public:
  int          id;
  ParticleType type;
  double       re;
  Vector2d     position;
  Vector2d     velocity, nextVelocity;
  double       height, nextHeight, minHeight;
  double       bottomHeight, minSurfaceHeight;
  double       volume;
  Vector2d     gradFlux, artificialViscosity;
  Vector2d     gradBottomHeight, friction;

  // volume-summation formula
  double initialDensity2D;
  double initialHeight;

  // continuity equation
  double divVelocity, dudx, dvdy;

  std::vector<Neighbor> neighbors;

  Particle(
      int             id,
      ParticleType    type,
      Eigen::Vector2d pos,
      double          h,
      Eigen::Vector2d vel,
      double          V,
      double          re,
      double          bottomHeight) {
    this->id            = id;
    this->type          = type;
    this->position      = pos;
    this->height        = h;
    this->initialHeight = h;
    this->velocity      = vel;
    this->volume        = V;
    this->re            = re;
    this->bottomHeight  = bottomHeight;

    this->gradFlux            = Vector2d::Zero();
    this->artificialViscosity = Vector2d::Zero();
  }

  double area() const {
    return volume / height;
  }

  double surfaceHeight() const {
    return height + bottomHeight;
  }
};
