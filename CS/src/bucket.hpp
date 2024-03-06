#pragma once

#include <vector>

#include "common.hpp"
#include "domain.hpp"
#include "particle.hpp"

class Bucket {
 private:
 public:
  int              num, numX, numY;
  double           length;
  std::vector<int> next, first, last;

  void generate(const int& particleNum) {
    next.resize(particleNum);
  }

  void resize(
      const std::vector<Particle>& particles,
      const double                 domainLengthX,
      const double                 domainLengthY) {
    length = 0.0;
    for (auto& p : particles) {
      if (p.type != ParticleType::Fluid) continue;

      if (length < 1.5 * p.re) {
        length = 1.5 * p.re;
      }
    }

    numX = (int)(domainLengthX / length) + 3;
    numY = (int)(domainLengthY / length) + 3;
    num  = numX * numY + 1;

    first.resize(num);
    last.resize(num);
    next.resize(particles.size());
  }

  void storeParticles(
      const vector<Particle>& particles,
      const Domain&           domain) {
    rep(i, 0, num) {
      first[i] = -1;
      last[i]  = -1;
    }
    rep(i, 0, particles.size()) {
      next[i] = -1;
    }

    for (auto& p : particles) {
      if (p.type == ParticleType::Ghost) continue;

      int ix = (int)((p.position.x() - domain.xMin) / length) + 1;
      int iy = (int)((p.position.y() - domain.yMin) / length) + 1;
      int ib = ix + iy * numX;

      if ((p.position.x() < domain.xMin || domain.xMax < p.position.x()) ||
          (p.position.y() < domain.yMin || domain.yMax < p.position.y())) {
        std::cerr << "ERROR: particle " << p.id << " is out of domain." << std::endl;
        std::cerr << "x = " << p.position.x() << " y = " << p.position.y() << std::endl;
        std::exit(-1);
      }

      if (last[ib] == -1)
        first[ib] = p.id;
      else
        next[last[ib]] = p.id;
      last[ib] = p.id;
    }
  }
};
