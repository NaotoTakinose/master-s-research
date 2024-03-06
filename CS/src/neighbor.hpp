#pragma once

#include "bucket.hpp"
#include "common.hpp"
#include "domain.hpp"
#include "math.hpp"
#include "particle.hpp"

namespace neighbor {
void        setNeighbors(std::vector<Particle>& particles, Domain& domain, Bucket& bucket) {
#pragma omp parallel for
  for (auto& pi : particles) {
    if (pi.type != ParticleType::Fluid) continue;

    pi.neighbors.clear();

    int ix = (int)((pi.position.x() - domain.xMin) / bucket.length) + 1;
    int iy = (int)((pi.position.y() - domain.yMin) / bucket.length) + 1;

    for (int jx = ix - 1; jx <= ix + 1; jx++) {
      for (int jy = iy - 1; jy <= iy + 1; jy++) {
        int jb = jx + jy * bucket.numX;
        int j  = bucket.first[jb];

        while (j != -1) {
          double dis = (particles[j].position - pi.position).norm();

          if (j != pi.id && dis < (pi.re + particles[j].re) / 2.0) {
            pi.neighbors.emplace_back(j);
          }
          j = bucket.next[j];
        }
      }
    }
  }
}
}  // namespace neighbor