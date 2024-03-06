#include <omp.h>

#include <fstream>
#include <iostream>
#include <vector>

#include "bucket.hpp"
#include "common.hpp"
#include "domain.hpp"
#include "neighbor.hpp"
#include "output.hpp"
#include "particle.hpp"
#include "system.hpp"

// 計算条件
constexpr double l0                  = 0.025;
constexpr double re0                 = 3.1 * l0;
double           dt                  = 0.001;  // 変更可能にした
constexpr bool   changeDtFlag        = false;   // 最終結果を出すときはきれいな時間で出力するために時間幅を変えない
constexpr double finishTime          = 5;
constexpr double outputPeriod        = 0.01;
constexpr double gravity             = 9.8;
constexpr bool   frictionFlag        = false;
constexpr double roughnessCoeff      = 0.01;
constexpr bool   volumeSummationFlag = true;
constexpr bool   restartFlag         = false;  // trueだとvolumeSummation時にinitialDensity2Dをprofから読み取る
constexpr double minCourant          = 0.3;
constexpr double maxCourant          = 0.4;

// 粒子分裂
constexpr bool   splitFlag       = true;
constexpr double refHeight       = 0.25;
constexpr double heightThreshold = 0.1;
constexpr double alphaSplit      = 0.8;

// 人工粘性
constexpr bool   viscosityFlag  = false;
constexpr double alphaViscosity = 0.2;

// 保存形
constexpr bool conservativeFluxFlag = true;

// デバッグ
constexpr bool debugModeFlag = false;
// constexpr bool debugModeFlag = true;

enum class NablaFlag {
  velocity,
  velocityX,
  velocityY,
  flux,
  artificialViscosity,
  bottomHeight,
  height,
};

std::vector<Particle> particles;
Bucket                bucket;
Domain                domain;
System                sys(finishTime, outputPeriod, debugModeFlag);

void   input();
void   calInitialDensity2D();
double calCourant();

void            calEffectiveRadius();
void            calHeight();
void            calVelocity();
Eigen::Vector2d calFriction(const Eigen::Vector2d& velocity, const double& height);
void            update();
double          calBottomHeight(const Eigen::Vector2d& position);

void            split();
Eigen::Vector2d calSplitDirection(const Particle& p);

void collision();

double          div(const Particle& pi, const NablaFlag& flag);
Eigen::Vector2d grad(const Particle& pi, const NablaFlag& flag);
Eigen::Matrix2d Jacobian(const Particle& pi, const NablaFlag& flag);
double          weight(const double& dis, const double& re);
double          weightSPH(const double& dis, const double& re);

int main(int argc, char** argv) {
  omp_set_num_threads(4);

  sys.startSimulation();

  input();
  bucket.generate(particles.size());

  if (volumeSummationFlag && !restartFlag) {
    calInitialDensity2D();
  }

  while (sys.time <= sys.finishTime + 0.1 * dt) {
    sys.timestepStartTime = std::clock();

    calEffectiveRadius();
    bucket.resize(particles, domain.xLength, domain.yLength);

    bucket.storeParticles(particles, domain);
    neighbor::setNeighbors(particles, domain, bucket);
    calHeight();
    calVelocity();

    sys.dataOutput(particles);  // 微分情報を出力するためupdateの直前で呼ぶ

    update();

    bucket.storeParticles(particles, domain);
    neighbor::setNeighbors(particles, domain, bucket);
    collision();

    if (splitFlag) {
      bucket.storeParticles(particles, domain);
      neighbor::setNeighbors(particles, domain, bucket);
      split();
    }

    double next_dt         = calCourant();
    sys.timestepFinishTime = std::clock();
    sys.logOutput(dt);

    sys.time += dt;
    sys.timestep++;

    if (changeDtFlag) dt = next_dt;
  }

  sys.endSimulation();

  return 0;
}

void input() {
  std::stringstream ss;
  std::ifstream     ifs;

  ss << "../input/input.prof";
  ifs.open(ss.str());
  if (ifs.fail()) {
    std::cerr << "cannnot read" << ss.str() << std::endl;
    std::exit(-1);
  }

  int particleSize;
  ifs >> sys.time;
  sys.initialTime = sys.time;
  ifs >> particleSize;
  int numParticle = 0;
  rep(i, 0, particleSize) {
    int          id, type_int;
    ParticleType type;
    Vector2d     position, velocity;
    double       height, volume;
    double       initialHeight, initialDensity2D;  // volume-summationのrestart時のみ使う

    ifs >> id;  // 使わない
    ifs >> type_int;
    ifs >> position.x() >> position.y() >> height;
    ifs >> velocity.x() >> velocity.y();
    ifs >> volume;
    ifs >> initialHeight >> initialDensity2D;

    type = static_cast<ParticleType>(type_int);
    if (type != ParticleType::Ghost) {
      // 直後でreは計算するので適当な値を入れておく
      particles.emplace_back(numParticle, type, position, height, velocity, volume, re0, calBottomHeight(position));
      if (volumeSummationFlag && restartFlag) {
        particles.back().initialHeight    = initialHeight;
        particles.back().initialDensity2D = initialDensity2D;
      }
      numParticle++;
    }
  }
  ifs.close();
  ifs.clear();

  ss.str("../input/input.domain");
  ifs.open(ss.str());
  if (ifs.fail()) {
    std::cerr << "cannnot read" << ss.str() << std::endl;
    std::exit(-1);
  }

  std::string dumstr;
  double      xMin, xMax, yMin, yMax;
  ifs >> dumstr >> xMin;
  ifs >> dumstr >> xMax;
  ifs >> dumstr >> yMin;
  ifs >> dumstr >> yMax;
  domain.generate(xMin, xMax, yMin, yMax);
  ifs.close();
  ifs.clear();
}

// volume-summation formulaの場合のみ
void        calInitialDensity2D() {
#pragma omp parallel for
  for (auto& pi : particles) {
    if (pi.type != ParticleType::Fluid)
      continue;

    double density3D    = 1.0;
    pi.initialDensity2D = 0.0;
    for (int ix = -5; ix <= 5; ix++) {
      for (int iy = -5; iy <= 5; iy++) {
        double dis = sqrt(ix * ix + iy * iy) * l0;
        pi.initialDensity2D += density3D * pi.volume * weightSPH(dis, pi.re);
      }
    }
  }
}

// split後粒子のreを設定するために計算の最初で実施
void        calEffectiveRadius() {
#pragma omp parallel for
  for (auto& p : particles) {
    if (p.type != ParticleType::Fluid) continue;

    p.re = (sqrt(p.area()) / l0) * re0;
    if (p.re < re0) p.re = re0;
  }
}

void        calHeight() {
#pragma omp parallel for
  for (auto& pi : particles) {
    if (pi.type != ParticleType::Fluid) continue;

    // volume-summation formula
    if (volumeSummationFlag) {
      double density3D = 1.0;
      double density2D = 0.0;
      density2D += density3D * pi.volume * weightSPH(0.0, pi.re);
      for (auto& neighbor : pi.neighbors) {
        Particle& pj  = particles[neighbor.id];
        double    dis = (pj.position - pi.position).norm();
        double    re  = (pi.re + pj.re) / 2.0;
        double    volume;
        if (pj.type == ParticleType::Fluid) {
          volume = pj.volume;
        } else if (pj.type == ParticleType::Wall) {
          volume = pi.height * l0 * l0;
        }
        density2D += density3D * volume * weightSPH(dis, re);
      }
      pi.nextHeight = (density2D / pi.initialDensity2D) * pi.initialHeight;
      // 要検討
      if (splitFlag) {
        // if (pi.nextHeight < refHeight * heightThreshold) {
        //   pi.nextHeight = refHeight * heightThreshold;
        // }
        // if (pi.volume / pi.nextHeight > l0 * l0 * 3.0) {
        //   pi.nextHeight = pi.volume / (l0 * l0 * 3.0);
        // }
      }

      // continuity equation
    } else {
      pi.dudx         = grad(pi, NablaFlag::velocityX).x();
      pi.dvdy         = grad(pi, NablaFlag::velocityY).y();
      pi.divVelocity  = pi.dudx + pi.dvdy;
      double nextArea = pi.area() + pi.area() * pi.divVelocity * dt;
      pi.nextHeight   = pi.volume / nextArea;
    }
  }
}

void calVelocity() {
  // #pragma omp parallel for
  //   for (auto& pi : particles) {
  //     if (pi.type != ParticleType::Fluid) continue;
  //
  //     pi.minHeight = pi.height;
  //     for (auto& neighbor : pi.neighbors) {
  //       Particle& pj = particles[neighbor.id];
  //       if (pj.type == ParticleType::Wall) continue;
  //
  //       if (pj.height < pi.minHeight) pi.minHeight = pj.height;
  //     }
  //   }

#pragma omp parallel for
  for (auto& p : particles) {
    if (p.type != ParticleType::Fluid) continue;

    if (conservativeFluxFlag) {
      p.gradFlux = grad(p, NablaFlag::flux);  // m^2/s^2
    } else {
      p.gradFlux = gravity * p.height * grad(p, NablaFlag::height);
    }
    if (viscosityFlag) {
      p.artificialViscosity = grad(p, NablaFlag::artificialViscosity);
    } else {
      p.artificialViscosity = Vector2d::Zero();
    }
    Eigen::Vector2d fluxTerm = p.gradFlux + p.artificialViscosity;

    p.friction          = calFriction(p.velocity, p.height);
    p.gradBottomHeight  = grad(p, NablaFlag::bottomHeight);
    Vector2d sourceTerm = -gravity * p.height * (p.gradBottomHeight + p.friction);

    Vector2d momentum = p.volume * p.velocity;  // m^3 * m/s = m^4/s
    momentum += (-p.area() * fluxTerm + p.area() * sourceTerm) * dt;
    p.nextVelocity = momentum / p.volume;
  }
}

Eigen::Vector2d calFriction(const Eigen::Vector2d& vel, const double& height) {
  if (frictionFlag) {
    double n = roughnessCoeff;
    return n * n * vel * vel.norm() / pow(height, 4.0 / 3.0);

  } else {
    return Eigen::Vector2d::Zero();
  }
}

void        update() {
#pragma omp parallel for
  for (auto& p : particles) {
    if (p.type != ParticleType::Fluid) continue;

    p.position += p.velocity * dt;
    p.velocity     = p.nextVelocity;
    p.height       = p.nextHeight;
    p.bottomHeight = calBottomHeight(p.position);
  }
}

double calBottomHeight(const Vector2d& position) {
  return 0.0;
}

void split() {
  vector<int>      ids_forTargetParticles;
  vector<Vector2d> newPos1, newPos2, newVel1, newVel2;
  vector<double>   newHeight1, newHeight2, newVolume1, newVolume2;
  vector<double>   initialHeight, initialDensity2D;

  for (auto& p : particles) {
    if (p.type != ParticleType::Fluid) continue;
    if (p.area() < 1.5 * l0 * l0) continue;

    Eigen::Vector2d splitDirection = calSplitDirection(p);

    double          alpha = alphaSplit;
    Eigen::Vector2d pos1  = p.position + (alpha / 2.0) * sqrt(p.area()) * splitDirection;
    Eigen::Vector2d pos2  = p.position - (alpha / 2.0) * sqrt(p.area()) * splitDirection;

    Eigen::Matrix2d velJacobian = Jacobian(p, NablaFlag::velocity);
    Eigen::Vector2d vel1        = p.velocity + velJacobian * (pos1 - p.position);
    Eigen::Vector2d vel2        = p.velocity + velJacobian * (pos2 - p.position);

    Eigen::Vector2d gradHeight = grad(p, NablaFlag::height);
    double          height1    = p.height + gradHeight.dot(pos1 - p.position);
    double          height2    = p.height + gradHeight.dot(pos2 - p.position);

    double volume1 = (height1 / (height1 + height2)) * p.volume;
    double volume2 = (height2 / (height1 + height2)) * p.volume;

    if (height1 < refHeight * heightThreshold || height2 < refHeight * heightThreshold)
      continue;

    ids_forTargetParticles.emplace_back(p.id);
    newPos1.emplace_back(pos1);
    newPos2.emplace_back(pos2);
    newVel1.emplace_back(vel1);
    newVel2.emplace_back(vel2);
    newHeight1.emplace_back(height1);
    newHeight2.emplace_back(height2);
    newVolume1.emplace_back(volume1);
    newVolume2.emplace_back(volume2);
    if (volumeSummationFlag) {
      initialHeight.emplace_back(p.initialHeight);
      initialDensity2D.emplace_back(p.initialDensity2D);
    }
  }

  rep(i, 0, ids_forTargetParticles.size()) {
    Particle& p = particles[ids_forTargetParticles[i]];

    // 影響半径は直後に決めるので変える必要なし
    p.position     = newPos1[i];
    p.velocity     = newVel1[i];
    p.height       = newHeight1[i];
    p.volume       = newVolume1[i];
    p.bottomHeight = calBottomHeight(newPos1[i]);
    if (volumeSummationFlag) {
      p.initialHeight    = initialHeight[i];
      p.initialDensity2D = initialDensity2D[i];
    }

    // 影響半径は直後に決めるので適当な値を入れておく
    particles.emplace_back(particles.size(), ParticleType::Fluid, newPos2[i], newHeight2[i], newVel2[i], newVolume2[i], sqrt(p.area()), calBottomHeight(newPos2[i]));
    if (volumeSummationFlag) {
      particles.back().initialHeight    = initialHeight[i];
      particles.back().initialDensity2D = initialDensity2D[i];
    }
  }
}

Eigen::Vector2d calSplitDirection(const Particle& pi) {
  Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();  // covariance matrix（共分散行列）
  for (auto& neighbor : pi.neighbors) {
    Particle& pj  = particles[neighbor.id];
    Vector2d  rij = pj.position - pi.position;
    double    dis = rij.norm();
    double    xij = rij.x();
    double    yij = rij.y();
    double    re  = (pi.re + pj.re) / 2.0;

    Eigen::Matrix2d mat;
    mat << xij * xij, xij * yij, yij * xij, yij * yij;
    cov += mat * weight(dis, re) / pow(dis, 2);
  }

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigenSolver(cov);

  Eigen::Vector2d splitDirection;
  double          lambda_0 = eigenSolver.eigenvalues()(0);
  double          lambda_1 = eigenSolver.eigenvalues()(1);
  if (1.5 * abs(lambda_0) < abs(lambda_1)) {
    splitDirection = eigenSolver.eigenvectors().col(0);

  } else if (1.5 * abs(lambda_1) < abs(lambda_0)) {
    splitDirection = eigenSolver.eigenvectors().col(1);

  } else {
    splitDirection = pi.velocity;
    splitDirection.normalize();
  }

  return splitDirection;
}

void collision() {
  double alpha_wall  = 1.0;
  double alpha_fluid = 0.5;
  double e_wall      = 0.0;
  double e_fluid     = 0.2;

#pragma omp parallel for
  for (auto& pi : particles) {
    if (pi.type != ParticleType::Fluid) continue;

    pi.nextVelocity = pi.velocity;

    Eigen::Vector2d xi = pi.position;
    Eigen::Vector2d ui = pi.velocity;
    double          li = sqrt(pi.area());

    for (auto& neighbor : pi.neighbors) {
      Particle& pj = particles[neighbor.id];
      Vector2d  xj = pj.position;
      Vector2d  uj = pj.velocity;
      double    lj = sqrt(pj.area());

      double alpha;
      if (pj.type == ParticleType::Fluid)
        alpha = alpha_fluid;
      else
        alpha = alpha_wall;

      double dis          = (xj - xi).norm();
      double collisionVel = (uj - ui).dot(xj - xi) / dis;
      double collisionDis = alpha * (li + lj) / 2.0;

      double   mi      = pi.height;
      double   mInv_i  = 1.0 / mi;
      Vector2d impulse = Vector2d::Zero();
      double   n       = 0.0;
      bool     flag    = false;
      if (dis < collisionDis && collisionVel < 0.0) {
        double mInv_j, e;
        if (pj.type == ParticleType::Fluid) {
          mInv_j = 1.0 / pj.height;
          e      = e_fluid;

        } else if (pj.type == ParticleType::Wall) {
          mInv_j = 0.0;
          e      = e_wall;
        }

        double mu = 1.0 / (mInv_i + mInv_j);
        double re = (pi.re + pj.re) / 2.0;
        impulse -= ((1.0 + e) * mu * collisionVel * (xj - xi) / dis) * weight(dis, re);
        n += weight(dis, re);
        flag = true;
      }

      if (flag) {
        impulse /= n;
        pi.nextVelocity -= impulse / mi;
      }
    }
  }

#pragma omp parallel for
  for (auto& p : particles) {
    if (p.type != ParticleType::Fluid) continue;

    p.position += (p.nextVelocity - p.velocity) * dt;
    p.velocity = p.nextVelocity;
  }
}

double calCourant() {
  double courant = 0.0;
  double next_dt = dt;
  for (auto& p : particles) {
    if (p.type != ParticleType::Fluid) continue;

    double d        = sqrt(p.area());
    double c        = sqrt(gravity * p.height);
    double iCourant = dt * (c + p.velocity.norm()) / d;
    if (iCourant > courant) {
      courant = iCourant;
    }
  }
  sys.courant = courant;

  if (courant > maxCourant)
    next_dt = (maxCourant / courant) * dt;
  if (courant < minCourant)
    next_dt = (minCourant / courant) * dt;

  return next_dt;
}

double div(const Particle& pi, const NablaFlag& flag) {
  double val = 0.0;
  if (flag == NablaFlag::velocity) {
    val = grad(pi, NablaFlag::velocityX).x() + grad(pi, NablaFlag::velocityY).y();

  } else {
    std::cerr << "ERROR: wrong NablaFlag at divergence calculation. NablaFlag = " << static_cast<int>(flag) << std::endl;
    exit(-1);
  }

  return val;
}

Eigen::Matrix2d Jacobian(const Particle& pi, const NablaFlag& flag) {
  Eigen::Matrix2d val;

  if (flag == NablaFlag::velocity) {
    val << grad(pi, NablaFlag::velocityX).x(), grad(pi, NablaFlag::velocityY).y(),
        grad(pi, NablaFlag::velocityY).x(), grad(pi, NablaFlag::velocityY).y();
  } else {
    std::cerr << "ERROR: wrong NablaFlag at Jacobian calculation. NablaFlag = " << static_cast<int>(flag) << std::endl;
    exit(-1);
  }

  return val;
}

Eigen::Vector2d grad(const Particle& pi, const NablaFlag& flag) {
  Eigen::Vector2d val          = Eigen::Vector2d::Zero();
  double          n            = 0.0;
  bool            devisionFlag = false;
  // 粒子数密度で割るかのフラッグ．nの加算を一度も行わなかった場合はゼロ除算を防ぐため割らない．
  // 例えば全近傍粒子が壁面粒子かつbottomHeightの微分を計算するときなど．

  for (auto& neighbor : pi.neighbors) {
    Particle& pj = particles[neighbor.id];
    double    hj;
    if (pj.type == ParticleType::Fluid)
      hj = pj.height;
    else if (pj.type == ParticleType::Wall)
      hj = pi.height;

    double phi_i, phi_j;
    if (flag == NablaFlag::velocityX) {
      phi_i = pi.velocity.x();
      phi_j = pj.velocity.x();

    } else if (flag == NablaFlag::velocityY) {
      phi_i = pi.velocity.y();
      phi_j = pj.velocity.y();

    } else if (flag == NablaFlag::flux) {
      // phi_i = gravity * pi.minHeight * pi.minHeight / 2.0;
      phi_i = -gravity * pi.height * pi.height / 2.0;  // i + j型
      // phi_i = gravity * pi.height * pi.height / 2.0;  // j-i型
      phi_j = gravity * hj * hj / 2.0;

    } else if (flag == NablaFlag::artificialViscosity) {
      double alpha = alphaViscosity;

      phi_i = 0.0;

      Vector2d xi = pi.position;
      Vector2d xj = pj.position;
      Vector2d ui = pi.velocity;
      Vector2d uj = pj.velocity;
      if ((uj - ui).dot(xj - xi) >= 0.0) {
        phi_j = 0.0;

      } else {
        double hi  = pi.height;
        double dis = (xj - xi).norm();
        double re  = (pi.re + pj.re) / 2.0;

        double hij = (hi + hj) / 2.0;
        // double lij = (sqrt(pi.area()) + sqrt(pj.area())) / 2.0;
        double lij = re;
        double cij = (sqrt(gravity * hi) + sqrt(gravity * hj)) / 2.0;
        phi_j      = -alpha * (hij * lij * cij * (uj - ui).dot(xj - xi)) / (dis * dis);
      }

    } else if (flag == NablaFlag::bottomHeight) {
      if (pj.type == ParticleType::Wall) continue;  // 壁粒子とは計算しない

      phi_i = pi.bottomHeight;
      phi_j = pj.bottomHeight;

    } else if (flag == NablaFlag::height) {
      phi_i = -pi.height;  // i+j型
      phi_j = hj;

    } else {
      std::cerr << "ERROR: wrong NablaFlag at gradient calculation. NablaFlag = " << static_cast<int>(flag) << std::endl;
      exit(-1);
    }

    Vector2d xi  = pi.position;
    Vector2d xj  = pj.position;
    double   Aj  = pj.area();
    double   dis = (xj - xi).norm();
    double   re  = (pi.re + pj.re) / 2.0;
    n += Aj * weight(dis, re);
    val += (phi_j - phi_i) * (xj - xi) * Aj * weight(dis, re) / pow(dis, 2);
    devisionFlag = true;
  }

  if (devisionFlag) {
    val *= 2.0 / n;
  }

  return val;
}

double weight(const double& dis, const double& re) {
  double w = 0.0;
  if (dis <= re)
    w = pow(dis / re - 1.0, 2);

  return w;
}

double weightSPH(const double& dis, const double& re) {
  double w = 0.0;
  if (dis <= re)
    w = (5.0 / (4.0 * 3.141592 * re * re)) * pow(2.0 * (1 - dis / re), 3);

  return w;
}
