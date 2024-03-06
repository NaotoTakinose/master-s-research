#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../src/output.hpp"
#include "../src/particle.hpp"

using namespace Eigen;
using namespace std;

constexpr double XMIN = 0.0;
constexpr double XMAX = 9.0;
constexpr double YMIN = -0.5;
constexpr double YMAX = 0.5;
constexpr double l0   = 0.025;

void check_fluid_range(
    std::vector<double> &x_range,
    std::vector<double> &y_range,
    const double &l0, double &eps);

bool isInside(
    Eigen::Vector2d     &position,
    std::vector<double> &x_range,
    std::vector<double> &y_range,
    double              &eps);

void set_fluid(
    std::vector<double>   &fluid_x_range,
    std::vector<double>   &fluid_y_range,
    const ParticleType    &type,
    const double          &h,
    const double          &l0,
    std::vector<Particle> &particles);

void writeDomain();

int main(int argc, char **argv) {
  std::vector<Particle> particles;

  double         eps = 0.01 * l0;
  vector<double> fluid_x_range(2), fluid_y_range(2);

  // Fluid
  fluid_x_range = {XMIN + l0 / 2.0, 4.65 - l0 / 2.0};
  fluid_y_range = {YMIN + l0 / 2.0, YMAX - l0 / 2.0};
  check_fluid_range(fluid_x_range, fluid_y_range, l0, eps);
  set_fluid(fluid_x_range, fluid_y_range, ParticleType::Fluid, 0.25, l0, particles);

  // 下辺
  fluid_x_range = {XMIN - l0 / 2.0 - 4.0 * l0, XMAX + l0 / 2.0 + 4.0 * l0};
  fluid_y_range = {YMIN - l0 / 2.0 - 4.0 * l0, YMIN - l0 / 2.0};
  check_fluid_range(fluid_x_range, fluid_y_range, l0, eps);
  set_fluid(fluid_x_range, fluid_y_range, ParticleType::Wall, 1.0, l0, particles);

  // // 上辺
  fluid_x_range = {XMIN - l0 / 2.0 - 4.0 * l0, XMAX + l0 / 2.0 + 4.0 * l0};
  fluid_y_range = {YMAX + l0 / 2.0, YMAX + l0 / 2.0 + 4.0 * l0};
  check_fluid_range(fluid_x_range, fluid_y_range, l0, eps);
  set_fluid(fluid_x_range, fluid_y_range, ParticleType::Wall, 1.0, l0, particles);

  // 右辺
  fluid_x_range = {XMAX + l0 / 2.0, XMAX + l0 / 2.0 + 4.0 * l0};
  fluid_y_range = {YMIN + l0 / 2.0, YMAX - l0 / 2.0};
  check_fluid_range(fluid_x_range, fluid_y_range, l0, eps);
  set_fluid(fluid_x_range, fluid_y_range, ParticleType::Wall, 1.0, l0, particles);

  // 左辺
  fluid_x_range = {XMIN - l0 / 2.0 - 4.0 * l0, XMIN - l0 / 2.0};
  fluid_y_range = {YMIN + l0 / 2.0, YMAX - l0 / 2.0};
  check_fluid_range(fluid_x_range, fluid_y_range, l0, eps);
  set_fluid(fluid_x_range, fluid_y_range, ParticleType::Wall, 1.0, l0, particles);

  std::stringstream ss;
  ss.str("input.prof");
  writeData::writeProf(ss, 0.0, particles);
  ss.str("input.vtu");
  writeData::writeVtu(ss, 0.0, particles);
  writeDomain();
}

void check_fluid_range(
    std::vector<double> &x_range,
    std::vector<double> &y_range,
    const double &l0, double &eps) {
  std::sort(x_range.begin(), x_range.end());
  std::sort(y_range.begin(), y_range.end());

  double nx = (x_range.at(1) - x_range.at(0)) / l0;
  if (abs(nx - std::round(nx)) > eps) {
    cerr << "x_range is not divisible by particle length." << endl;
    cerr << "x_range: {" << x_range[0] << ", " << x_range[1] << "}" << endl;
    cerr << "y_range: {" << y_range[0] << ", " << y_range[1] << "}" << endl;
    exit(-1);
  }

  double ny = (y_range.at(1) - y_range.at(0)) / l0;
  if (abs(ny - std::round(ny)) > eps) {
    cerr << "y_range is not divisible by particle length." << endl;
    cerr << "x_range: {" << x_range[0] << ", " << x_range[1] << "}" << endl;
    cerr << "y_range: {" << y_range[0] << ", " << y_range[1] << "}" << endl;
    exit(-1);
  }
}

void set_fluid(
    std::vector<double>   &fluid_x_range,
    std::vector<double>   &fluid_y_range,
    const ParticleType    &type,
    const double          &h,
    const double          &l0,
    std::vector<Particle> &particles) {
  int ix_end = round((fluid_x_range.at(1) - fluid_x_range.at(0)) / l0);
  int iy_end = round((fluid_y_range.at(1) - fluid_y_range.at(0)) / l0);

  for (int ix = 0; ix <= ix_end; ix++) {
    for (int iy = 0; iy <= iy_end; iy++) {
      Eigen::Vector2d pos(fluid_x_range.at(0) + (double)(ix)*l0, fluid_y_range.at(0) + (double)(iy)*l0);
      Eigen::Vector2d vel          = Eigen::Vector2d::Zero();
      double          V            = l0 * l0 * h;
      double          re           = l0;  // 使わないので適当な値を入れておく
      double          bottomHeight = 0.0;

      particles.emplace_back(particles.size(), type, pos, h, vel, V, re, bottomHeight);
    }
  }
}

bool isInside(Eigen::Vector2d &pos, std::vector<double> &x_range,
              std::vector<double> &y_range, double &eps) {
  std::sort(x_range.begin(), x_range.end());
  std::sort(y_range.begin(), y_range.end());

  if (((x_range.at(0) - eps < pos.x()) && (pos.x() < x_range.at(1) + eps)) &&
      ((y_range.at(0) - eps < pos.y()) && (pos.y() < y_range.at(1) + eps)))
    return true;
  else
    return false;
}

void writeDomain() {
  std::stringstream ss;
  ss << "input.domain";

  std::ofstream ofs(ss.str());
  if (ofs.fail()) {
    std::cerr << "cannot write " << ss.str() << std::endl;
    std::exit(-1);
  }

  ofs << "xMin " << XMIN - 5.0 * l0 << std::endl;
  ofs << "xMax " << XMAX + 5.0 * l0 << std::endl;
  ofs << "yMin " << YMIN - 5.0 * l0 << std::endl;
  ofs << "yMax " << YMAX + 5.0 * l0 << std::endl;
}