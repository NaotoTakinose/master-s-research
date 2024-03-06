#pragma once
#include <iomanip>
#include <sstream>
#include <vector>

#include "particle.hpp"

namespace writeData {

void writeProf(
    const std::stringstream     &ss,
    const double                &time,
    const std::vector<Particle> &particles) {
  std::ofstream ofs(ss.str());
  if (ofs.fail()) {
    std::cerr << "cannot write " << ss.str() << std::endl;
  }

  ofs << time << std::endl;
  ofs << particles.size() << std::endl;
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.id << " ";
    ofs << static_cast<int>(p.type) << " ";
    ofs << p.position.x() << " " << p.position.y() << " ";
    ofs << p.height << " ";
    ofs << p.velocity.x() << " " << p.velocity.y() << " ";
    ofs << p.volume << " ";                                      // splitのせいで各粒子の体積が違うのでリスタート時に必要
    ofs << p.initialHeight << " " << p.initialDensity2D << " ";  // volumeSummationのリスタート時に必要
    ofs << std::endl;
  }
}  // namespace void writeProf()

void dataArrayBegin(
    std::ofstream     &ofs,
    const std::string &numberOfComponents,
    const std::string &type,
    const std::string &name) {
  ofs << "<DataArray NumberOfComponents='" << numberOfComponents << "' type='"
      << type << "' Name='" << name << "' format='ascii'>" << std::endl;
}

void dataArrayEnd(std::ofstream &ofs) {
  ofs << "</DataArray>" << std::endl;
}

void writeVtu(
    const std::stringstream    &ss,
    const double               &time,
    const std::vector<Particle> particles) {
  std::ofstream ofs(ss.str());
  if (ofs.fail()) {
    std::cerr << "cannot write " << ss.str() << std::endl;
    std::exit(-1);
  }

  // --------------
  // --- Header ---
  // --------------
  int numParticles = 0;
  for (auto &pi : particles) {
    if (pi.type == ParticleType::Ghost) continue;
    numParticles++;
  }
  ofs << "<?xml version='1.0' encoding='UTF-8'?>" << endl;
  ofs << "<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>" << endl;
  ofs << "<UnstructuredGrid>" << endl;
  ofs << "<Piece NumberOfCells='" << numParticles << "' ";
  ofs << "NumberOfPoints='" << numParticles << "'>" << endl;

  /// ------------------
  /// ----- Points -----
  /// ------------------
  ofs << "<Points>" << endl;
  ofs << "<DataArray NumberOfComponents='3' type='Float64' Name='position' format='ascii'>"
      << endl;
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.position.x() << " ";
    ofs << p.position.y() << " ";
    if (p.type == ParticleType::Fluid)
      ofs << p.height + p.bottomHeight << std::endl;
    else
      ofs << 0.0 << std::endl;
  }
  ofs << "</DataArray>" << std::endl;
  ofs << "</Points>" << std::endl;

  // ---------------------
  // ----- PointData -----
  // ---------------------
  ofs << "<PointData>" << std::endl;

  dataArrayBegin(ofs, "1", "Int32", "Particle Type");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << static_cast<int>(p.type) << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "3", "Float64", "Velocity");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.velocity.x() << " ";
    ofs << p.velocity.y() << " ";
    ofs << 0.0 << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "Volume");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.volume << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "Height");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.height << endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "Area");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.area() << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "Particle Size");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << sqrt(p.area()) << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "Effective Radius");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.re << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Int32", "Neighboring Particles Number");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.neighbors.size() << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "3", "Float64", "grad(Flux)");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.gradFlux.x() << " ";
    ofs << p.gradFlux.y() << " ";
    ofs << 0.0 << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "3", "Float64", "Artificial Viscosity");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.artificialViscosity.x() << " ";
    ofs << p.artificialViscosity.y() << " ";
    ofs << 0.0 << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "3", "Float64", "grad(Bottom Height)");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.gradBottomHeight.x() << " ";
    ofs << p.gradBottomHeight.y() << " ";
    ofs << 0.0 << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "3", "Float64", "Friction");
  for (auto &p : particles) {
    if (p.type == ParticleType::Ghost) continue;
    ofs << p.friction.x() << " ";
    ofs << p.friction.y() << " ";
    ofs << 0.0 << std::endl;
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "dudx");
  for (auto &pi : particles) {
    if (pi.type == ParticleType::Ghost) continue;
    if (pi.dudx<=1e-100){
      ofs << 0 << endl;
    }
    else{
      ofs << pi.dudx << endl;
    }
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "dvdy");
  for (auto &pi : particles) {
    if (pi.type == ParticleType::Ghost) continue;
    if (pi.dvdy<=1e-100){
      ofs << 0 << endl;
    }
    else{
      ofs << pi.dvdy << endl;
    }
  }
  dataArrayEnd(ofs);

  dataArrayBegin(ofs, "1", "Float64", "div(Velocity)");
  for (auto &pi : particles) {
    if (pi.type == ParticleType::Ghost) continue;
    if (pi.divVelocity<=1e-100){
      ofs << 0 << endl;
    }
    else{
      ofs << pi.divVelocity << endl;
    }
  }
  dataArrayEnd(ofs);

  ofs << "</PointData>" << std::endl;

  // -----------------
  // ----- Cells -----
  // -----------------
  ofs << "<Cells>" << std::endl;
  ofs << "<DataArray type='Int32' Name='connectivity' format='ascii'>"
      << std::endl;
  for (int i = 0; i < particles.size(); i++) {
    ofs << i << std::endl;
  }
  ofs << "</DataArray>" << std::endl;
  ofs << "<DataArray type='Int32' Name='offsets' format='ascii'>" << std::endl;
  for (int i = 0; i < particles.size(); i++) {
    ofs << i + 1 << std::endl;
  }
  ofs << "</DataArray>" << std::endl;
  ofs << "<DataArray type='UInt8' Name='types' format='ascii'>" << std::endl;
  for (int i = 0; i < particles.size(); i++) {
    ofs << "1" << std::endl;
  }
  ofs << "</DataArray>" << std::endl;
  ofs << "</Cells>" << std::endl;

  // ------------------
  // ----- Footer -----
  // ------------------
  ofs << "</Piece>" << std::endl;
  ofs << "</UnstructuredGrid>" << std::endl;
  ofs << "</VTKFile>" << std::endl;
}

}  // namespace writeData
