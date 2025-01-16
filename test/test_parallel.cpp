#include <colormap.h>
#include <particle.h>
#include <potential.h>
#include <potentialparallel.h>
#include <uqam/tp.h>

#include <catch2/catch_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cmath>

#include "experiments.h"

#define DATA_DIR SOURCE_DIR "/test/data"

/*
 * Dominique Elias
 * ELID14019800
 * 01/10/2024
 *
 */

using namespace Catch;

static const double abstol = 1e-6;

TEST_CASE("PotentialParallelBasicExperiment") {
  int resol = 512;
  double dt = 1e-9;
  int substeps = 10;
  std::string colormap_name(SOURCE_DIR "/data/colormap_parula.png");
  ColorMap cmap;
  cmap.load(colormap_name);

  IPotential* serial = new PotentialSerial(resol, resol);
  IPotential* parallel = new PotentialParallel(resol, resol);

  std::vector<Particle> particles_serial;
  experiment_basic(particles_serial);
  double lo_serial, hi_serial;

  std::vector<Particle> particles_parallel;
  experiment_basic(particles_parallel);
  double lo_parallel, hi_parallel;

  SECTION("ComputeField") {
    serial->compute_field(particles_serial, lo_serial, hi_serial);
    parallel->compute_field(particles_parallel, lo_parallel, hi_parallel);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }

    CHECK(std::abs(lo_serial - lo_parallel) < abstol);
    CHECK(std::abs(hi_serial - hi_parallel) < abstol);
  }

  SECTION("MoveParticles") {
    serial->move_particles(particles_serial, dt, substeps);
    parallel->move_particles(particles_parallel, dt, substeps);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }
  }

  SECTION("SaveSolution") {
    std::ostringstream oss_serial, oss_parallel;
    serial->save_solution(oss_serial, cmap);
    parallel->save_solution(oss_parallel, cmap);
    CHECK(oss_serial.str() == oss_parallel.str());
  }
  delete serial;
  delete parallel;
}

TEST_CASE("PotentialParallelCrystalExperiment") {
  int resol = 512;
  double dt = 1e-9;
  int substeps = 10;
  std::string colormap_name(SOURCE_DIR "/data/colormap_parula.png");
  ColorMap cmap;
  cmap.load(colormap_name);
  int numpart = 25;

  IPotential* serial = new PotentialSerial(resol, resol);
  IPotential* parallel = new PotentialParallel(resol, resol);

  std::vector<Particle> particles_serial;
  experiment_crystal(numpart, particles_serial);
  double lo_serial, hi_serial;

  std::vector<Particle> particles_parallel;
  experiment_crystal(numpart, particles_parallel);
  double lo_parallel, hi_parallel;

  SECTION("ComputeField") {
    serial->compute_field(particles_serial, lo_serial, hi_serial);
    parallel->compute_field(particles_parallel, lo_parallel, hi_parallel);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }

    CHECK(std::abs(lo_serial - lo_parallel) < abstol);
    CHECK(std::abs(hi_serial - hi_parallel) < abstol);
  }

  SECTION("MoveParticles") {
    serial->move_particles(particles_serial, dt, substeps);
    parallel->move_particles(particles_parallel, dt, substeps);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }
  }

  SECTION("SaveSolution") {
    std::ostringstream oss_serial, oss_parallel;
    serial->save_solution(oss_serial, cmap);
    parallel->save_solution(oss_parallel, cmap);
    CHECK(oss_serial.str() == oss_parallel.str());
  }

  delete serial;
  delete parallel;
}

TEST_CASE("PotentialParallelCollisionExperiment") {
  int resol = 512;
  double dt = 1e-9;
  int substeps = 10;
  std::string colormap_name(SOURCE_DIR "/data/colormap_parula.png");
  ColorMap cmap;
  cmap.load(colormap_name);
  int numpart = 25;

  IPotential* serial = new PotentialSerial(resol, resol);
  IPotential* parallel = new PotentialParallel(resol, resol);

  std::vector<Particle> particles_serial;
  experiment_collision(numpart, particles_serial);
  double lo_serial, hi_serial;

  std::vector<Particle> particles_parallel;
  experiment_collision(numpart, particles_parallel);
  double lo_parallel, hi_parallel;

  SECTION("ComputeField") {
    serial->compute_field(particles_serial, lo_serial, hi_serial);
    parallel->compute_field(particles_parallel, lo_parallel, hi_parallel);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }

    CHECK(std::abs(lo_serial - lo_parallel) < abstol);
    CHECK(std::abs(hi_serial - hi_parallel) < abstol);
  }

  SECTION("MoveParticles") {
    serial->move_particles(particles_serial, dt, substeps);
    parallel->move_particles(particles_parallel, dt, substeps);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }
  }

  SECTION("SaveSolution") {
    std::ostringstream oss_serial, oss_parallel;
    serial->save_solution(oss_serial, cmap);
    parallel->save_solution(oss_parallel, cmap);
    CHECK(oss_serial.str() == oss_parallel.str());
  }

  delete serial;
  delete parallel;
}

TEST_CASE("PotentialParallelRandomExperiment") {
  int resol = 512;
  double dt = 1e-9;
  int substeps = 10;
  std::string colormap_name(SOURCE_DIR "/data/colormap_parula.png");
  ColorMap cmap;
  cmap.load(colormap_name);
  int numpart = 25;

  IPotential* serial = new PotentialSerial(resol, resol);
  IPotential* parallel = new PotentialParallel(resol, resol);

  std::vector<Particle> particles_serial;
  experiment_random(numpart, particles_serial);
  double lo_serial, hi_serial;

  std::vector<Particle> particles_parallel;
  experiment_random(numpart, particles_parallel);
  double lo_parallel, hi_parallel;

  SECTION("ComputeField") {
    serial->compute_field(particles_serial, lo_serial, hi_serial);
    parallel->compute_field(particles_parallel, lo_parallel, hi_parallel);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }

    CHECK(std::abs(lo_serial - lo_parallel) < abstol);
    CHECK(std::abs(hi_serial - hi_parallel) < abstol);
  }

  SECTION("MoveParticles") {
    serial->move_particles(particles_serial, dt, substeps);
    parallel->move_particles(particles_parallel, dt, substeps);

    CHECK(particles_serial.size() == particles_parallel.size());

    for (size_t i = 0; i < particles_serial.size(); ++i) {
      CHECK((particles_serial[i].m_x - particles_parallel[i].m_x).norm() < abstol);
      CHECK((particles_serial[i].m_p - particles_parallel[i].m_p).norm() < abstol);
      CHECK((particles_serial[i].m_v - particles_parallel[i].m_v).norm() < abstol);
      CHECK((particles_serial[i].m_f - particles_parallel[i].m_f).norm() < abstol);
      CHECK(std::abs(particles_serial[i].m_q - particles_parallel[i].m_q) < abstol);
    }
  }

  SECTION("SaveSolution") {
    std::ostringstream oss_serial, oss_parallel;
    serial->save_solution(oss_serial, cmap);
    parallel->save_solution(oss_parallel, cmap);
    CHECK(oss_serial.str() == oss_parallel.str());
  }

  delete serial;
  delete parallel;
}
