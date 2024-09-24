#include <cmath>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_predicate.hpp>
#include <catch2/matchers/catch_matchers_quantifiers.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include <colormap.h>
#include <particle.h>
#include <potential.h>
#include <potentialparallel.h>
#include <uqam/tp.h>

#include "experiments.h"

#define DATA_DIR SOURCE_DIR "/test/data"

using namespace Catch;

static const double abstol = 1e-6;

TEST_CASE("PotentialParallel") {

  /*
   * Mettre les variables communes aux tests ici
   */
  // compute var
  int resol = 512;
  // move var
  int dt = 1e-9;
  int substeps = 10;
  // save var
  std::string colormap_name(SOURCE_DIR "/data/colormap_parula.png");
  ColorMap cmap;
  cmap.load(colormap_name);
  IPotential *serial = new PotentialSerial(resol, resol);
  IPotential *parallel = new PotentialParallel(resol, resol);

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
    CHECK(particles_serial == particles_parallel);
    
  }

  SECTION("MoveParticles") {
    serial->move_particles(particles_serial, dt, substeps);
    parallel->move_particles(particles_parallel, dt, substeps);
    CHECK(particles_serial.size() == particles_parallel.size());
    CHECK(particles_serial == particles_parallel);
  }

  SECTION("SaveSolution") {
    // on utilise std::ostringstream pour écrire en mémoire plutôt que dans un
    // fichier, on ne s'embête pas à créer un fichier temporaire et le relire,
    // etc.
    std::ostringstream oss_serial, oss_parallel;
    serial->save_solution(oss_serial, cmap);
    parallel->save_solution(oss_parallel, cmap);
    CHECK(oss_serial.str() == oss_parallel.str());
  }
}
