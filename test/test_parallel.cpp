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

#define DATA_DIR SOURCE_DIR "/test/data"

using namespace Catch;

static const double abstol = 1e-6;

TEST_CASE("PotentialParallel") {

  /*
   * Mettre les variables communes aux tests ici
   */

  SECTION("ComputeField") {
    CHECK(1 == 0);
    // À COMPLETER
  }

  SECTION("MoveParticles") {
    CHECK(1 == 0);
    // À COMPLETER
  }

  SECTION("SaveSolution") {
    // on utilise std::ostringstream pour écrire en mémoire plutôt que dans un
    // fichier, on ne s'embête pas à créer un fichier temporaire et le relire,
    // etc.
    std::ostringstream neant;
    CHECK(1 == 0);

    // À COMPLETER
  }
}
