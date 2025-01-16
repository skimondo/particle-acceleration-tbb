#include <colormap.h>
#include <particle.h>
#include <potential.h>
#include <uqam/tp.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#define DATA_DIR SOURCE_DIR "/test/data"

using namespace Catch;

static const double abstol = 1e-6;
static const double reltol = 1e-2;

TEST_CASE("ChargeForce") {
  Vector2d z0 = Vector2d::Zero();

  SECTION("AttireAngle0") {
    Particle c1({-1, 0}, z0, 1e-6);
    Particle c2({1, 0}, z0, 1e-6);
    Vector2d f1 = c1.force(c2);
    Vector2d f2 = c2.force(c1);
    CHECK_THAT((f1 + f2).norm(), Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
    Vector2d orient(-1, 0);
    CHECK_THAT(std::acos(f1.dot(orient) / (f1.norm() * orient.norm())),
               Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
  }

  SECTION("AttireAngle45") {
    Particle c1({0, 0}, z0, 1e-6);
    Particle c2({1, 1}, z0, 1e-6);
    Vector2d f1 = c1.force(c2);
    Vector2d f2 = c2.force(c1);
    CHECK_THAT((f1 + f2).norm(), Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
    Vector2d orient(-1, -1);
    CHECK_THAT(std::acos(f1.dot(orient) / (f1.norm() * orient.norm())),
               Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
  }

  SECTION("RepulseAngle0") {
    Particle c1({-1, 0}, z0, 1e-6);
    Particle c2({1, 0}, z0, -1e-6);
    Vector2d f1 = c1.force(c2);
    Vector2d f2 = c2.force(c1);
    CHECK_THAT((f1 + f2).norm(), Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
    Vector2d orient(1, 0);
    CHECK_THAT(std::acos(f1.dot(orient) / (f1.norm() * orient.norm())),
               Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
  }

  SECTION("Square") {
    std::vector<Particle> charges = {
        Particle({-1, -1}, z0, 1e-6),
        Particle({1, -1}, z0, 1e-6),
        Particle({1, 1}, z0, 1e-6),
        Particle({-1, 1}, z0, 1e-6),
    };

    for (int i = 0; i < charges.size(); i++) {
      Particle& c = charges[i];
      Vector2d f = Vector2d::Zero();
      for (int j = 0; j < charges.size(); j++) {
        if (i != j) {
          f += c.force(charges[j]);
        }
      }
      CHECK_THAT(std::acos(c.m_x.dot(f) / (c.m_x.norm() * f.norm())),
                 Matchers::WithinRel(0, reltol) || Matchers::WithinAbs(0, abstol));
    }
  }
}
