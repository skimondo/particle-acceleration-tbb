#include "experiments.h"

#include <random>

void experiment_basic(std::vector<Particle> &particles) {
  const double delta = 0.01;
  particles.push_back(Particle({0.5 - delta, 0.5}, {0, 0}, 1));
  particles.push_back(Particle({0.5 + delta, 0.5}, {0, 0}, -1));
}

void experiment_random(int n, std::vector<Particle> &particles) {
  std::mt19937 rnd(0);
  double margin = 0.2;
  std::uniform_real_distribution<double> pos(0.0 + margin, 1.0 - margin);
  std::uniform_real_distribution<double> q(-1, 1);
  for (int i = 0; i < n; i++) {
    particles.push_back(Particle({pos(rnd), pos(rnd)}, {0, 0}, q(rnd)));
  }
}

void experiment_crystal(int n, std::vector<Particle> &particles) {
  // regular pattern of particles in the center
  int nrows = std::ceil(std::sqrt(n));
  int ncols = n / nrows;
  double x1 = 0.3;
  double x2 = 0.7;
  double step = (x2 - x1) / (n - 1);
  double off = step / 2;
  for (int i = 0; i < nrows; i++) {
    double y = x1 + i * step;
    double begin = 0.0;
    int cols = ncols;
    double q = (i % 2) ? -1 : 1;
    if (i % 2 == 0) {
      begin = off;
      cols = n - 1;
    }
    for (int j = 0; j < cols; j++) {
      particles.push_back(Particle({x1 + begin + j * step, y}, {0, 0}, q));
    }
  }
}

void experiment_collision(int n, std::vector<Particle> &particles) {
  // colliding particles
  double step = 0.5 / (n - 1);

  // ----- >      < +++++
  for (int i = 0; i < n; i++) {
    Particle p({0.1 - i * step, 0.499}, {1e7, 0}, 10);
    Particle n({0.9 + i * step, 0.501}, {-1e7, 0}, 10);
    particles.push_back(p);
    particles.push_back(n);
  }
}
