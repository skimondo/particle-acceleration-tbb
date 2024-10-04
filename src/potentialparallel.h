#pragma once

#include "potential.h"

class PotentialParallel : public IPotential {
public:
  PotentialParallel(int width_, int height_)
      : m_width(width_), m_height(height_), m_sol(m_width * m_height), m_img(m_width, m_height) {
  }

  void compute_field(std::vector<Particle>& charge, double& lo, double& hi) override;
  void move_particles(std::vector<Particle>& charge, double dt, int substeps) override;
  void save_solution(std::ostream& ofs, ColorMap& cmap) override;

  int m_width;
  int m_height;
  std::vector<double> m_sol;
  png::image<png::rgb_pixel> m_img;
};
