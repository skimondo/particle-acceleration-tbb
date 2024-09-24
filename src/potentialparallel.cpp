#include "potentialparallel.h"

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>

/*
 * Dominique Elias
 * ELID14019800
 * 24/09/2024
 *
 * ATTENTION, TO BE OPTIMIZED!!
 */

void PotentialParallel::compute_field(std::vector<Particle> &charge, double &lo,
                                      double &hi) {
  
  lo = std::numeric_limits<double>::max();
  hi = std::numeric_limits<double>::min();
  int n = charge.size();
  tbb::parallel_for(0, m_height, [&](int i) {
    tbb::parallel_for(0, m_width, [&](int j) {
      double x = 1.0 * j / m_width;
      double y = 1.0 * i / m_height;

      double v = 0.0;
      for (int k = 0; k < n; k++) {
        v += charge[k].potential_at({x, y});
      }
      m_sol[IDX2(i, j, m_width)] = v;
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    });
  });
}

void PotentialParallel::move_particles(std::vector<Particle> &charge, double dt,
                                       int substeps) {
  
  int n = charge.size();
  double ssdt = dt / substeps;
  tbb::parallel_for(0, substeps, [&](int ss){

    // Calculer les forces entre les charges
    tbb::parallel_for(0, n, [&](int i) {
      Particle &c = charge[i];
      Vector2d f = Vector2d::Zero();
      for (int j = 0; j < n; j++) {
        if (i != j) {
          f += c.force(charge[j]);
        }
      }
      c.m_f = f;
    });

    // On déplace ensuite les charges
    tbb::parallel_for(0, n, [&](int i) {

      Particle &c = charge[i];
      // mettre à jour la vitesse en fonction de l'accélération
      c.m_v += c.m_f * ssdt;

      // savegarder la position courante
      c.m_p = c.m_x;

      // déplacer la particule en fonction de sa vitesse
      c.m_x += c.m_v * ssdt;
    });
  });
}

void PotentialParallel::save_solution(std::ostream &ofs, ColorMap &cmap) {
  
  tbb::parallel_for(0, m_height, [&](int i) {
    tbb::parallel_for(0, m_width, [&](int j) {
      double v = m_sol[IDX2(i, j, m_width)];
      png::rgb_pixel pix = cmap.get_color(v);
      // inversion axe y pour correspondre au plan cartésien
      m_img.set_pixel(j, m_height - i - 1, pix);
    });
  });
  m_img.write_stream(ofs);
}
