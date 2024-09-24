#include "potentialparallel.h"

/*
 * IDENTIFICATION : indiquer le ou les codes permanents des auteurs ici.
 *
 *
 */

void PotentialParallel::compute_field(std::vector<Particle> &charge, double &lo,
                                      double &hi) {
  // À COMPLETER
  lo = std::numeric_limits<double>::max();
  hi = std::numeric_limits<double>::min();
  int n = charge.size();
  for (int i = 0; i < m_height; i++) {
    for (int j = 0; j < m_width; j++) {
      double x = 1.0 * j / m_width;
      double y = 1.0 * i / m_height;

      double v = 0.0;
      for (int k = 0; k < n; k++) {
        v += charge[k].potential_at({x, y});
      }
      m_sol[IDX2(i, j, m_width)] = v;
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
  }
}

void PotentialParallel::move_particles(std::vector<Particle> &charge, double dt,
                                       int substeps) {
  // À COMPLETER
  int n = charge.size();
  double ssdt = dt / substeps;
  for (int ss = 0; ss < substeps; ss++) {

    // Calculer les forces entre les charges
    for (int i = 0; i < n; i++) {
      Particle &c = charge[i];
      Vector2d f = Vector2d::Zero();
      for (int j = 0; j < n; j++) {
        if (i != j) {
          f += c.force(charge[j]);
        }
      }
      c.m_f = f;
    }

    // On déplace ensuite les charges
    for (int i = 0; i < n; i++) {

      Particle &c = charge[i];
      // mettre à jour la vitesse en fonction de l'accélération
      c.m_v += c.m_f * ssdt;

      // savegarder la position courante
      c.m_p = c.m_x;

      // déplacer la particule en fonction de sa vitesse
      c.m_x += c.m_v * ssdt;
    }
  }
}

void PotentialParallel::save_solution(std::ostream &ofs, ColorMap &cmap) {
  // À COMPLETER
  for (int i = 0; i < m_height; i++) {
    for (int j = 0; j < m_width; j++) {
      double v = m_sol[IDX2(i, j, m_width)];
      png::rgb_pixel pix = cmap.get_color(v);
      // inversion axe y pour correspondre au plan cartésien
      m_img.set_pixel(j, m_height - i - 1, pix);
    }
  }
  m_img.write_stream(ofs);
}
