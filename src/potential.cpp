#include <format>

#include "potential.h"

void IPotential::run(std::vector<Particle> &particles, int max_iter, double dt,
                     int substeps, bool update_scale, ColorMap &cmap,
                     std::string outfmt, bool verbose) {

  double time = 0.0;
  int iter = 0;
  double lo;
  double hi;

  // Définir l'échelle de couleurs et sauvegarder la solution initiale
  compute_field(particles, lo, hi);
  cmap.update_scale(lo, hi);
  {
    std::string fname = std::vformat(outfmt, std::make_format_args(iter));
    std::ofstream ofs(fname);
    save_solution(ofs, cmap);
  }

  if (verbose) {
    std::cout << std::scientific;
    std::cout << "lo: " << lo << "\n";
    std::cout << "ho: " << hi << "\n";
  }

  while (iter < max_iter) {
    std::cout << "iter: " << ++iter << " time: " << time << std::endl;

    // Déplacement des charges
    move_particles(particles, dt, substeps);

    // Mettre à jour la solution
    compute_field(particles, lo, hi);

    // Mettre à jour l'échelle de couleurs si demandé
    if (update_scale) {
      cmap.update_scale(lo, hi);
    }

    {
      std::string fname = std::vformat(outfmt, std::make_format_args(iter));
      std::ofstream ofs(fname);
      save_solution(ofs, cmap);
    }

    // Incrément du temps absolu de la solution
    time = time + dt;
  }
}

void PotentialSerial::compute_field(std::vector<Particle> &charge, double &lo,
                                    double &hi) {
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

void PotentialSerial::move_particles(std::vector<Particle> &charge, double dt,
                                     int substeps) {
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

void PotentialSerial::save_solution(std::ostream &ofs, ColorMap &cmap) {
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
