#pragma once

#include <Eigen/Dense>

#include "colormap.h"
#include "particle.h"

// Macro qui calcule l'indice dans un tableau
// linéaire qui correspond à un indice 2D.
#define IDX2(i, j, w) (i * w + j)

class IPotential {
public:
  // Le destructeur doit être virtuel pour une classe de base
  virtual ~IPotential() = default;

  // Calculer le champs du potentiel électrique
  virtual void compute_field(std::vector<Particle> &charge, double &lo,
                             double &hi) = 0;

  // Déplacer les particules en fonction des forces entre elles
  virtual void move_particles(std::vector<Particle> &charge, double dt,
                              int substeps) = 0;

  // Sauvegarder le champs du potentiel électrique courant dans une image
  virtual void save_solution(std::ostream &ofs, ColorMap &cmap) = 0;

  // Exécuter la simulation au complet
  // Oui, il y a beaucoup d'arguments, ce serait sans doute mieux de faire une
  // structure, mais bon, ça marche ;-)
  void run(std::vector<Particle> &particles, int max_iter, double dt,
           int substeps, bool update_scale, ColorMap &cmap, std::string outfmt,
           bool verbose);
};

class PotentialSerial : public IPotential {
public:
  PotentialSerial(int width_, int height_)
      : m_width(width_), m_height(height_), m_sol(m_width * m_height),
        m_img(m_width, m_height) {}

  void compute_field(std::vector<Particle> &charge, double &lo,
                     double &hi) override;
  void move_particles(std::vector<Particle> &charge, double dt,
                      int substeps) override;
  void save_solution(std::ostream &ofs, ColorMap &cmap) override;

  int m_width;
  int m_height;
  std::vector<double> m_sol;
  png::image<png::rgb_pixel> m_img;
};
