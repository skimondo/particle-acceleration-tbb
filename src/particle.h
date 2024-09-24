#pragma once

#include <Eigen/Dense>

using namespace Eigen;

extern const double k;
extern const double eps;

struct Particle {

  // Constructeur: position, vitesse et charge
  Particle(const Vector2d &p_, const Vector2d &v_, double q_)
      : m_x(p_), m_p(0, 0), m_v(v_), m_f(0, 0), m_q(q_) {}

  // Contribution au potentiel de la particule courant à la distance passée en
  // argument.
  inline double potential_at(const Vector2d &loc) const {
    /*
     * le champs du potentiel contient une discontinuité à la position de la
     * charge elle-même On ajouter une petite constante pour éviter une
     * division par zéro. N'est pas physiquement réaliste.
     *
     * la norme norm() est l'équivalent de la distance euclidienne
     */
    double r = (loc - m_x).norm() + 1e-2;
    return k * m_q / r;
  };

  // Calcule la force que l'autre particule exerice sur celle-ci
  inline Vector2d force(const Particle &o) const {
    Vector2d dir = m_x - o.m_x;           // direction
    double r = dir.norm() + eps;          // distance
    double f = k * o.m_q * m_q / (r * r); // loi de coulomb
    return f * dir.normalized();
  };

  // Added by Dominique
  bool operator==(const Particle &other) const {
    return m_x == other.m_x &&
           m_p == other.m_p &&
           m_v == other.m_v &&
           m_f == other.m_f &&
           m_q == other.m_q;
  }

  Vector2d m_x; // position courante
  Vector2d m_p; // position precedente
  Vector2d m_v; // vitesse
  Vector2d m_f; // force
  double m_q;   // charge
};
