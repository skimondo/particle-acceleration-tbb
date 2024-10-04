#pragma once

#include <algorithm>
#include <png++/png.hpp>

struct ColorMap {
  ColorMap() {
    set_scale(std::numeric_limits<double>::max(), std::numeric_limits<double>::min());
  }

  void load(const std::string& filename) {
    png::image<png::rgb_pixel> img;
    img.read(filename);
    int w = img.get_width() - 4;
    m_colors.resize(w);
    for (int i = 0; i < w; i++) {
      m_colors[i] = img.get_pixel(i + 2, 2);
    }
  }

  void set_scale(double lo, double hi) {
    m_lo = lo;
    m_hi = hi;
    m_scale_inv = 1.0 / (m_hi - m_lo);
  }

  void update_scale(double lo, double hi) {
    m_lo = std::min(lo, lo);
    m_hi = std::max(hi, hi);
    m_scale_inv = 1.0 / (m_hi - m_lo);
  }

  inline png::rgb_pixel get_color(double val) {
    val = std::clamp(val, m_lo, m_hi);
    // Astuce: comme ici le dénominateur est constant,
    // il est beaucoup plus rapide de multiplier par
    // l'inverse que de diviser

    // FIXME: utiliser une autre échelle que
    // linéaire pour atténuer les valeurs extrêmes (tanh?)
    // double norm = (1 + std::tanh(6.5 * (val - m_lo) * m_scale_inv)) / 2;
    double norm = (val - m_lo) * m_scale_inv;
    int idx = norm * (m_colors.size() - 1);
    return m_colors[idx];
  }

  std::vector<png::rgb_pixel> m_colors;
  double m_lo;
  double m_hi;
  double m_scale_inv;
};
