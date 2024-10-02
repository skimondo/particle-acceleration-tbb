#include "potentialparallel.h"

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range2d.h>

/*
 * Dominique Elias
 * ELID14019800
 * 01/10/2024
 *
 */

/*
 * since parallel reduce can only return one variable,
 * this struct is essential to combine the two variables
 * that are dependent inside the parallel reduce
 */
struct LoHi {
  double lo;
  double hi;

  // Initial value for parallel reduce
  LoHi() :
      lo(std::numeric_limits<double>::max()),
      hi(std::numeric_limits<double>::min()) {}

  // used inside combine
  LoHi(double low, double high) : lo(low), hi(high) {}

  // reduce function
  LoHi combine(const LoHi& other) const {
    return LoHi(std::min(lo, other.lo), std::max(hi, other.hi));
  }
};

void PotentialParallel::compute_field(std::vector<Particle> &charge, double &lo,
    double &hi) {
  int n = charge.size();
  auto lohi = tbb::parallel_reduce(
    tbb::blocked_range<int>(0, m_height),
        LoHi(),
        [&](tbb::blocked_range<int> r, LoHi local_lohi) {
      for(int i = r.begin(); i < r.end(); i++) {
        for (int j = 0; j < m_width; j++) {
          double x = 1.0 * j / m_width;
          double y = 1.0 * i / m_height;

          double v = 0.0;
          // since n is relativly low, parallalizing here is not worth it
          for (int k = 0; k < n; k++) {
            v += charge[k].potential_at({x, y});
          }
          m_sol[IDX2(i, j, m_width)] = v;
          local_lohi.lo = std::min(local_lohi.lo, v);
          local_lohi.hi = std::max(local_lohi.hi, v);
        }
      }
      return local_lohi;
    }, [](const LoHi& a, const LoHi& b) {
      return a.combine(b);
    }
  );
  lo = lohi.lo;
  hi = lohi.hi;
}

void PotentialParallel::move_particles(std::vector<Particle> &charge, double dt,
    int substeps) {
  int n = charge.size();
  double ssdt = dt / substeps;
  for (int ss = 0; ss < substeps; ss++) {
    tbb::parallel_for(tbb::blocked_range<int>(0, n),
        [&](const tbb::blocked_range<int>& r) {
      for (int i = r.begin(); i < r.end(); ++i) {
        Particle &c = charge[i];
        Vector2d f = Vector2d::Zero();
        for (int j = 0; j < n; j++) {
          if (i != j) {
            f += c.force(charge[j]);
          }
        }
        c.m_f = f;
      }
    });
    tbb::parallel_for(tbb::blocked_range<int>(0, n),
        [&](const tbb::blocked_range<int>& r) {
      for (int i = r.begin(); i < r.end(); ++i) {
        Particle &c = charge[i];
        c.m_v += c.m_f * ssdt;
        c.m_p = c.m_x;
        c.m_x += c.m_v * ssdt;
      }
    });
  }
}

void PotentialParallel::save_solution(std::ostream &ofs, ColorMap &cmap) {
  tbb::parallel_for(
      tbb::blocked_range2d<int>(0, m_height, 0, m_width),
      [&](const tbb::blocked_range2d<int>& r) {
        for (int i = r.rows().begin(); i < r.rows().end(); ++i) {
          for (int j = r.cols().begin(); j < r.cols().end(); ++j) {
            double v = m_sol[IDX2(i, j, m_width)];
            png::rgb_pixel pix = cmap.get_color(v);
            m_img.set_pixel(j, m_height - i - 1, pix);
          }
        }
      });
  m_img.write_stream(ofs);
}
