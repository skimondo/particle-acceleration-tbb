#include <chrono>
#include <tbb/global_control.h>
#include <tbb/tbb.h>

#include "colormap.h"
#include "particle.h"
#include "uqam/tp.h"

#include "experiments.h"
#include "potential.h"
#include "potentialparallel.h"

// Raccourcis utilisé dans lab-01
using ns = std::chrono::nanoseconds;

int main(int argc, char **argv) {

  const int num_repetition = 10;

  auto ncpus = std::thread::hardware_concurrency();
  int resolution = 200;
  int max_iter = 500;
  double dt = 2e-9;
  int substeps = 10;
  int numpart = 25;
  std::string outfmt("results/potential-{:06d}.png");
  bool update_scale = false;
  int repeat = 10;
  ColorMap cmap;
  cmap.load(SOURCE_DIR "/data/colormap_parula.png");
  bool verbose = false;

  std::ostringstream neant;
  std::ofstream log("bench-random.dat");

  std::vector<Particle> particles_serial;
  std::vector<Particle> particles_parallel;

  experiment_crystal(numpart, particles_serial);
  experiment_crystal(numpart, particles_parallel);

  IPotential *serial = new PotentialSerial(resolution, resolution);
  IPotential *parallel = new PotentialParallel(resolution, resolution);

  // SÉRIE

  // SÉRIE
  long long total_elapsed_time_serial = 0;
    for (int i = 0; i < num_repetition; ++i) {
      auto t1 = std::chrono::high_resolution_clock::now();
      serial->run(particles_serial, max_iter, dt, substeps, update_scale, cmap, outfmt, verbose);
      auto t2 = std::chrono::high_resolution_clock::now();

      auto elapsed_time = std::chrono::duration_cast<ns>(t2 - t1).count();
      total_elapsed_time_serial += elapsed_time;
    }

    auto avg_elapsed_time_serial = total_elapsed_time_serial / num_repetition;

  log << "# ncpu temps acceleration" << "\n";

  // PARALLEL
  {
    ncpus = ncpus > 0 ? ncpus : 1;

    for (int num_threads = 1; num_threads <= ncpus; ++num_threads) {
      long long total_elapsed_time_parallel = 0; // Reset total elapsed time for each thread count

              // Set the number of threads globally for TBB
      tbb::global_control gc(tbb::global_control::max_allowed_parallelism, num_threads);

      for (int i = 0; i < num_repetition; ++i) {
        auto t1 = std::chrono::high_resolution_clock::now();

                // Run the parallel function
        parallel->run(particles_parallel, max_iter, dt, substeps, update_scale, cmap, outfmt, verbose);

        auto t2 = std::chrono::high_resolution_clock::now();

        auto elapsed_time = std::chrono::duration_cast<ns>(t2 - t1).count();
        total_elapsed_time_parallel += elapsed_time;
      }

      auto avg_elapsed_time_parallel = total_elapsed_time_parallel / num_repetition; // Calculate average time
      double acceleration = static_cast<double>(avg_elapsed_time_serial) / avg_elapsed_time_parallel;

      log << num_threads << " " << avg_elapsed_time_parallel << " " << acceleration << "\n";
    }
  }
  std::cout << "DONE!" << std::endl;
}
