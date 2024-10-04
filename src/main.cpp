#include <uqam/tp.h>

#include <Eigen/Dense>
#include <filesystem>
#include <format>
#include <iostream>
#include <vector>

#include "colormap.h"
#include "experiments.h"
#include "optparser.hpp"
#include "particle.h"
#include "potential.h"
#include "potentialparallel.h"

using namespace Eigen;

// Calcule l'indice linéaire d'un élément
// (i,j) pour un tableau 2D avec une largeur w
#define IDX2(i, j, w) (i * w + j)

int main(int argc, char** argv) {
  // Charger la carte des couleurs
  std::string colormap_name(SOURCE_DIR "/data/colormap_parula.png");
  bool verbose = false;
  int max_iter = 1000;
  double dt = 1e-9;
  int resol = 512;
  int experiment = 0;
  std::string outfmt("results/potential-{:06d}.png");
  bool update_scale = false;
  int parallel = 0;
  int numpart = 10;
  int substeps = 10;

  OptionsParser args(argc, argv);
  args.AddOption(&colormap_name, "-c", "--color-map", "color map file");
  args.AddOption(&max_iter, "-i", "--iter", "number of time step");
  args.AddOption(&dt, "-dt", "--timestep", "timestep");
  args.AddOption(&verbose, "-v", "--verbose", "-q", "--quiet", "timestep");
  args.AddOption(&resol, "-r", "--resolution", "image resolution");
  args.AddOption(&experiment, "-e", "--experiment", "experiment id");
  args.AddOption(&outfmt, "-o", "--output", "output file template");
  args.AddOption(&update_scale, "-us", "--update-scale", "-nus", "--no-update-scale",
                 "update color scale at each timestep");
  args.AddOption(&parallel, "-p", "--parallel", "parallel engine (0: serial, 1: tbb)");
  args.AddOption(&numpart, "-n", "--num-particles", "number of particles to create");
  args.AddOption(&substeps, "-s", "--substeps", "number of substeps per timestep");

  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(std::cout);
    return 1;
  }
  // args.PrintOptions(std::cout);

  std::filesystem::path path(outfmt);
  std::filesystem::path basedir = path.parent_path();
  if (!basedir.empty() && !std::filesystem::exists(basedir)) {
    std::filesystem::create_directories(basedir);
  }

  // Charger la carte de couleurs
  ColorMap cmap;
  cmap.load(colormap_name);

  // Générer des charges
  std::vector<Particle> particles;
  if (experiment == 0) {
    experiment_random(numpart, particles);
  } else if (experiment == 1) {
    experiment_crystal(numpart, particles);
  } else if (experiment == 2) {
    experiment_collision(numpart, particles);
  } else if (experiment == 3) {
    experiment_basic(particles);
  }

  // Afficher les particules
  if (verbose) {
    std::cout << "*** PARTICULES DEBUT ***" << std::endl;
    std::cout << std::fixed << std::setprecision(3) << std::setw(6);
    for (int i = 0; i < particles.size(); i++) {
      const Particle& c = particles[i];
      std::cout << "particule " << i << " " << c.m_x.transpose() << " " << c.m_v.transpose() << " " << c.m_q << "\n";
    }
    std::cout << "*** PARTICULES FIN ***" << std::endl;
  }

  // Calculer le potentiel électrique sur la grille
  IPotential* simulator;
  if (parallel == 0) {
    simulator = new PotentialSerial(resol, resol);
  } else {
    simulator = new PotentialParallel(resol, resol);
  }

  simulator->run(particles, max_iter, dt, substeps, update_scale, cmap, outfmt, verbose);

  delete simulator;
  // std::cout << "Fin normale du programme" << std::endl;
  return 0;
}
