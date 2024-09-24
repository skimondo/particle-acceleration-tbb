#include <chrono>
#include <tbb/global_control.h>

#include "colormap.h"
#include "particle.h"
#include "uqam/tp.h"

int main(int argc, char **argv) {

  // Paramètres (à ajuster pour votre situation)
  int resolution = 2000;
  std::vector<Particle> particles;
  double dt = 1e-9;
  int substeps = 10;
  int numpart = 100;

  int repeat = 10;
  ColorMap cmap;
  cmap.load(SOURCE_DIR "/data/colormap_parula.png");

  std::ostringstream neant;

  // SÉRIE
  {
    // À COMPLETER

  }

  // PARALLEL
  {
    // À COMPLETER
  }

  // RAPPORT
  {
    // À COMPLETER
  }

  std::cout << "DONE!" << std::endl;
}
