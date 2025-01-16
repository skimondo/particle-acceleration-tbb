#pragma once

#include "particle.h"

void experiment_basic(std::vector<Particle>& particles);
void experiment_random(int n, std::vector<Particle>& particles);
void experiment_crystal(int n, std::vector<Particle>& particles);
void experiment_collision(int n, std::vector<Particle>& particles);
