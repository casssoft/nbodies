#ifndef STEP_PARTICLES_H
#define STEP_PARTICLES_H

#include <vector>
#include "Particle.h"

void init(Particles&, double step, double softening);
void stepParticles(Particles&, double, double);

#endif
