#include "Particle.h"

#include <iostream>
#include <cstdlib>

#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)

static int signal;
static int *signalTag = &signal;

__attribute__((target(mic)))
  void stepFunctionMic(
      double *xPos, double *yPos, double *zPos, 
      double *xVel, double *yVel, double *zVel,
      double *xAcc, double *yAcc, double *zAcc,
      double *masses,
      unsigned int length, double step, double softening) {

    double xDiff[length];
    double yDiff[length];
    double zDiff[length];
    double bottoms[length];

#pragma omp parallel for shared(xAcc, yAcc, zAcc, xPos, yPos, zPos)
    for (unsigned int i = 0; i < length; ++i) {
      double xDiff[length];
      double yDiff[length];
      double zDiff[length];
#pragma simd
      for (unsigned int j = 0; j < i; ++j) {
        xDiff[j] = xPos[j] - xPos[i];  
        yDiff[j] = yPos[j] - yPos[i];
        zDiff[j] = zPos[j] - zPos[i];
      }
#pragma simd
      for (unsigned int j = i + 1; j < length; ++j) {
        xDiff[j] = xPos[j] - xPos[i];  
        yDiff[j] = yPos[j] - yPos[i];
        zDiff[j] = zPos[j] - zPos[i];
      }

      // bottom is scaling factor we divide by, we don't need to separate it out
      double bottoms[length];
#pragma simd
      for (unsigned int j = 0; j < i; ++j) {
        double bottom = xDiff[j] * xDiff[j] + 
          yDiff[j] * yDiff[j] + 
          zDiff[j] * zDiff[j] + 
          softening;

        bottoms[j] = sqrt(bottom * bottom * bottom); // bottom ^(3/2) TODO: MKL
      }
#pragma simd
      for (unsigned int j = i + 1; j < length; ++j) {
        double bottom = xDiff[j] * xDiff[j] + 
          yDiff[j] * yDiff[j] + 
          zDiff[j] * zDiff[j] + 
          softening;

        bottoms[j] = sqrt(bottom * bottom * bottom); // bottom ^(3/2) TODO: MKL
      }
#pragma simd
      for (unsigned int j = 0; j < i; ++j) {
        xDiff[j] /= bottoms[j];
        yDiff[j] /= bottoms[j];
        zDiff[j] /= bottoms[j];
      }
#pragma simd
      for (unsigned int j = i + 1; j < length; ++j) {
        xDiff[j] /= bottoms[j];
        yDiff[j] /= bottoms[j];
        zDiff[j] /= bottoms[j];
      }

      for (unsigned int j = 0; j < length; ++j) {
        if (j == i) continue;
        xAcc[i] += masses[j] * xDiff[j];
        yAcc[i] += masses[j] * yDiff[j];
        zAcc[i] += masses[j] * zDiff[j];
      }
    }

#pragma omp parallel for shared(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc)
    for (unsigned int i = 0; i < length; ++i) {
      // Integrate with Symplectic euler.
      // Important to update velocity and use updated velocity to update position
      xVel[i] += step * xAcc[i];
      yVel[i] += step * yAcc[i];
      zVel[i] += step * zAcc[i];

      // Zero out acceleration now instead because of access patterns? 
      xAcc[i] = 0;
      yAcc[i] = 0;
      zAcc[i] = 0;

      xPos[i] += step * xVel[i];
      yPos[i] += step * yVel[i];
      zPos[i] += step * zVel[i];
    }

  }


void stepParticles(Particles &particles, double step, double softening) {
  double *xPos, *yPos, *zPos, *xVel, *yVel, *zVel, *xAcc, *yAcc, *zAcc, *masses;
  unsigned int length;

  xPos = &particles.positions.xs[0];
  yPos = &particles.positions.ys[0];
  zPos = &particles.positions.zs[0];

  xVel = &particles.velocities.xs[0];
  yVel = &particles.velocities.ys[0];
  zVel = &particles.velocities.zs[0];

  xAcc = &particles.accelerations.xs[0];
  yAcc = &particles.accelerations.ys[0];
  zAcc = &particles.accelerations.zs[0];

  masses = &particles.mass[0];

  length = particles.length;

#pragma offload_transfer target(mic:0)\
  out(xPos, yPos, zPos : length(length) REUSE) wait(signalTag)

#pragma offload target(mic:0)\
  in(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(0) REUSE)\
  in(length, step, softening)\
  signal(signalTag)
  stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses,
      length, step, softening);

}

void init(Particles &particles, double step, double softening) {
  double *xPos, *yPos, *zPos, *xVel, *yVel, *zVel, *xAcc, *yAcc, *zAcc, *masses;
  unsigned int length;

  xPos = &particles.positions.xs[0];
  yPos = &particles.positions.ys[0];
  zPos = &particles.positions.zs[0];

  xVel = &particles.velocities.xs[0];
  yVel = &particles.velocities.ys[0];
  zVel = &particles.velocities.zs[0];

  xAcc = &particles.accelerations.xs[0];
  yAcc = &particles.accelerations.ys[0];
  zAcc = &particles.accelerations.zs[0];

  masses = &particles.mass[0];

  length = particles.length;

#pragma offload target(mic:0)\
  in(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(length) ALLOC)\
  in(length, step, softening)\
  signal(signalTag)
  stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, 
      masses, length, step, softening);

}
