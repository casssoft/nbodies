#include "Particle.h"

#include <iostream>
#include <cstdlib>

#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)

#define _SIMD

static int signal;
static int *signalTag = &signal;

inline void calculateAcceleration(
    double *xPos, double *yPos, double *zPos,
    double *xAcc, double *yAcc, double *zAcc,
    double *masses, unsigned int length, double softening,
    unsigned int i, unsigned int j) {

  // r_i_j is the distance vector
  //Vector3d r_i_j = particles[j]->getPosition() - particles[i]->getPosition();
  double r_i_j_x = xPos[j] - xPos[i];
  double r_i_j_y = yPos[j] - yPos[i];
  double r_i_j_z = zPos[j] - zPos[i];

  // bottom is scaling factor we divide by, we don't need to separate it out
  double bottom = r_i_j_x * r_i_j_x + r_i_j_y * r_i_j_y + r_i_j_z * r_i_j_z + softening;

  bottom = sqrt(bottom * bottom * bottom); // bottom ^(3/2) TODO: MKL

  //Vector3d f_i_j = r_i_j/ bottom;
  // Resuse r_i_j as f_i_j cause fuck it's verbose otherwise
  r_i_j_x /= bottom;
  r_i_j_y /= bottom;
  r_i_j_z /= bottom;

  // distvector = j pos - i pos
  // particles[i].acceleration accumlator = (m of j/ (dist^2 + e2)) * distvector

  // so f_i_j is the shared part of the calculation between the pair
  // which = distvector/(dist^2 + e2) but I multiply f_i_j by negative 1 to
  // reverse the direction so that it works for particle i too

  // Notice we are just adding to the accelerator
  //particles[i]->setAcceleration(particles[j]->getMass() * f_i_j + particles[i]->getAcceleration());
  //particles[j]->setAcceleration(particles[i]->getMass() * -1 * f_i_j + particles[j]->getAcceleration());
  xAcc[i] += masses[j] * r_i_j_x;
  yAcc[i] += masses[j] * r_i_j_y;
  zAcc[i] += masses[j] * r_i_j_z;

}

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
#if 0
      for (unsigned int j = 0; j < length; ++j) {
        if (j == i) { continue; }
        // r_i_j is the distance vector
        //Vector3d r_i_j = particles[j]->getPosition() - particles[i]->getPosition();
        xDiff[j] = xPos[j] - xPos[i];
        yDiff[j] = yPos[j] - yPos[i];
        zDiff[j] = zPos[j] - zPos[i];

        // bottom is scaling factor we divide by, we don't need to separate it out
        bottoms[j] = xDiff[j] * xDiff[j] + yDiff[j] * yDiff[j] + zDiff[j] * zDiff[j] + softening;

        bottoms[j] = sqrt(bottoms[j] * bottoms[j] * bottoms[j]); // bottom ^(3/2) TODO: MKL

        //Vector3d f_i_j = r_i_j/ bottom;
        // Resuse r_i_j as f_i_j cause fuck it's verbose otherwise
        xDiff[j] /= bottoms[j];
        yDiff[j] /= bottoms[j];
        zDiff[j] /= bottoms[j];

        // distvector = j pos - i pos
        // particles[i].acceleration accumlator = (m of j/ (dist^2 + e2)) * distvector

        // so f_i_j is the shared part of the calculation between the pair
        // which = distvector/(dist^2 + e2) but I multiply f_i_j by negative 1 to
        // reverse the direction so that it works for particle i too

        // Notice we are just adding to the accelerator
        //particles[i]->setAcceleration(particles[j]->getMass() * f_i_j + particles[i]->getAcceleration());
        //particles[j]->setAcceleration(particles[i]->getMass() * -1 * f_i_j + particles[j]->getAcceleration());
        xAcc[i] += masses[j] * xDiff[j];
        yAcc[i] += masses[j] * yDiff[j];
        zAcc[i] += masses[j] * zDiff[j];

      }
#endif

#if 1
      double xDiff[length];
      double yDiff[length];
      double zDiff[length];
#ifdef _SIMD
#pragma simd
#endif
      for (unsigned int j = 0; j < i; ++j) {
        xDiff[j] = xPos[j] - xPos[i];  
        yDiff[j] = yPos[j] - yPos[i];
        zDiff[j] = zPos[j] - zPos[i];
      }
#ifdef _SIMD
#pragma simd
#endif
      for (unsigned int j = i + 1; j < length; ++j) {
        xDiff[j] = xPos[j] - xPos[i];  
        yDiff[j] = yPos[j] - yPos[i];
        zDiff[j] = zPos[j] - zPos[i];
      }

      // bottom is scaling factor we divide by, we don't need to separate it out
      double bottoms[length];
#ifdef _SIMD
#pragma simd
#endif
      for (unsigned int j = 0; j < i; ++j) {
        double bottom = xDiff[j] * xDiff[j] + 
          yDiff[j] * yDiff[j] + 
          zDiff[j] * zDiff[j] + 
          softening;

        bottoms[j] = sqrt(bottom * bottom * bottom); // bottom ^(3/2) TODO: MKL
      }
#ifdef _SIMD
#pragma simd
#endif
      for (unsigned int j = i + 1; j < length; ++j) {
        double bottom = xDiff[j] * xDiff[j] + 
          yDiff[j] * yDiff[j] + 
          zDiff[j] * zDiff[j] + 
          softening;

        bottoms[j] = sqrt(bottom * bottom * bottom); // bottom ^(3/2) TODO: MKL
      }
#ifdef _SIMD
#pragma simd
#endif
      for (unsigned int j = 0; j < i; ++j) {
        xDiff[j] /= bottoms[j];
        yDiff[j] /= bottoms[j];
        zDiff[j] /= bottoms[j];
      }
#ifdef _SIMD
#pragma simd
#endif
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
#endif

#pragma omp parallel for
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

#if 0
  std::cerr << "Starting function with 1st particle: " << xPos[1] << ", " << yPos[1] << ", " << zPos[1] << std::endl;
#endif

#if 1
  xVel = &particles.velocities.xs[0];
  yVel = &particles.velocities.ys[0];
  zVel = &particles.velocities.zs[0];

  xAcc = &particles.accelerations.xs[0];
  yAcc = &particles.accelerations.ys[0];
  zAcc = &particles.accelerations.zs[0];

  masses = &particles.mass[0];

#endif
  length = particles.length;

#if 0
  //#pragma offload_transfer target(mic : 0) wait(signalTag) out(xPos, yPos, zPos\
  : length(length) alloc_if(0) free_if(0))

    //#pragma offload target(mic : 0)\
    nocopy(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses \
        : length(length) alloc_if(0) free_if(0))\
    nocopy(length, step, softening)
#endif
#if 0
#pragma offload target(mic : 0)\
    nocopy(xPos : length(length) REUSE into(xPosD))\
    nocopy(yPos : length(length) REUSE into(yPosD))\
    nocopy(zPos : length(length) REUSE into(zPosD))\
    nocopy(xVel : length(length) REUSE into(xVelD))\
    nocopy(yVel : length(length) REUSE into(yVelD))\
    nocopy(zVel : length(length) REUSE into(zVelD))\
    nocopy(xAcc : length(length) REUSE into(xAccD))\
    nocopy(yAcc : length(length) REUSE into(yAccD))\
    nocopy(zAcc : length(length) REUSE into(zAccD))\
    nocopy(masses : length(length) REUSE into(massesD))\
    nocopy(length : REUSE into(lengthD))\
    nocopy(step : REUSE into(stepD))\
    nocopy(softening: REUSE into(softeningD))
    stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses,
        length, step, softening);

#pragma offload_transfer target(mic : 0)\
  out(xPos : length(length) REUSE into(xPosD))\
  out(yPos : length(length) REUSE into(yPosD))\
  out(zPos : length(length) REUSE into(zPosD))
#endif
#if 0
#pragma offload target(mic : 0)\
  inout(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses \
      : length(length) alloc_if(1) free_if(1))\
  in(length, step, softening: alloc_if(0) free_if(0))
  stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses,
      length, step, softening);
#endif
#if 1
#pragma offload target(mic : 0)\
  in(xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(0) REUSE)\
  in(length, step, softening)\
  out(xPos, yPos, zPos : length(length) REUSE)
  stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses,
      length, step, softening);

#endif
#if 0
#pragma offload target(mic:0) out(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(0) REUSE)
  {
    xPos[1] -= 1;
    yPos[1] -= 1;
    zPos[1] -= 1;
    xVel[1] -= 1;
    yVel[1] -= 1;
    zVel[1] -= 1;
    xAcc[1] -= 1;
    yAcc[1] -= 1;
    zAcc[1] -= 1;
    masses[1] -= 1;
  }
#endif
#if 0
#pragma offload_transfer target(mic:0) out(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(length) REUSE)
  std::cerr << std::endl << "After step, Position: " << xPos[1] << ", " << yPos[1] << ", " << zPos[1] << std::endl;
  std::cerr << "Velocity: " << xVel[1] << ", " << yVel[1] << ", " << zVel[1] << std::endl;
  std::cerr << "Acceleration: " << xAcc[1] << ", " << yAcc[1] << ", " << zAcc[1] << std::endl;
  std::cerr << "Mass: " << masses[1] << std::endl << std::endl;
#endif

#if 0
  std::cerr << "Ending function with 1st particle: " << xPos[1] << ", " << yPos[1] << ", " << zPos[1] << std::endl;
#endif
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

#if 0
  std::cerr << "Initial first position: " << xPos[1] << ", " << yPos[1] << ", " << zPos[1] << std::endl;
  std::cerr << "Length: " << length << std::endl;
#endif

#if 0
  //#pragma offload_transfer target(mic : 0)\
  in(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(length) alloc_if(1) free_if(0))\
    in(length, step, softening)\
    signal(signalTag)

    //#pragma offload target(mic : 0)\
    nocopy(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses \
        : length(length) alloc_if(0) free_if(0))\
    nocopy(length, step, softening)
    stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses,
        length, step, softening);
#endif
#if 0
#pragma offload_transfer target(mic : 0)\
  in(xPos : length(length) ALLOC into(xPosD))\
  in(yPos : length(length) ALLOC into(yPosD))\
  in(zPos : length(length) ALLOC into(zPosD))\
  in(xVel : length(length) ALLOC into(xVelD))\
  in(yVel : length(length) ALLOC into(yVelD))\
  in(zVel : length(length) ALLOC into(zVelD))\
  in(xAcc : length(length) ALLOC into(xAccD))\
  in(yAcc : length(length) ALLOC into(yAccD))\
  in(zAcc : length(length) ALLOC into(zAccD))\
  in(masses : length(length) ALLOC into(massesD))\
  in(length : ALLOC into(lengthD))\
  in(step : ALLOC into(stepD))\
  in(softening : ALLOC into(softeningD))
#endif
#if 1
#pragma offload_transfer target(mic:0)\
  in(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(length) ALLOC)
#endif
#if 0
#pragma offload target(mic:0)\
  in(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(length) ALLOC)
  {
    xPos[1]++;
    yPos[1]++;
    zPos[1]++;
    xVel[1]++;
    yVel[1]++;
    zVel[1]++;
    xAcc[1]++;
    yAcc[1]++;
    zAcc[1]++;
    masses[1]++;
  }
#endif

#if 0
#pragma offload_transfer target(mic:0) out(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses : length(length) REUSE)
  std::cerr << "After initial offload, position: " << xPos[1] << ", " << yPos[1] << ", " << zPos[1] << std::endl;
  std::cerr << "Velocity: " << xVel[1] << ", " << yVel[1] << ", " << zVel[1] << std::endl;
  std::cerr << "Acceleration: " << xAcc[1] << ", " << yAcc[1] << ", " << zAcc[1] << std::endl;
  std::cerr << "mass: " << masses[1] << std::endl;
#endif

}
