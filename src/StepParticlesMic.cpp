#include "Particle.h"

__attribute__((target(mic)))
void stepFunctionMic(
    double *xPos, double *yPos, double *zPos, 
    double *xVel, double *yVel, double *zVel,
    double *xAcc, double *yAcc, double *zAcc,
    double *masses,
    unsigned int length, double step, double softening
) {

  for (unsigned int i = 0; i < length; ++i) {
    for (unsigned int j = i + 1; j < length; ++j) {
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

      xAcc[j] -= masses[i] * r_i_j_x;
      yAcc[j] -= masses[i] * r_i_j_y;
      zAcc[j] -= masses[i] * r_i_j_z;
    }
  }
  for (unsigned int i = 0; i < length; ++i) {
    // Integrate with Symplectic euler.
    // Important to update velocity and use updated velocity to update position
    //particles[i]->setVelocity(particles[i]->getVelocity() + h * particles[i]->getAcceleration());
    //particles[i]->setPosition(particles[i]->getPosition() + h * particles[i]->getVelocity());
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


#pragma offload target(mic) \
 inout(xPos : length(length)) inout(yPos : length(length)) inout(zPos : length(length))\
 inout(xVel : length(length)) inout(yVel : length(length)) inout(zVel : length(length))\
 inout(xAcc : length(length)) inout(yAcc : length(length)) inout(zAcc : length(length))\
 in(length, step, softening) in(masses : length(length))
    stepFunctionMic(xPos, yPos, zPos, xVel, yVel, zVel, xAcc, yAcc, zAcc, masses,
        length, step, softening);
 
}
