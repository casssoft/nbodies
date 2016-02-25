#include <memory>
#include <vector>

#include "Particle.h"


void stepParticles(Particles &particles, double step, double softening) {

  //for (unsigned int i = 0; i < particles.length; ++i) {
  //  // Not actually acceleration
  //  // This is the acceleration accumulator
  //  //particles[i]->setAcceleration(Vector3d::Zero());
  //  particles.acceleration[i*3] = 0;
  //  particles.acceleration[i*3 + 1] = 0;
  //  particles.acceleration[i*3 + 2] = 0;
  //}

  for (unsigned int i = 0; i < particles.length; ++i) {
    for (unsigned int j = i + 1; j < particles.length; ++j) {
      // r_i_j is the distance vector
      //Vector3d r_i_j = particles[j]->getPosition() - particles[i]->getPosition();
      double r_i_j_x = particles.positions.xs[j] - particles.positions.xs[i];
      double r_i_j_y = particles.positions.ys[j] - particles.positions.ys[i];
      double r_i_j_z = particles.positions.zs[j] - particles.positions.zs[i];

      // bottom is scaling factor we divide by, we don't need to separate it out
      //double bottom = r_i_j.squaredNorm() + e2;
      double bottom = r_i_j_x * r_i_j_x + r_i_j_y * r_i_j_y + r_i_j_z * r_i_j_z + softening;

      bottom = sqrt(bottom * bottom * bottom); // bottom ^(3/2)

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
      particles.accelerations.xs[i] += particles.mass[j] * r_i_j_x;
      particles.accelerations.ys[i] += particles.mass[j] * r_i_j_y;
      particles.accelerations.zs[i] += particles.mass[j] * r_i_j_z;

      particles.accelerations.xs[j] -= particles.mass[i] * r_i_j_x;
      particles.accelerations.ys[j] -= particles.mass[i] * r_i_j_y;
      particles.accelerations.zs[j] -= particles.mass[i] * r_i_j_z;
    }
  }
  for (unsigned int i = 0; i < particles.length; ++i) {
    // Integrate with Symplectic euler.
    // Important to update velocity and use updated velocity to update position
    //particles[i]->setVelocity(particles[i]->getVelocity() + h * particles[i]->getAcceleration());
    //particles[i]->setPosition(particles[i]->getPosition() + h * particles[i]->getVelocity());
    particles.velocities.xs[i] += step * particles.accelerations.xs[i];
    particles.velocities.ys[i] += step * particles.accelerations.ys[i];
    particles.velocities.zs[i] += step * particles.accelerations.zs[i];

    // Zero out acceleration now instead because of access patterns? 
    particles.accelerations.xs[i] = 0;
    particles.accelerations.ys[i] = 0;
    particles.accelerations.zs[i] = 0;

    particles.positions.xs[i] += step * particles.velocities.xs[i];
    particles.positions.ys[i] += step * particles.velocities.ys[i];
    particles.positions.zs[i] += step * particles.velocities.zs[i];
  }

}
