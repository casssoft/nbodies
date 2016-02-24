#pragma once
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <memory>
#include <vector>
#include <list>
#include <cstdlib>
#include <cmath>
#include <limits>

#define GLEW_STATIC
#include <GL/glew.h>

double randRange(double l, double h);
double generateGaussianNoise(double mu, double sigma);

class MatrixStack;
class Program;
class Texture;

class Particles
{
public:
	Particles();
	virtual ~Particles();

	// OpenGL methods
	void init();
	void draw(std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> MV, unsigned int index) const;

  unsigned int length;
  std::vector<double> mass;
  	
  // 3 * length
  std::vector<double> position;

  // 3 * length
  std::vector<double> velocity;

  // 3 * length
  // Accumulates acceleration from multiple
  // particles
  std::vector<double> acceleration;

  // DISPLAY ONLY
  // 4 * length
  std::vector<float> color;

  // 1 * length
  std::vector<float> radius;


private:
	
	// For display only
	std::vector<float> posBuf;
	std::vector<float> texBuf;
	std::vector<unsigned int> indBuf;
	GLuint posBufID;
	GLuint texBufID;
	GLuint indBufID;
};

#endif
