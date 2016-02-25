#pragma once
#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <memory>
#include <vector>
#include <list>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <cstdint>

//#define GLEW_STATIC
//#include <GL/glew.h>

double randRange(double l, double h);
double generateGaussianNoise(double mu, double sigma);

class MatrixStack;
class Program;
class Texture;

template <typename T>
struct vectors3d {
  std::vector<T> xs;
  std::vector<T> ys;
  std::vector<T> zs;
};

template <typename T>
struct vector4d {
  T x, y, z, w;
};

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
  	
  vectors3d<double> positions;
  vectors3d<double> velocities;

  // Accumulates acceleration from multiple
  // particles
  vectors3d<double> accelerations;

  // DISPLAY ONLY
  std::vector<vector4d<float> > color;

  std::vector<float> radius;


private:
	
	// For display only
	std::vector<float> posBuf;
	std::vector<float> texBuf;
	std::vector<unsigned int> indBuf;
	uint_least32_t posBufID;
	uint_least32_t texBufID;
	uint_least32_t indBufID;
};

#endif
