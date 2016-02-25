#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "Particle.h"
#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Texture.h"

using namespace std;

//
// Utility functions for random numbers
//

// Uniform random double between l and h
double randRange(double l, double h)
{
	double r = rand() / double(RAND_MAX);
	return (1.0 - r) * l + r * h;
}

// Random Gaussian double with mean mu and std sigma
// http://en.wikipedia.org/wiki/Box–Muller_transform
double generateGaussianNoise(double mu, double sigma)
{
	const double epsilon = std::numeric_limits<double>::min();
	const double two_pi = 2.0*3.14159265358979323846;
	
	static double z0, z1;
	static bool generate;
	generate = !generate;
	
	if (!generate)
		return z1 * sigma + mu;
	
	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	}
	while ( u1 <= epsilon );
	
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

Particles::Particles() :
  length(0),
	posBufID(0),
	texBufID(0),
	indBufID(0)
{
}

Particles::~Particles()
{
}

void Particles::init()
{
	// Load geometry
	// 0
	posBuf.push_back(-1.0);
	posBuf.push_back(-1.0);
	posBuf.push_back(0.0);
	texBuf.push_back(0.0);
	texBuf.push_back(0.0);
	// 1
	posBuf.push_back(1.0);
	posBuf.push_back(-1.0);
	posBuf.push_back(0.0);
	texBuf.push_back(1.0);
	texBuf.push_back(0.0);
	// 2
	posBuf.push_back(-1.0);
	posBuf.push_back(1.0);
	posBuf.push_back(0.0);
	texBuf.push_back(0.0);
	texBuf.push_back(1.0);
	// 3
	posBuf.push_back(1.0);
	posBuf.push_back(1.0);
	posBuf.push_back(0.0);
	texBuf.push_back(1.0);
	texBuf.push_back(1.0);
	// indices
	indBuf.push_back(0);
	indBuf.push_back(1);
	indBuf.push_back(2);
	indBuf.push_back(3);
	
	// Send the position array to the GPU
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STATIC_DRAW);
	
	// Send the texture coordinates array to the GPU
	glGenBuffers(1, &texBufID);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glBufferData(GL_ARRAY_BUFFER, texBuf.size()*sizeof(float), &texBuf[0], GL_STATIC_DRAW);
	
	// Send the index array to the GPU
	glGenBuffers(1, &indBufID);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indBuf.size()*sizeof(unsigned int), &indBuf[0], GL_STATIC_DRAW);
	
	// Unbind the arrays
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	
	assert(glGetError() == GL_NO_ERROR);
}

void Particles::draw(shared_ptr<Program> prog, shared_ptr<MatrixStack> MV, unsigned int index) const
{
	// Enable and bind position array for drawing
	GLint h_pos = prog->getAttribute("vertPos");
	GLSL::enableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glVertexAttribPointer(h_pos, 3, GL_FLOAT, GL_FALSE, 0, 0);
	
	// Enable and bind texcoord array for drawing
	GLint h_tex = prog->getAttribute("vertTex");
	GLSL::enableVertexAttribArray(h_tex);
	glBindBuffer(GL_ARRAY_BUFFER, texBufID);
	glVertexAttribPointer(h_tex, 2, GL_FLOAT, GL_FALSE, 0, 0);
	
	// Bind index array for drawing
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indBufID);
	
	// Transformation matrix
	MV->pushMatrix();
	MV->translate(Eigen::Vector3f(positions.xs[index], positions.ys[index], positions.zs[index]));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, MV->topMatrix().data());
	MV->popMatrix();
	
	// Color and scale
	glUniform4fv(prog->getUniform("color"), 1, &(color[index].w));
	glUniform1f(prog->getUniform("radius"), radius[index]);
	
	// Draw
	glDrawElements(GL_TRIANGLE_STRIP, (int)indBuf.size(), GL_UNSIGNED_INT, 0);
	
	// Disable and unbind
	GLSL::disableVertexAttribArray(h_tex);
	GLSL::disableVertexAttribArray(h_pos);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
