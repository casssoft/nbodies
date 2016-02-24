#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "Camera.h"
#include "GLSL.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Particle.h"
#include "Texture.h"

using namespace std;
using namespace Eigen;

bool keyToggles[256] = {false}; // only for English keyboards!

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from

shared_ptr<Program> progSimple;
shared_ptr<Program> prog;
shared_ptr<Camera> camera;
Particles particles;
shared_ptr<Texture> texture;
double t, h, e2;

static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
	if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
}

static void char_callback(GLFWwindow *window, unsigned int key)
{
	keyToggles[key] = !keyToggles[key];
}

static void cursor_position_callback(GLFWwindow* window, double xmouse, double ymouse)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if(state == GLFW_PRESS) {
		camera->mouseMoved(xmouse, ymouse);
	}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	// Get the current mouse position.
	double xmouse, ymouse;
	glfwGetCursorPos(window, &xmouse, &ymouse);
	// Get current window size.
	int width, height;
	glfwGetWindowSize(window, &width, &height);
	if(action == GLFW_PRESS) {
		bool shift = mods & GLFW_MOD_SHIFT;
		bool ctrl  = mods & GLFW_MOD_CONTROL;
		bool alt   = mods & GLFW_MOD_ALT;
		camera->mouseClicked(xmouse, ymouse, shift, ctrl, alt);
	}
}

static void initGL()
{
	GLSL::checkVersion();
	
	// Set background color
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	// Enable z-buffer test
	glEnable(GL_DEPTH_TEST);
	// Enable alpha blending
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	progSimple = make_shared<Program>();
	progSimple->setShaderNames(RESOURCE_DIR + "simple_vert.glsl", RESOURCE_DIR + "simple_frag.glsl");
	progSimple->setVerbose(false); // Set this to true when debugging.
	progSimple->init();
	progSimple->addUniform("P");
	progSimple->addUniform("MV");
	
	prog = make_shared<Program>();
	prog->setVerbose(true); // Set this to true when debugging.
	prog->setShaderNames(RESOURCE_DIR + "particle_vert.glsl", RESOURCE_DIR + "particle_frag.glsl");
	prog->init();
	prog->addUniform("P");
	prog->addUniform("MV");
	prog->addAttribute("vertPos");
	prog->addAttribute("vertTex");
	prog->addUniform("radius");
	prog->addUniform("alphaTexture");
	prog->addUniform("color");
	
	texture = make_shared<Texture>();
	texture->setFilename(RESOURCE_DIR + "alpha.jpg");
	texture->init();
	
	camera = make_shared<Camera>();
	
	// Initialize OpenGL for particles.
  particles.init();
	
	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

// Sort particles by their z values in camera space
class ParticleSorter {
public:
	bool operator()(size_t i0, size_t i1) const
	{
		// Particle positions in world space
		//const Vector3d &x0 = particles[i0]->getPosition();
    double x0_x = particles.position[i0 * 3];
    double x0_y = particles.position[i0 * 3 + 1];
    double x0_z = particles.position[i0 * 3 + 2];
    double x1_x = particles.position[i1 * 3];
    double x1_y = particles.position[i1 * 3 + 1];
    double x1_z = particles.position[i1 * 3 + 2];
		// Particle positions in camera space
		float z0 = V.row(2) * Vector4f(x0_x, x0_y, x0_z, 1.0f);
		float z1 = V.row(2) * Vector4f(x1_x, x1_y, x1_z, 1.0f);
		return z0 < z1;
	}
	
	Matrix4f V; // current view matrix
};

ParticleSorter sorter;

// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
vector<size_t> sortIndices(const T &v) {
	// initialize original index locations
	vector<size_t> idx(v.length);
	for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
	// sort indexes based on comparing values in v
	sort(idx.begin(), idx.end(), sorter);
	return idx;
}

void renderGL()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	glViewport(0, 0, width, height);
	
	// Use the window size for camera.
	glfwGetWindowSize(window, &width, &height);
	camera->setAspect((float)width/(float)height);
	
	// Clear buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(keyToggles[(unsigned)'c']) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}
	if(keyToggles[(unsigned)'l']) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	
	auto P = make_shared<MatrixStack>();
	auto MV = make_shared<MatrixStack>();
	
	// Apply camera transforms
	P->pushMatrix();
	camera->applyProjectionMatrix(P);
	MV->pushMatrix();
	camera->applyViewMatrix(MV);
	// Set view matrix for the sorter
	sorter.V = MV->topMatrix();
	
	// Draw particles
	prog->bind();
	texture->bind(prog->getUniform("alphaTexture"), 0);
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, P->topMatrix().data());
	// Sort particles by Z for transparency rendering.
	// Since we don't want to modify the contents of the vector, we compute the
	// sorted indices and traverse the particles in this sorted order.
	
  auto indices = sortIndices(particles);
  for (unsigned int i = 0; i < indices.size(); ++i) {
    particles.draw(prog, MV, indices[i]);
  }
  //for(auto i : sortIndices(particles)) {
	//	particles[i]->draw(prog, MV);
	//}

	texture->unbind(0);
	prog->unbind();
	
	//////////////////////////////////////////////////////
	// Cleanup
	//////////////////////////////////////////////////////
	
	// Pop stacks
	MV->popMatrix();
	P->popMatrix();
	
	GLSL::checkError(GET_FILE_LINE);
}

void loadParticles(const char *filename)
{
	ifstream in;
	in.open(filename);
	if(!in.good()) {
		cout << "Cannot read " << filename << endl;
		return;
	}

	// 1st line:
	// <n> <h> <e2>
	unsigned int n;
	in >> n;
	in >> h;
	in >> e2;

  for (unsigned int i = 0; i < n; ++i) {
    double mass;
    double px, py, pz;
    double vx, vy, vz;
    float cr, cg, cb;
    float r;
    in >> mass >> px >> py >> pz >> vx >> vy >> vz >> cr >> cg >> cb >> r;
    particles.mass.push_back(mass);
    particles.position.push_back(px);
    particles.position.push_back(py);
    particles.position.push_back(pz);
    particles.velocity.push_back(vx);
    particles.velocity.push_back(vy);
    particles.velocity.push_back(vz);
    particles.acceleration.push_back(0);
    particles.acceleration.push_back(0);
    particles.acceleration.push_back(0);

    particles.color.push_back(cr);
    particles.color.push_back(cg);
    particles.color.push_back(cb);
    particles.color.push_back(1.0f);

    particles.radius.push_back(r);
  }
  particles.length = n;

	in.close();
	cout << "Loaded galaxy from " << filename << endl;
}


void stepParticles()
{

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
      double r_i_j_x = particles.position[j * 3] - particles.position[i * 3];
      double r_i_j_y = particles.position[j * 3 + 1] - particles.position[i * 3 + 1];
      double r_i_j_z = particles.position[j * 3 + 2] - particles.position[i * 3 + 2];
      
      // bottom is scaling factor we divide by, we don't need to separate it out
      //double bottom = r_i_j.squaredNorm() + e2;
      double bottom = r_i_j_x * r_i_j_x + r_i_j_y * r_i_j_y + r_i_j_z * r_i_j_z + e2;
      
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
      particles.acceleration[i * 3] += particles.mass[j] * r_i_j_x;
      particles.acceleration[i * 3 + 1] += particles.mass[j] * r_i_j_y;
      particles.acceleration[i * 3 + 2] += particles.mass[j] * r_i_j_z; 
      
      particles.acceleration[j * 3] -= particles.mass[i] * r_i_j_x;
      particles.acceleration[j * 3 + 1] -= particles.mass[i] * r_i_j_y; 
      particles.acceleration[j * 3 + 2] -= particles.mass[i] * r_i_j_z; 
    }
  }
  for (unsigned int i = 0; i < particles.length; ++i) {
    // Integrate with Symplectic euler.
    // Important to update velocity and use updated velocity to update position
    //particles[i]->setVelocity(particles[i]->getVelocity() + h * particles[i]->getAcceleration());
    //particles[i]->setPosition(particles[i]->getPosition() + h * particles[i]->getVelocity());
    particles.velocity[i * 3] += h * particles.acceleration[i * 3];
    particles.velocity[i * 3 + 1] += h * particles.acceleration[i * 3 + 1];
    particles.velocity[i * 3 + 2] += h * particles.acceleration[i * 3 + 2];

    // Zero out acceleration now instead because of access patterns? 
    particles.acceleration[i * 3] = 0;
    particles.acceleration[i * 3 + 1] = 0;
    particles.acceleration[i * 3 + 2] = 0;

    particles.position[i * 3] += h * particles.velocity[i * 3];
    particles.position[i * 3 + 1] += h * particles.velocity[i * 3 + 1];
    particles.position[i * 3 + 2] += h * particles.velocity[i * 3 + 2];
  }
  t += h;
}

int main(int argc, char **argv)
{
	if(argc != 2 && argc != 3) {
		// Wrong number of arguments
		cout << "Usage: Lab09 <RESOURCE_DIR> <(OPTIONAL) INPUT FILE>" << endl;
		cout << "   or: Lab09 <#steps>       <(OPTIONAL) INPUT FILE>" << endl;
		exit(0);
	}
	// Create the particles...
	if(argc == 2) {
		// ... without input file
		//createParticles();
    cout << "Please specify input file" << endl;
    exit(1);
	} else {
		// ... with input file
		loadParticles(argv[2]);
	}
  
  int steps;
  
  // Try parsing `steps`
  if (sscanf(argv[1], "%i", &steps) == 1) {
		// Success!
		cout << "Running without OpenGL for " << steps << " steps" << endl;
		// Run without OpenGL
		for(int k = 0; k < steps; ++k) {
			stepParticles();
		}
	} else {
		// `steps` could not be parsed
		cout << "Running with OpenGL" << endl;
		// Run with OpenGL until the window is closed
		RESOURCE_DIR = argv[1] + string("/");
		// Set error callback.
		glfwSetErrorCallback(error_callback);
		// Initialize the library.
		if(!glfwInit()) {
			return -1;
		}

		// Create a windowed mode window and its OpenGL context.
		window = glfwCreateWindow(640, 480, "YOUR NAME", NULL, NULL);
		if(!window) {
			glfwTerminate();
			return -1;
		}
		// Make the window's context current.
		glfwMakeContextCurrent(window);
		// Initialize GLEW.
		glewExperimental = true;
		if(glewInit() != GLEW_OK) {
			cerr << "Failed to initialize GLEW" << endl;
			return -1;
		}
		glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
		cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
		cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
		// Set vsync.
		glfwSwapInterval(1);
		// Set keyboard callback.
		glfwSetKeyCallback(window, key_callback);
		// Set char callback.
		glfwSetCharCallback(window, char_callback);
		// Set cursor position callback.
		glfwSetCursorPosCallback(window, cursor_position_callback);
		// Set mouse button callback.
		glfwSetMouseButtonCallback(window, mouse_button_callback);
		// Initialize scene.
		initGL();
		// Loop until the user closes the window.
		while(!glfwWindowShouldClose(window)) {
			// Step simulation.
			stepParticles();
			// Render scene.
			renderGL();
			// Swap front and back buffers.
			glfwSwapBuffers(window);
			// Poll for and process events.
			glfwPollEvents();
		}
		// Quit program.
		glfwDestroyWindow(window);
		glfwTerminate();
	}
  cout << "First particle: \n" << particles.position[0] << "\n" <<
    particles.position[1] << "\n" <<
    particles.position[2] << "\n";
  cout << "Second particle: \n" << particles.position[1] << "\n" <<
    particles.position[2] << "\n" <<
    particles.position[3] << "\n";
	cout << "Elapsed time: " << (t*3.261539827498732e6) << " years" << endl;
	return 0;
}
