ICPC=icpc
NVCC=nvcc

GLFW_LIBRARIES=/usr/lib64/librt.so /usr/lib64/libm.so /usr/lib64/libX11.so -lpthread\
 /usr/lib64/libXrandr.so /usr/lib64/libXinerama.so /usr/lib64/libXi.so\
 /usr/lib64/libXxf86vm.so /usr/lib64/libXcursor.so /usr/lib64/libGL.so

CXX=$(ICPC)
CXXFLAGS+=-Wall -pedantic -std=c++0x $(INC) -g
ICFLAGS+=
NVFLAGS+=-ccbin=clang++ -std=c++11
LDFLAGS+=-L$(GLFW_DIR)/release/src -lX11 -lglfw3 $(GLFW_LIBRARIES) $(GLEW_DIR)/lib/libGLEW.a -lGL
INC+=-I$(EIGEN3_INCLUDE_DIR) -I$(GLFW_DIR)/include -I$(GLEW_DIR)/include

REGULAR_STEP=src/StepParticlesSerial.cpp
MIC_STEP=src/StepParticlesMic.cpp
CUDA_STEP=src/StepParticlesCuda.cu

BASE_SOURCES=$(filter-out $(REGULAR_STEP) $(MIC_STEP), $(wildcard src/*.cpp))
BASE_OBJS=$(BASE_SOURCES:.cpp=.o)

.PHONY: clean all help

all: lab09-serial lab09-mic lab09-cuda

lab09-serial: $(BASE_OBJS) $(REGULAR_STEP:.cpp=.o)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

lab09-mic: $(BASE_OBJS) $(MIC_STEP:.cpp=.o)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

lab09-cuda: $(BASE_OBJS) $(CUDA_STEP:.cu=.o)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -lcuda -lcudart -o $@

$(MIC_STEP:.cpp=.o): $(MIC_STEP) src/Particle.h
	$(ICPC) $(CXXFLAGS) $(ICFLAGS) -c $< -o $@

$(CUDA_STEP:.cu=.o): $(CUDA_STEP) src/Particle.h
	$(NVCC) $(filter-out -std=c++0x -Wall -pedantic, $(CXXFLAGS)) $(NVFLAGS) -c $< -o $@

Camera.o: src/Camera.cpp src/Camera.h src/MatrixStack.h
GLSL.o: src/GLSL.cpp src/GLSL.h
main.o: src/main.cpp src/Camera.h src/GLSL.h src/Program.h \
 src/MatrixStack.h src/Particle.h src/Texture.h src/stepParticles.h
MatrixStack.o: src/MatrixStack.cpp src/MatrixStack.h
Particle.o: src/Particle.cpp src/Particle.h src/GLSL.h src/MatrixStack.h \
 src/Program.h src/Texture.h
Program.o: src/Program.cpp src/Program.h src/GLSL.h
stb_image.o: src/stb_image.h
stepParticlesSerial.o: src/stepParticlesSerial.cpp src/Particle.h
Texture.o: src/Texture.cpp src/Texture.h src/stb_image.h

clean:
	rm -rf src/*.o lab09-serial lab09-cuda lab09-mic
