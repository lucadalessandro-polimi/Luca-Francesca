# LIBRARIES PATHS

EIGEN_DIR := /home/lucad/eigen
FDAPDE_DIR := fdaPDE

SYS_INCLUDE_DIR := /usr/include
ARCH_INCLUDE_DIR := /usr/include/x86_64-linux-gnu
NLOHMANN_DIR := /usr/include/nlohmann

TRIANGLE_SRC := Meshes/Test_triangle/triangle.c
TRIANGLE_NAME := skyline
TRIANGLE_OUT := Meshes/Test_triangle/triangle_standalone
TRIANGLE_POLY := Meshes/Test_triangle/$(TRIANGLE_NAME).poly

FMESHER_SCRIPT := Meshes/Test_fmesher/star.r

# FLAGS                  
CXXFLAGS := -std=c++20 -g -march=native -DFDAPDE_NO_DEBUG 
#CXXFLAGS := g++ -fsanitize=address -g -O1
INCLUDES := -I$(EIGEN_DIR) -I$(FDAPDE_DIR) -I$(SYS_INCLUDE_DIR) -I$(ARCH_INCLUDE_DIR) -I$(NLOHMANN_DIR)

# Run all targets
.PHONY: all
all: delaunay conflict refinement triangle fmesher 

# Plot a Delaunay mesh
.PHONY: delaunay
delaunay:
	g++ $(CXXFLAGS) $(INCLUDES) -O2 -o Meshes/Test_fdaPDEmesher/main Meshes/Test_fdaPDEmesher/main.cpp
	Meshes/Test_fdaPDEmesher/main
	python3 Meshes/Test_fdaPDEmesher/plot_mesh.py

# Efficiency of Delaunay conflict graph algorithm 
.PHONY: conflict
conflict: 
	g++ $(CXXFLAGS) $(INCLUDES) -o Meshes/Test_fdaPDEmesher/conflict Meshes/Test_fdaPDEmesher/conflict_efficiency.cpp
	Meshes/Test_fdaPDEmesher/conflict
	python3 Meshes/Test_fdaPDEmesher/plot_timing.py

# Efficiency of delaunay refinement algorithm
.PHONY: refinement
refinement:
	g++ $(CXXFLAGS) $(INCLUDES) -O2 -o Meshes/Test_fdaPDEmesher/refinement Meshes/Test_fdaPDEmesher/refinement_efficiency.cpp
	./Meshes/Test_fdaPDEmesher/refinement
	python3 ./Meshes/Test_fdaPDEmesher/plot_timing.py


# Plot a Triangle mesh, and calculate Triangle efficiency
.PHONY: triangle
triangle:
	gcc $(TRIANGLE_SRC) -o $(TRIANGLE_OUT) -lm
	time $(TRIANGLE_OUT) -pq20a1 $(TRIANGLE_POLY)
	python3 Meshes/Test_triangle/plot_mesh.py $(TRIANGLE_NAME)
	Meshes/Test_triangle/timing_triangle.sh

# Plot a fmesher mesh, and calculate fmesher efficiency
.PHONY: fmesher
fmesher:
	docker run --rm -v $(PWD):/mnt rocker/geospatial:latest /bin/bash -c '\
		R --vanilla -e "install.packages(\"fmesher\", repos=\"https://cloud.r-project.org\")" && \
		Rscript /mnt/$(FMESHER_SCRIPT) && \
		Rscript /mnt/Meshes/Test_fmesher/rectangle_timing.r'	
	python3 Meshes/Test_fmesher/plot_timing.py

.PHONY: solve_PDE
solve_PDE:
	docker run --rm -v $(shell pwd):/root/progetto -ti aldoclemente/fdapde-docker /bin/bash -c '\
		mkdir -p usr/include/nlohmann && \
		curl -L https://github.com/nlohmann/json/releases/latest/download/json.hpp -o usr/include/nlohmann/json.hpp && \
		cd /root/progetto/workingdir/ && cpp=/root/progetto/workingdir/fdaPDE-cpp && core=$$cpp/fdaPDE/core && cd /root/progetto/workingdir/test && \
		g++ -w -o script depde_fire.cpp -I$$cpp -I$$core -I/usr/include/eigen3 -O2 -std=c++20 -march=native -s && \
		./script'



# Cleaning directories
.PHONY: clean
clean:
	@rm -f Meshes/Test_fdaPDEmesher/main  Meshes/Test_fdaPDEmesher/conflict Meshes/Test_fdaPDEmesher/refinement
	@rm -f Meshes/Test_fdaPDEmesher/*.csv
	@rm -f Meshes/Test_triangle/*.csv
	@rm -f $(TRIANGLE_OUT)
	@rm -f Meshes/Test_triangle/*.1.node
	@rm -f Meshes/Test_triangle/*.1.ele
	@rm -f Meshes/Test_triangle/*.1.poly
	@rm -f Meshes/Test_fmesher/*.csv
	@rm -f workingdir/test/script


