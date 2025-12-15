#include <Kokkos_Core.hpp>
#include "geometry.hpp"
#include "ldc3dt.hpp"

int main(int argc, char** argv) {
    Kokkos::initialize(argc, argv);
    {
        Geometry3D* g;
        g = new LDC3D(256, 2);

	g->arrayAllocation();

    	g->getCoordinates();

    	g->getNeighbors();

        delete g;
    }
    Kokkos::finalize();

}

