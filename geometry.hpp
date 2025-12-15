#pragma once
#include "settings.hpp"
#include "kokkos_defs.hpp"

class Geometry {
public:
    int dim; //dimension
    int nNeighbors; // number of neighbors per cell

    int nGhostLayers; // number of Ghost Layers
    int hHaloLayers; // number of Halo Layers
    int nBulkNodes; // number of bulk nodes
    int nGhostNodes1; // number of ghost nodes (without outer layer)
    int nGhostNodes2; // number of ghost nodes of outer layer
    int nNodes; // total number of nodes

    HostView<int8_t> map; // geometry map

    HostView<int> coord; // coordinates (in lattice unit)

    HostView<int> index; // cell index in sparse geometry domain

    HostView<int> neighbor_h;

    View<int> neighbor; // (3^DIM-1)*(nBulk+nGhostNodes1), neighbor connection array

    int nLevels; // number of levels

    int level; // level ID

    HostView<FloatL> gridSpacing; //grid spacing according to level ID

    int nBoundaries; // number of boundaries

    HostView<int> boundaryOffset, boundaryLength; // boundary ID offsets and boundary array lengths

    virtual void arrayAllocation() = 0;

    virtual void getCoordinates() = 0;

    virtual void getNeighbors() = 0;
};


class Geometry3D:public Geometry {
public:
    int nx, ny, nz;
    Geometry3D() { 
        dim = 3; 
        nNeighbors = 26;
    }
};

class Geometry2D :public Geometry {
public:
    int nx, ny;
    Geometry2D() {
        dim = 2;
        nNeighbors = 8;
    }
};
