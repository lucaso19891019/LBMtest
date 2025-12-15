#pragma once
#include "geometry.hpp"

class LDC3D :public Geometry3D {
private:
    FloatL size;

public:
    LDC3D(int n, int nGhostLayers_);

    void arrayAllocation();

    void getCoordinates();

    void getNeighbors();


};
