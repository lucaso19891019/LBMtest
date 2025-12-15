#include "ldc3dt.hpp"
#include <cmath>

LDC3D::LDC3D(int n, int nGhostLayers_): size(FloatL(n)) {
    // geometry info
    nGhostLayers = nGhostLayers_;
    nx = n; ny = n; nz = n;
    map = HostView<int8_t>("map", (nx + 2 * nGhostLayers) * (ny + 2 * nGhostLayers)
        * (nx + 2 * nGhostLayers));
    index = HostView<int>("index", (nx + 2 * nGhostLayers) * (ny + 2 * nGhostLayers)
        * (nx + 2 * nGhostLayers));

    // multi-grid info
    nLevels = 1;
    level = 0;
    gridSpacing = HostView<FloatL>("gridSpacing", nLevels);
    for (int i = 0; i < nLevels; i++)
        gridSpacing(i) = 1.0 / std::pow(2, i);

    // boundary info
    nBoundaries = 1;
    boundaryOffset = HostView<int>("boundaryOffset", nBoundaries);
    boundaryLength = HostView<int>("boundaryLength", nBoundaries);
}

void LDC3D::arrayAllocation() {
    nBulkNodes = nx * ny * nz;
    nGhostNodes1 = (nx + 2 * (nGhostLayers - 1)) * (ny + 2 * (nGhostLayers - 1))
        * (nz + 2 * (nGhostLayers - 1)) - nBulkNodes;
    nGhostNodes2 = (nx + 2 * nGhostLayers) * (ny + 2 * nGhostLayers)
        * (nz + 2 * nGhostLayers) - nBulkNodes - nGhostNodes1;
    nNodes = nBulkNodes + nGhostNodes1 + nGhostNodes2;

    coord = HostView<int>("coordinates", nNodes * dim);

    neighbor = View<int>("neighbor", (nBulkNodes + nGhostNodes1) * nNeighbors);
    neighbor_h = HostView<int>("neighbor_h", (nBulkNodes + nGhostNodes1) * nNeighbors);
}

void LDC3D::getCoordinates() {
    int countGhost1 = 0;
    int countGhost2 = 0;
    for (int k = 0; k < nz + 2 * nGhostLayers; k++) {
        for (int j = 0; j < ny + 2 * nGhostLayers; j++) {
            for (int i = 0; i < nx + 2 * nGhostLayers; i++) {
                int ii, jj, kk;
                ii = i - nGhostLayers;
                jj = j - nGhostLayers;
                kk = k - nGhostLayers;
                int cellID;
                if ((ii >= 0 && ii < nx) &&
                    (jj >= 0 && jj < ny) &&
                    (kk >= 0 && kk < nz))
                    cellID = ii + jj * nx + kk * (nx * ny);
                else if ((ii >= -nGhostLayers + 1 && ii < nx + nGhostLayers - 1) &&
                    (jj >= -nGhostLayers + 1 && jj < ny + nGhostLayers - 1) &&
                    (kk >= -nGhostLayers + 1 && kk < nz + nGhostLayers - 1)){
                    cellID = nBulkNodes + countGhost1;
                    countGhost1++;
                }
                else {
                    cellID = nBulkNodes + nGhostNodes1 + countGhost2;
                    countGhost2++;
                }

                int N = i + j * (nx + 2 * nGhostLayers)
                    + k * (nx + 2 * nGhostLayers) * (ny + 2 * nGhostLayers);
                index(N) = cellID;
                coord(INDEX(cellID, 0, nNodes, dim)) = ii;
                coord(INDEX(cellID, 1, nNodes, dim)) = jj;
                coord(INDEX(cellID, 2, nNodes, dim)) = kk;
            }
        }
    }
}

void LDC3D::getNeighbors() {
    for (int cellID = 0; cellID < nBulkNodes + nGhostNodes1; cellID++) {
        int neighborID, i, j, k;
        i = coord(INDEX(cellID, 0, nNodes, dim)) + nGhostLayers;
        j = coord(INDEX(cellID, 1, nNodes, dim)) + nGhostLayers;
        k = coord(INDEX(cellID, 2, nNodes, dim)) + nGhostLayers;

        int N = i + j * (nx + 2 * nGhostLayers)
            + k * (nx + 2 * nGhostLayers) * (ny + 2 * nGhostLayers);
        int dx = 1;
        int dy = nx + 2 * nGhostLayers;
        int dz = (nx + 2 * nGhostLayers) * (ny + 2 * nGhostLayers);

        neighbor_h(INDEX(cellID, 0, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx);

        neighbor_h(INDEX(cellID, 1, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx);

        neighbor_h(INDEX(cellID, 2, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dy);

        neighbor_h(INDEX(cellID, 3, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dy);

        neighbor_h(INDEX(cellID, 4, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dz);

        neighbor_h(INDEX(cellID, 5, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dz);

        neighbor_h(INDEX(cellID, 6, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx + dy);

        neighbor_h(INDEX(cellID, 7, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx - dy);

        neighbor_h(INDEX(cellID, 8, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx + dy);

        neighbor_h(INDEX(cellID, 9, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx - dy);

        neighbor_h(INDEX(cellID, 10, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx + dz);

        neighbor_h(INDEX(cellID, 11, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx - dz);

        neighbor_h(INDEX(cellID, 12, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx + dz);

        neighbor_h(INDEX(cellID, 13, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx - dz);

        neighbor_h(INDEX(cellID, 14, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dy + dz);

        neighbor_h(INDEX(cellID, 15, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dy - dz);

        neighbor_h(INDEX(cellID, 16, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dy + dz);

        neighbor_h(INDEX(cellID, 17, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dy - dz);

        neighbor_h(INDEX(cellID, 18, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx + dy + dz);

        neighbor_h(INDEX(cellID, 19, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx + dy - dz);

        neighbor_h(INDEX(cellID, 20, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx - dy + dz);

        neighbor_h(INDEX(cellID, 21, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N + dx - dy - dz);

        neighbor_h(INDEX(cellID, 22, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx + dy + dz);

        neighbor_h(INDEX(cellID, 23, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx + dy - dz);

        neighbor_h(INDEX(cellID, 24, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx - dy + dz);

        neighbor_h(INDEX(cellID, 25, nBulkNodes + nGhostNodes1, nNeighbors))
            = index(N - dx - dy - dz);
    }

    Kokkos::deep_copy(neighbor_h, neighbor);

}
