#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>

float h;

struct Cell {
    float p;
    std::vector<float> u;
    std::vector<int> pos;
    std::string cellType;
    int layer;

    Cell* neighbor1; // (x + 1, y, z)
    Cell* neighbor2; // (x, y + 1, z)
    Cell* neighbor3; // (x + 1, y + 1, z)
    Cell* neighbor4; // (x, y, z + 1)
    Cell* neighbor5; // (x + 1, y, z + 1)
    Cell* neighbor6; // (x, y + 1, z + 1)
    Cell* neighbor7; // (x + 1, y + 1, z + 1)

    bool neighbor1Exists = false;
    bool neighbor2Exists = false;
    bool neighbor3Exists = false;
    bool neighbor4Exists = false;
    bool neighbor5Exists = false;
    bool neighbor6Exists = false;
    bool neighbor7Exists = false;

    std::vector<float> temp;

    Cell(std::vector<int> pos, std::string cellType, int layer, float p = 0.0f, std::vector<float> u = {0.0f, 0.0f, 0.0f}) : p(p), u(u), pos(pos), cellType(cellType), layer(layer)
    {}
};

struct CellHashFunction {
    size_t operator()(const Cell& c) const {
        return 541 * c.u[0] + 79 * c.u[1] + 31 * c.u[2];
    }
};

struct MarkerParticle {
    std::vector<float> pos;

    MarkerParticle(float pos) : pos(pos)
    {}
};

// Trace a particle from point (x, y, z) for t time using RK2
std::vector<float> traceParticle(float x, float y, float z, float t) {
    std::vector<float> V = getVelocity(x, y, z);
    V = getVelocity(x + 0.5 * t * V[0], y + 0.5 * t * V[1], y + 0.5 * t * V[2]);
    return {x + V[0] * t, y + V[1] * t, z + V[2] * t};
}

// Get the interpolated velocity at a point in space.
std::vector<float> getVelocity(float x, float y, float z) {
    std::vector<float> V;
    V[0] = getInterpolatedValue(x / h, y / h - 0.5, z / h - 0.5, 0);
    V[1] = getInterpolatedValue(x / h - 0.5, y / h, z / h - 0.5, 1);
    V[2] = getInterpolatedValue(x / h - 0.5, y / h - 0.5, z / h, 2);
}

float getInterpolatedValue(float x, float y, float z, int index) {
    
}

int main(void) {
    // Equation 11 or Section 4.1 Step 1
    std::unordered_map<int, Cell, CellHashFunction> grid;

    float h = 0.1;

    float tMax = 9999999999; // Arbitrarily large number
    float tMin = 0.0; // No limit on min
    float kCFL = 1.0;

    float uMax = 0;

    for (auto cell : grid) {
        float cellVelocity = sqrt(pow(cell.second.u[0], 2) + pow(cell.second.u[1], 2) + pow(cell.second.u[2], 2));

        if (cellVelocity > uMax)
            uMax = cellVelocity;
    }

    float deltaT = kCFL * h / abs(uMax);

    // Dynamic Grid Update or Section 4.1 Step 2
    for (auto cell : grid) {
        cell.second.layer = -1;
    }

    std::vector<Cell> markerParticles;

    int xMin = -99999999; // Arbitrarily large limits
    int xMax = 99999999;
    int yMin = -99999999;
    int yMax = 99999999; 
    int zMin = -99999999;
    int zMax = 99999999;

    // Update cells that currently have fluid in them
    for (auto particle : markerParticles) {
        int particleCellX = int(particle.pos[0]);
        int particleCellY = int(particle.pos[1]);
        int particleCellZ = int(particle.pos[2]);

        int cellKey = 541 * particleCellX + 79 * particleCellY + 31 * particleCellZ;

        if (!grid.count(cellKey)) {
            if (xMin < particleCellX < xMax && yMin < particleCellY < yMax && zMin < particleCellZ < zMax) {
                Cell newCell({particleCellX, particleCellY, particleCellZ}, "fluid", 0);

                // TODO: Add cell to hash table
            }
        }
        else if (true) { // TODO: Later change to check if it's part of a solid object
            grid[cellKey].cellType = "fluid";
            grid[cellKey].layer = 0;
        }
    }

    // Create a buffer zone around the fluid
    for (int i = 0; i < std::max(2.0f, kCFL) + 1; i++) { // TODO: Check if symbols around kCFL do anything later and fix loop
        for (auto cell : grid) {
            if ((cell.second.cellType == "fluid" || cell.second.cellType == "air") && cell.second.layer == i - 1) {
                for (int j = 0; j < 6; j++) {
                    int neighborX = cell.second.pos[0];
                    int neighborY = cell.second.pos[1];
                    int neighborZ = cell.second.pos[2];

                    switch (j)
                    {
                        case 0:
                            neighborX++;
                            break;
                        case 1:
                            neighborY++;
                            break;
                        case 2:
                            neighborZ++;
                            break;
                        case 3:
                            neighborX--;
                            break;
                        case 4:
                            neighborY--;
                            break;
                        case 5:
                            neighborZ--;
                            break;
                        
                        default:
                            break;
                    }

                    int neighborKey = 541 * neighborX + 79 * neighborY + 31 * neighborZ;

                    if (grid.count(neighborKey)) {
                        if (grid[neighborKey].layer == -1 && grid[neighborKey].cellType != "solid") {
                            grid[neighborKey].cellType = "air";
                            grid[neighborKey].layer = i;
                        }
                    }
                    else {
                        Cell newNeighborCell({neighborX, neighborY, neighborZ}, "", i);

                        if (xMin < neighborX < xMax && yMin < neighborY < yMax && zMin < neighborZ < zMax) {
                            newNeighborCell.cellType = "air";
                        }
                        else {
                            newNeighborCell.cellType = "solid";
                        }

                        // TODO: Add cell to hash table
                    }
                }
            }
        }
    }

    for (auto cell : grid) {
        if (cell.second.layer == -1) {
            grid.erase(cell.first);
        }
    }
}