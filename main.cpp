#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cmath>
#include "libraries/eigen-5.0.1/Eigen/IterativeLinearSolvers"

struct Cell {
    float p;
    std::vector<float> u;
    std::vector<int> pos;
    std::string cellType; // fluid, air, or solid
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
    int tempInt; // For pressure solve

    Cell(){}

    Cell(std::vector<int> pos, std::string cellType, int layer, float p = 0.0f, std::vector<float> u = {0.0f, 0.0f, 0.0f}) : p(p), u(u), pos(pos), cellType(cellType), layer(layer)
    {}
};

struct MarkerParticle {
    std::vector<float> pos;

    MarkerParticle(float pos) : pos(pos)
    {}
};

// Equation 11 or Section 4.1 Step 1
std::unordered_map<int, Cell> grid;

float h = 1.0; // For now

float tMax = 9999999999; // Arbitrarily large number
float tMin = 0.0; // No limit on min
float kCFL = 1.0;

float uMax = 0;

std::vector<MarkerParticle> markerParticles;

int xMin = -99999999; // Arbitrarily large limits
int xMax = 99999999;
int yMin = -99999999;
int yMax = 99999999;
int zMin = -99999999;
int zMax = 99999999;

float g = 9.81;

float deltaT = 0.1;
float kinematicViscosity = 0.0;

float density = 1.0;

float pAtm = 1.0; // Change later

// Hash function from Section 3.1, Equation 10
int hashFunction(int x, int y, int z) {
    return 541 * x + 79 * y + 31 * z;
}


// Assuming z is up/down
std::vector<float> calculateViscosity(int x, int y, int z, float deltaT, float kinematicViscosity) {
    bool currentIsFluid = grid[hashFunction(x, y, z)].cellType == "fluid";
    float sumX = 0.0;

    // Optimize this later
    // Calculate Laplacian
    if (grid.count(hashFunction(x + 1, y, z))) {
        if (currentIsFluid || grid[hashFunction(x + 1, y, z)].cellType == "fluid") {
            sumX += grid[hashFunction(x + 1, y, z)].u[0];
        }
    }
    if (grid.count(hashFunction(x - 1, y, z))) {
        if (grid[hashFunction(x - 1, y, z)].cellType == "fluid" || grid.count(hashFunction(x - 2, y, z)) && grid[hashFunction(x - 2, y, z)].cellType == "fluid") {
            sumX += grid[hashFunction(x - 1, y, z)].u[0];
        }
    }
    if (grid.count(hashFunction(x, y + 1, z))) {
        if (grid[hashFunction(x, y + 1, z)].cellType == "fluid" || grid.count(hashFunction(x - 1, y + 1, z)) && grid[hashFunction(x - 1, y + 1, z)].cellType == "fluid") {
            sumX += grid[hashFunction(x, y + 1, z)].u[0];
        }
    }
    if (grid.count(hashFunction(x, y - 1, z))) {
        if (grid[hashFunction(x, y - 1, z)].cellType == "fluid" || grid.count(hashFunction(x - 1, y - 1, z)) && grid[hashFunction(x - 1, y - 1, z)].cellType == "fluid") {
            sumX += grid[hashFunction(x, y - 1, z)].u[0];
        }
    }
    if (grid.count(hashFunction(x, y, z + 1))) {
        if (grid[hashFunction(x, y, z + 1)].cellType == "fluid" || grid.count(hashFunction(x - 1, y, z + 1)) && grid[hashFunction(x - 1, y, z + 1)].cellType == "fluid") {
            sumX += grid[hashFunction(x, y, z + 1)].u[0];
        }
    }
    if (grid.count(hashFunction(x, y, z - 1))) {
        if (grid[hashFunction(x, y, z - 1)].cellType == "fluid" || grid.count(hashFunction(x - 1, y, z - 1)) && grid[hashFunction(x - 1, y, z - 1)].cellType == "fluid") {
            sumX += grid[hashFunction(x, y, z - 1)].u[0];
        }
    }

    sumX -= 6 * grid[hashFunction(x, y, z)].u[0];

    float sumY = 0.0;

    if (grid.count(hashFunction(x, y + 1, z))) {
        if (currentIsFluid || grid[hashFunction(x, y + 1, z)].cellType == "fluid") {
            sumY += grid[hashFunction(x, y + 1, z)].u[1];
        }
    }
    if (grid.count(hashFunction(x, y - 1, z))) {
        if (grid[hashFunction(x, y - 1, z)].cellType == "fluid" || grid.count(hashFunction(x, y - 2, z)) && grid[hashFunction(x, y - 2, z)].cellType == "fluid") {
            sumY += grid[hashFunction(x, y - 1, z)].u[1];
        }
    }
    if (grid.count(hashFunction(x + 1, y, z))) {
        if (grid[hashFunction(x + 1, y, z)].cellType == "fluid" || grid.count(hashFunction(x + 1, y - 1, z)) && grid[hashFunction(x + 1, y - 1, z)].cellType == "fluid") {
            sumY += grid[hashFunction(x + 1, y, z)].u[1];
        }
    }
    if (grid.count(hashFunction(x - 1, y, z))) {
        if (grid[hashFunction(x - 1, y, z)].cellType == "fluid" || grid.count(hashFunction(x - 1, y - 1, z)) && grid[hashFunction(x - 1, y - 1, z)].cellType == "fluid") {
            sumY += grid[hashFunction(x - 1, y, z)].u[1];
        }
    }
    if (grid.count(hashFunction(x, y, z + 1))) {
        if (grid[hashFunction(x, y, z + 1)].cellType == "fluid" || grid.count(hashFunction(x, y - 1, z + 1)) && grid[hashFunction(x, y - 1, z + 1)].cellType == "fluid") {
            sumY += grid[hashFunction(x, y, z + 1)].u[1];
        }
    }
    if (grid.count(hashFunction(x, y, z - 1))) {
        if (grid[hashFunction(x, y, z - 1)].cellType == "fluid" || grid.count(hashFunction(x, y - 1, z - 1)) && grid[hashFunction(x, y - 1, z - 1)].cellType == "fluid") {
            sumY += grid[hashFunction(x, y, z - 1)].u[1];
        }
    }

    sumY -= 6 * grid[hashFunction(x, y, z)].u[1];

    float sumZ = 0.0;

    if (grid.count(hashFunction(x, y, z + 1))) {
        if (currentIsFluid || grid[hashFunction(x, y, z + 1)].cellType == "fluid") {
            sumZ += grid[hashFunction(x, y, z + 1)].u[2];
        }
    }
    if (grid.count(hashFunction(x, y, z - 1))) {
        if (grid[hashFunction(x, y, z - 1)].cellType == "fluid" || grid.count(hashFunction(x, y, z - 2)) && grid[hashFunction(x, y, z - 2)].cellType == "fluid") {
            sumZ += grid[hashFunction(x, y, z - 1)].u[2];
        }
    }
    if (grid.count(hashFunction(x + 1, y, z))) {
        if (grid[hashFunction(x + 1, y, z)].cellType == "fluid" || grid.count(hashFunction(x + 1, y, z - 1)) && grid[hashFunction(x + 1, y, z - 1)].cellType == "fluid") {
            sumZ += grid[hashFunction(x + 1, y, z)].u[2];
        }
    }
    if (grid.count(hashFunction(x - 1, y, z))) {
        if (grid[hashFunction(x - 1, y, z)].cellType == "fluid" || grid.count(hashFunction(x - 1, y, z - 1)) && grid[hashFunction(x - 1, y, z - 1)].cellType == "fluid") {
            sumZ += grid[hashFunction(x - 1, y, z)].u[2];
        }
    }
    if (grid.count(hashFunction(x, y + 1, z))) {
        if (grid[hashFunction(x, y + 1, z)].cellType == "fluid" || grid.count(hashFunction(x, y + 1, z - 1)) && grid[hashFunction(x, y + 1, z - 1)].cellType == "fluid") {
            sumZ += grid[hashFunction(x, y + 1, z)].u[2];
        }
    }
    if (grid.count(hashFunction(x, y - 1, z))) {
        if (grid[hashFunction(x, y - 1, z)].cellType == "fluid" || grid.count(hashFunction(x, y - 1, z - 1)) && grid[hashFunction(x, y - 1, z - 1)].cellType == "fluid") {
            sumZ += grid[hashFunction(x, y - 1, z)].u[2];
        }
    }

    sumZ -= 6 * grid[hashFunction(x, y, z)].u[2];

    // Calculate added viscosity
    float addedViscosityX = deltaT * kinematicViscosity * sumX;
    float addedViscosityY = deltaT * kinematicViscosity * sumY;
    float addedViscosityZ = deltaT * kinematicViscosity * sumZ;

    return {addedViscosityX, addedViscosityY, addedViscosityZ};
}

std::vector<float> traceParticle(float x, float y, float z, float t);
std::vector<float> getVelocity(float x, float y, float z);
float getInterpolatedValue(float x, float y, float z, int index);

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
    return V;
}

float getInterpolatedValue(float x, float y, float z, int index) {
    int i = floor(x);
    int j = floor(y);
    int k = floor(z);

    Cell currentCell = grid[hashFunction(i, j, k)];

    // Assuming that currentCell exists
    float sum = (i + 1 - x) * (j + 1 - y) * (k + 1 - z) * currentCell.u[index];
    float totalWeight = (i + 1 - x) * (j + 1 - y) * (k + 1 - z);

    if (currentCell.neighbor1Exists) {
        sum += (x - i) * (j + 1 - y) * (k + 1 - z) * currentCell.neighbor1 -> u[index];
        totalWeight += (x - i) * (j + 1 - y) * (k + 1 - z);
    }

    if (currentCell.neighbor2Exists) {
        sum += (i + 1 - x) * (y - j) * (k + 1 - z) * currentCell.neighbor2 -> u[index];
        totalWeight += (i + 1 - x) * (y - j) * (k + 1 - z);
    }

    if (currentCell.neighbor3Exists) {
        sum += (x - i) * (y - j) * (k + 1 - z) * currentCell.neighbor3 -> u[index];
        totalWeight += (x - i) * (y - j) * (k + 1 - z);
    }

    if (currentCell.neighbor4Exists) {
        sum += (i + 1 - x) * (j + 1 - y) * (z - k) * currentCell.neighbor4 -> u[index];
        totalWeight += (i + 1 - x) * (j + 1 - y) * (z - k);
    }

    if (currentCell.neighbor5Exists) {
        sum += (x - i) * (j + 1 - y) * (z - k) * currentCell.neighbor5 -> u[index];
        totalWeight += (x - i) * (j + 1 - y) * (z - k);
    }

    if (currentCell.neighbor6Exists) {
        sum += (i + 1 - x) * (y - j) * (z - k) * currentCell.neighbor6 -> u[index];
        totalWeight += (i + 1 - x) * (y - j) * (z - k);
    }

    if (currentCell.neighbor7Exists) {
        sum += (x - i) * (y - j) * (z - k) * currentCell.neighbor7 -> u[index];
        totalWeight += (x - i) * (y - j) * (z - k);
    }

    return sum / totalWeight;
}

int main(void) {
    for (auto cell : grid) {
        float cellVelocity = sqrt(pow(cell.second.u[0], 2) + pow(cell.second.u[1], 2) + pow(cell.second.u[2], 2));

        if (cellVelocity > uMax)
            uMax = cellVelocity;
    }

    float deltaT = kCFL * h / abs(uMax);

    // TODO: fix pointer error or something
    // Dynamic Grid Update or Section 4.1 Step 2
    for (auto cell : grid) {
        cell.second.layer = -1;
    }

    // Update cells that currently have fluid in them
    for (auto particle : markerParticles) {
        int particleCellX = int(particle.pos[0]);
        int particleCellY = int(particle.pos[1]);
        int particleCellZ = int(particle.pos[2]);

        int cellKey = hashFunction(particleCellX, particleCellY, particleCellZ);

        if (!grid.count(cellKey)) {
            if (xMin < particleCellX && particleCellX < xMax && yMin < particleCellY && particleCellY < yMax && zMin < particleCellZ && particleCellZ < zMax) {
                Cell newCell({particleCellX, particleCellY, particleCellZ}, "fluid", 0);

                grid.insert({cellKey, newCell});

                for (int i = 0; i < 7; i++) {
                    int neighborX = particleCellX;
                    int neighborY = particleCellY;
                    int neighborZ = particleCellZ;

                    switch (i)
                    {
                        case 0:
                            neighborX--;
                            break;
                        case 1:
                            neighborY--;
                            break;
                        case 2:
                            neighborX--;
                            neighborY--;
                            break;
                        case 3:
                            neighborZ--;
                            break;
                        case 4:
                            neighborX--;
                            neighborZ--;
                            break;
                        case 5:
                            neighborY--;
                            neighborZ--;
                            break;
                        case 6:
                            neighborX--;
                            neighborY--;
                            neighborZ--;
                            break;
                        
                        default:
                            break;
                    }

                    int neighborKey = hashFunction(neighborX, neighborY, neighborZ);

                    if (grid.count(neighborKey)) {
                        switch (i)
                        {
                        case 0:
                            grid[neighborKey].neighbor1 = &newCell;
                            grid[neighborKey].neighbor1Exists = true;
                            break;
                        case 1:
                            grid[neighborKey].neighbor2 = &newCell;
                            grid[neighborKey].neighbor2Exists = true;
                        case 2:
                            grid[neighborKey].neighbor3 = &newCell;
                            grid[neighborKey].neighbor3Exists = true;
                        case 3:
                            grid[neighborKey].neighbor4 = &newCell;
                            grid[neighborKey].neighbor4Exists = true;
                        case 4:
                            grid[neighborKey].neighbor5 = &newCell;
                            grid[neighborKey].neighbor5Exists = true;
                        case 5:
                            grid[neighborKey].neighbor6 = &newCell;
                            grid[neighborKey].neighbor6Exists = true;
                        case 7:
                            grid[neighborKey].neighbor7 = &newCell;
                            grid[neighborKey].neighbor7Exists = true;
                        
                        default:
                            break;
                        }
                    }
                }
            }
        }
        else if (true) { // TODO: Later change to check if it's part of a solid object
            grid[cellKey].cellType = "fluid";
            grid[cellKey].layer = 0;
        }
    }

    // Create a buffer zone around the fluid
    // + 1 here might be problamatic depending on meaning of paper
    for (int i = 0; i < std::max(2, int(std::ceil(kCFL))) + 1; i++) {
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

                    int neighborKey = hashFunction(neighborX, neighborY, neighborZ);

                    if (grid.count(neighborKey)) {
                        if (grid[neighborKey].layer == -1 && grid[neighborKey].cellType != "solid") {
                            grid[neighborKey].cellType = "air";
                            grid[neighborKey].layer = i;
                        }
                    }
                    else {
                        Cell newNeighborCell({neighborX, neighborY, neighborZ}, "", i);

                        if (xMin < neighborX && neighborX < xMax && yMin < neighborY && neighborY < yMax && zMin < neighborZ && neighborZ < zMax) {
                            newNeighborCell.cellType = "air";
                        }
                        else {
                            newNeighborCell.cellType = "solid";
                        }

                        grid.insert({neighborKey, newNeighborCell});

                        for (int i = 0; i < 7; i++) {
                            int affectedX = neighborX;
                            int affectedY = neighborY;
                            int affectedZ = neighborZ;

                            switch (i)
                            {
                                case 0:
                                    affectedX--;
                                    break;
                                case 1:
                                    affectedY--;
                                    break;
                                case 2:
                                    affectedX--;
                                    affectedY--;
                                    break;
                                case 3:
                                    affectedZ--;
                                    break;
                                case 4:
                                    affectedX--;
                                    affectedZ--;
                                    break;
                                case 5:
                                    affectedY--;
                                    affectedZ--;
                                    break;
                                case 6:
                                    affectedX--;
                                    affectedY--;
                                    affectedZ--;
                                    break;
                                
                                default:
                                    break;
                            }

                            int neighborKey = hashFunction(affectedX, affectedY, affectedZ);

                            if (grid.count(neighborKey)) {
                                switch (i)
                                {
                                case 0:
                                    grid[neighborKey].neighbor1 = NULL;
                                    grid[neighborKey].neighbor1Exists = false;
                                    break;
                                case 1:
                                    grid[neighborKey].neighbor2 = NULL;
                                    grid[neighborKey].neighbor2Exists = false;
                                case 2:
                                    grid[neighborKey].neighbor3 = NULL;
                                    grid[neighborKey].neighbor3Exists = false;
                                case 3:
                                    grid[neighborKey].neighbor4 = NULL;
                                    grid[neighborKey].neighbor4Exists = false;
                                case 4:
                                    grid[neighborKey].neighbor5 = NULL;
                                    grid[neighborKey].neighbor5Exists = false;
                                case 5:
                                    grid[neighborKey].neighbor6 = NULL;
                                    grid[neighborKey].neighbor6Exists = false;
                                case 7:
                                    grid[neighborKey].neighbor7 = NULL;
                                    grid[neighborKey].neighbor7Exists = false;
                                
                                default:
                                    break;
                                }
                            }
                        }
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

    // Apply convection with backwards particle trace or Section 4.1 Step 3a
    for (auto cell : grid) {
        std::vector<float> newParticlePos = traceParticle(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2], -deltaT);
        cell.second.temp = getVelocity(newParticlePos[0], newParticlePos[1], newParticlePos[2]);
    }

    for (auto cell : grid) {
        cell.second.u = cell.second.temp;
    }
    
    // Apply external forces or Section 4.1 Step 3b
    // Assuming that z is the up/down direction and that velocity components in air or solid celltypes but bordering fluids are to be affected
    for (auto cell : grid) {
        if (cell.second.cellType == "fluid") {
            cell.second.u[2] -= g;
        }
        else {
            int neighborKey = hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] - 1);

            if (grid.count(neighborKey) && grid[neighborKey].cellType == "fluid") {
                cell.second.u[2] -= g;
            }
        }
    }
    
    // Apply viscosity or Section 4.1 Step 3c
    for (auto cell : grid) {
        int x = cell.second.pos[0];
        int y = cell.second.pos[1];
        int z = cell.second.pos[2];

        cell.second.temp = calculateViscosity(x, y, z, deltaT, kinematicViscosity);
    }

    for (auto cell : grid) {
        cell.second.u[0] += cell.second.temp[0];
        cell.second.u[1] += cell.second.temp[1];
        cell.second.u[2] += cell.second.temp[2];
    }

    // Calculate pressure or Section 4.1 Step 3d
    Eigen::SparseMatrix<double> A; // row, column assuming a NxN matrix where N = # fluid cells
    Eigen::VectorXd B, P; // idk if it needs B(something)

    int i = 0;
    for (auto cell : grid) {
        if (cell.second.cellType == "fluid") {
            int nonSolidNeighbors = 0;
            int airNeighbors = 0;

            int neighborX = cell.second.pos[0];
            int neighborY = cell.second.pos[1];
            int neighborZ = cell.second.pos[2];

            for (int j = 0; j < 6; j++) {
                switch (j)
                {
                    case 0:
                        neighborX++;
                        break;
                    case 1:
                        neighborX--;
                        break;
                    case 2:
                        neighborY++;
                        break;
                    case 3:
                        neighborY--;
                        break;
                    case 4:
                        neighborZ++;
                        break;
                    case 5:
                        neighborZ--;
                        break;
                    
                    default:
                        break;
                }

                A.coeffRef(i, grid[hashFunction(neighborX, neighborY, neighborZ)].tempInt) = 0;

                if (grid[hashFunction(neighborX, neighborY, neighborZ)].cellType != "solid") {
                    nonSolidNeighbors++;

                    if (grid[hashFunction(neighborX, neighborY, neighborZ)].cellType == "fluid") {
                        A.coeffRef(i, grid[hashFunction(neighborX, neighborY, neighborZ)].tempInt) = 1;
                    }
                    else { // if air
                        airNeighbors++;
                    }
                }
            }

            // No check for existing neighbors is needed because these are fluid cells
            int currentCellX = cell.second.pos[0];
            int currentCellY = cell.second.pos[1];
            int currentCellZ = cell.second.pos[2];

            float divergence = -cell.second.u[0];

            if (grid[hashFunction(currentCellX + 1, currentCellY, currentCellZ)].cellType != "solid") {
                divergence += grid[hashFunction(currentCellX + 1, currentCellY, currentCellZ)].u[0];
            }

            divergence -= cell.second.u[1];

            if (grid[hashFunction(currentCellX, currentCellY + 1, currentCellZ)].cellType != "solid") {
                divergence += grid[hashFunction(currentCellX, currentCellY + 1, currentCellZ)].u[1];
            }

            divergence -= cell.second.u[2];

            if (grid[hashFunction(currentCellX, currentCellY, currentCellZ + 1)].cellType != "solid") {
                divergence += grid[hashFunction(currentCellX, currentCellY, currentCellZ + 1)].u[2];
            }

            B[i] = density * h / deltaT * divergence - airNeighbors * pAtm;

            cell.second.tempInt = i;
            i++;
        }
    }

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    P = cg.solve(B);

    for (auto cell : grid) {
        if (cell.second.cellType == "fluid") {
            cell.second.p = P[cell.second.tempInt];
        }
    }

    //  Apply pressure or Section 4.1 Step 3e
    for (auto cell : grid) {
        if (cell.second.cellType == "fluid") {
            if (grid[hashFunction(cell.second.pos[0] - 1, cell.second.pos[1], cell.second.pos[2])].cellType != "solid") {
                cell.second.u[0] = cell.second.u[0] - deltaT / (density * h) * (
                    cell.second.p - grid[hashFunction(cell.second.pos[0] - 1, cell.second.pos[1], cell.second.pos[2])].p);
            }
            if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1] - 1, cell.second.pos[2])].cellType != "solid") {
                cell.second.u[1] = cell.second.u[1] - deltaT / (density * h) * (
                    cell.second.p - grid[hashFunction(cell.second.pos[0], cell.second.pos[1] - 1, cell.second.pos[2])].p);
            }
            if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] - 1)].cellType != "solid") {
                cell.second.u[2] = cell.second.u[2] - deltaT / (density * h) * (
                    cell.second.p - grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] - 1)].p);
            }
        }
    }

    // Extrapolate fluid velocities into surrounding cells or Section 4.1 Step 3f
    for (auto cell : grid) {
        if (cell.second.cellType == "fluid") {
            cell.second.layer = 0;
        }
        else {
            cell.second.layer = -1;
        }
    }

    for (int i = 0; i < std::max(2, int(std::ceil(kCFL))) + 1; i++) {
        for (auto cell : grid) {
            if (cell.second.layer == -1) {
                // Assuming neighbors means the six cells directly neighboring the current cell
                int neighborCount = 0;
                std::vector<float> totalNeighborVelocity = {0, 0, 0};
                
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

                    if (grid[hashFunction(neighborX, neighborY, neighborZ)].layer == i - 1) {
                        neighborCount++;
                        
                        totalNeighborVelocity[0] += grid[hashFunction(neighborX, neighborY, neighborZ)].u[0];
                        totalNeighborVelocity[1] += grid[hashFunction(neighborX, neighborY, neighborZ)].u[1];
                        totalNeighborVelocity[2] += grid[hashFunction(neighborX, neighborY, neighborZ)].u[2];
                    }
                }

                if (neighborCount > 0) {
                    if (grid[hashFunction(cell.second.pos[0] - 1, cell.second.pos[1], cell.second.pos[2])].cellType != "fluid") {
                        cell.second.u[0] = totalNeighborVelocity[0] / neighborCount;
                    }
                    if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1] - 1, cell.second.pos[2])].cellType != "fluid") {
                        cell.second.u[1] = totalNeighborVelocity[1] / neighborCount;
                    }
                    if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] - 1)].cellType != "fluid") {
                        cell.second.u[2] = totalNeighborVelocity[2] / neighborCount;
                    }

                    cell.second.layer = i;
                }
            }
        }
    }

    // Set the velocities of solid cells or Section 4.1 Step 3g
    // Assuming boundary cells always only have one side exposed to air or fluid
    for (auto cell : grid) {
        if (cell.second.cellType == "solid") {
            if (grid[hashFunction(cell.second.pos[0] - 1, cell.second.pos[1], cell.second.pos[2])].cellType != "solid" && cell.second.u[0] > 0) {
                cell.second.u[0] = 0;
            }
            else if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1] - 1, cell.second.pos[2])].cellType != "solid" && cell.second.u[1] > 0) {
                cell.second.u[1] = 0;
            }
            else if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] - 1)].cellType != "solid" && cell.second.u[2] > 0) {
                cell.second.u[2] = 0;
            }
            else if (grid[hashFunction(cell.second.pos[0] + 1, cell.second.pos[1], cell.second.pos[2])].cellType != "solid" && grid[hashFunction(cell.second.pos[0] + 1, cell.second.pos[1], cell.second.pos[2])].u[0] < 0) {
                grid[hashFunction(cell.second.pos[0] + 1, cell.second.pos[1], cell.second.pos[2])].u[1] = 0;
            }
            else if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1] + 1, cell.second.pos[2])].cellType != "solid" && grid[hashFunction(cell.second.pos[0], cell.second.pos[1] + 1, cell.second.pos[2])].u[1] < 0) {
                grid[hashFunction(cell.second.pos[0], cell.second.pos[1] + 1, cell.second.pos[2])].u[1] = 0;
            }
            else if (grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] + 1)].cellType != "solid" && grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] + 1)].u[2] < 0) {
                grid[hashFunction(cell.second.pos[0], cell.second.pos[1], cell.second.pos[2] + 1)].u[1] = 0;
            }
        }
    }

    // Move the marker particles through the velocity field or Section 4.1 Step 4
    // TODO: check for frame boundaries after update loop is implemented
    for (auto particle : markerParticles) {
        particle.pos = traceParticle(particle.pos[0], particle.pos[1], particle.pos[2], deltaT);
    }
}
