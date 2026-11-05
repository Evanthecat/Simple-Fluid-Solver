#include <vector>
#include <unordered_map>

float h;

struct Cell {
    float p;
    std::vector<float> u;
    std::vector<int> pos;
    int cellType;
    int layer;

    std::vector<float> temp;

    Cell(float p, std::vector<float> u, std::vector<int> pos, int cellType, int layer) : p(p), u(u), pos(pos), cellType(cellType), layer(layer)
    {}
};

struct CellHashFunction {
    size_t operator()(const Cell& c) const {
        return 541 * c.u[0] + 79 * c.u[1] + 31 * c.u[2];
    }
};

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
    
}