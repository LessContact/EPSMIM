#include "wave_solver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

void runSimulation(const int nx, const int ny, const int nt, const int sx, const int sy) {
    std::cout << "\nRunning simulation" << std::endl;
    std::cout << "Source position: SX=" << sx << ", SY=" << sy << std::endl;

    WaveSolver solver(nx, ny, nt, sx, sy);
    solver.solve();

    // Формируем имя файла для вывода
    std::stringstream ss;
    ss << "double" << std::setfill('0') << std::setw(5) << nt << ".dat";
    std::string filename = ss.str();

    solver.saveToFile(filename);

    std::cout << "Simulation results:" << std::endl;
    std::cout << "Total computation time: " << solver.getTotalTime() << " seconds" << std::endl;
    std::cout << "Final max U value: " << solver.getMaxU() << std::endl;
    std::cout << "Results saved to: " << filename << std::endl;
}

int main() {
    constexpr int NX = 600;
    constexpr int NY = 600;
    constexpr int NT = 1500;

    std::cout << "Grid size: " << NX << "x" << NY << std::endl;
    std::cout << "Time steps: " << NT << std::endl;

    // Случай а) SX = 1, SY = 1
    runSimulation(NX, NY, NT, 1, 1);
    // runSimulation(NX, NY, NT, NX-1, NY-1);

    // getc(stdin);
    // Случай б) SX = NX-1, SY = NY-1
    // runSimulation(NX, NY, NT, NX-1, NY-1, "b");

    return 0;
}