#include "affinity.h"
#include "wave_solver.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <likwid.h>

void runSimulation(const int nx, const int ny, const int nt, const int sx, const int sy) {
    std::cout << "\nRunning simulation" << std::endl;
    std::cout << "Source position: SX=" << sx << ", SY=" << sy << std::endl;

    WaveSolver solver(nx, ny, nt, sx, sy);

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_REGISTER("test");
    LIKWID_MARKER_START("test"); // for likwid benchmarking
    solver.solve();
    LIKWID_MARKER_STOP("test");
    LIKWID_MARKER_CLOSE;

    // Формируем имя файла для вывода
    std::stringstream ss;
    ss << "double" << std::setfill('0') << std::setw(5) << nt << ".dat";
    std::string filename = ss.str();

    // solver.saveToFile(filename);
    // std::cout << "Results saved to: " << filename << std::endl;

    std::cout << "Simulation results:" << std::endl;
    std::cout << "Total computation time: " << solver.getTotalTime() << " milliseconds" << std::endl;
    std::cout << "Final max U value: " << solver.getMaxU() << std::endl;
    // solver.~WaveSolver();
}

int main() {
    affinity::setAffinity(12);
    constexpr int NX = 8000;
    constexpr int NY = 8000;
    constexpr int NT = 100;

    std::cout << "Grid size: " << NX << "x" << NY << std::endl;
    std::cout << "Time steps: " << NT << std::endl;

    runSimulation(NX, NY, NT, 1, 1);

    return 0;
}