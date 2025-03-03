#ifndef WAVE_SOLVER_HPP
#define WAVE_SOLVER_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>

// #define PLOT

#ifdef PLOT
#include "../realtime_plot/GnuplotRealTime.h"
#endif

class WaveSolver {
public:
    WaveSolver(int nx, int ny, int nt, int sx, int sy);
    void saveToFile(const std::string& filename) const;
    void solve();
    size_t getTotalTime() const { return totalTime; }
    double getMaxU() const { return currentMaxU; }

private:
    // Размеры сетки
    const int NX, NY, NT;
    const double XA = 0.0, XB = 4.0;
    const double YA = 0.0, YB = 4.0;
    
    // Параметры источника
    const double f0 = 1.0;
    const double t0 = 1.5;
    const double gamma = 4.0;
    const int SX, SY;  // координаты источника

    double inv2hx2;
    double inv2hy2;

    // Шаги по пространству и времени
    const double hx, hy;
    const double tau;

    std::vector<double> data;      // массив со всеми U
    std::vector<double> P;         // фазовая скорость

    size_t totalTime = 0.0;    // полное время расчёта
    double currentMaxU = 0.0;  // текущее максимальное значение U

    // uint32_t prevGridIndex = 0;
    uint32_t currentGridIndex = 0;
    uint32_t nextGridIndex = 1;

    __inline __attribute__((always_inline)) size_t access(int x, int y) const { return x+NX*y; };
    void initializeArrays();
    void updateWaveField(int n);
#ifdef PLOT
    GnuplotRealTime plotter;
#endif
};

#endif