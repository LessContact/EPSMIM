#include "wave_solver.h"

#include <fstream>
#include <iostream>

WaveSolver::WaveSolver(const int nx, const int ny, const int nt, const int sx, const int sy)
#ifdef PLOT
    : NX(nx), NY(ny), NT(nt), SX(sx), SY(sy),
      hx((XB - XA) / (NX - 1)),
      hy((YB - YA) / (NY - 1)),
      tau((NX <= 1000 && NY <= 1000) ? 0.01 : 0.001),
      plotter(nx, ny, "wxt")
#else
    : NX(nx), NY(ny), NT(nt), SX(sx), SY(sy),
      hx((XB - XA) / (NX - 1)),
      hy((YB - YA) / (NY - 1)),
      tau((NX <= 1000 && NY <= 1000) ? 0.01 : 0.001)
#endif
{
    inv2hx2 = 1.0/(2.0 * hx*hx);
    inv2hy2 = 1.0/(2.0 * hy*hy);
    initializeArrays();
}

void WaveSolver::initializeArrays() {
    data.resize(NY*NX*2, 0.0);
    P.resize(NY*NX);

    // Инициализация фазовой скорости
    for (int i = 0; i < NY; ++i) {
        for (int j = 0; j < NX; ++j) {
            P[access(j,i)] = (j < NX / 2) ? 0.01 : 0.04;  // 0.1*0.1 или 0.2*0.2
        }
    }
}

void WaveSolver::saveToFile(const std::string& filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    out.write((reinterpret_cast<const char*>(data.data()) + NX*NY), NX * NY * sizeof(double));
    if (!out) {
        std::cerr << "Error: Failed to write data to file " << filename << std::endl;
        return;
    }
    out.close();
}

__inline __attribute__((always_inline)) void WaveSolver::updateWaveField(const int n) {
    const uint32_t gridStride = NX * NY;
    const double tau2 = tau*tau;

    double * __restrict__ current_data = &data[gridStride*currentGridIndex];
    double * __restrict__ next_data = &data[gridStride*nextGridIndex];
    const double * __restrict__ p_values = P.data();

    for (int i = 1; i < NY-1; ++i) {
        for (int j = 1; j < NX-1; ++j) {
            // Коэффициенты для производных по x
            const double px1 = (p_values[access(j,i-1)] + p_values[access(j,i)]) * inv2hx2;
            const double px2 = (p_values[access(j-1,i-1)] + p_values[access(j-1,i)]) * inv2hx2;
            // Коэффициенты для производных по y
            const double py1 = (p_values[access(j-1,i)] + p_values[access(j,i)]) * inv2hy2;
            const double py2 = (p_values[access(j-1,i-1)] + p_values[access(j,i-1)]) * inv2hy2;

            // Вычисление следующего временного слоя
            next_data[access(j,i)] = 2 * current_data[access(j,i)] - next_data[access(j,i)] +
                      tau2 * (
                          px1 * (current_data[access(j+1,i)] - current_data[access(j,i)]) +
                          px2 * (current_data[access(j-1,i)] - current_data[access(j,i)]) +
                          py1 * (current_data[access(j,i+1)] - current_data[access(j,i)]) +
                          py2 * (current_data[access(j,i-1)] - current_data[access(j,i)])
                      );
        }
    }
    const double t = n * tau;
    const double arg = 2 * std::numbers::pi * f0 * (t - t0);
    next_data[access(SY, SX)] += tau2*std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
}

void WaveSolver::solve() {
    const auto start = std::chrono::steady_clock::now();

    for (int n = 0; n < NT; ++n) {
        updateWaveField(n);

        currentGridIndex = ++currentGridIndex % 2;
        nextGridIndex = ++nextGridIndex % 2;

#ifdef PLOT
        std::string filename = "double" + std::string(5 - std::to_string(n).length(), '0') + std::to_string(n);
        plotter.updatePlot(&(data[NX*NY * currentGridIndex]), filename, false, NX*NY);
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
#endif
        // // Вывод прогресса каждые 10% итераций
        if (n % (NT/10) == 0) {
            std::cout << "Progress: " << (n * 100.0 / NT) << std::endl;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    totalTime = diff.count();
}