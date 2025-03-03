#include "wave_solver.h"

#include <fstream>
#include <iostream>

WaveSolver::WaveSolver(const int nx, const int ny, const int nt, const int sx, const int sy)
#ifdef PLOT
    : NX(nx), NY(ny), NT(nt), SX(sx), SY(sy),
      hx((XB - XA) / (NX - 1)),
      hy((YB - YA) / (NY - 1)),
      tau((NX <= 1000 && NY <= 1000) ? 0.01 : 0.001),
      // inv2hx2(1.0/(2.0 * hx*hx)),
      // inv2hy2(1.0/(2.0 * hy*hy)),
      plotter(nx, ny, "wxt")
#else
    : NX(nx), NY(ny), NT(nt), SX(sx), SY(sy),
      hx((XB - XA) / (NX - 1)),
      hy((YB - YA) / (NY - 1)),
      tau((NX <= 1000 && NY <= 1000) ? 0.01 : 0.001)
      // inv2hx2(1.0/(2.0 * hx*hx)),
      // inv2hy2(1.0/(2.0 * hy*hy))
#endif
{
    inv2hx2 = 1.0/(2.0 * hx*hx);
    inv2hy2 = 1.0/(2.0 * hy*hy);
    initializeArrays();
}

void WaveSolver::initializeArrays() {
    // Инициализация массивов нулями
    // U_curr.resize(NY*NX, 0.0);
    // U_prev.resize(NY*NX, 0.0);
    // U_next.resize(NY*NX, 0.0);
    data.resize(NY*NX*3, 0.0);
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

// double WaveSolver::calculateSource(const int n, const int i, const int j) const {
//     if (i == SY && j == SX) {
//         const double t = n * tau;
//         // const double t = n * tau;
//         const double arg = 2 * std::numbers::pi * f0 * (t - t0);
//         return std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
//     }
//     return 0.0;
// }

// double WaveSolver::findMaxAbsU() const {
//     double maxU = 0.0;
//     for (int i = 0; i < NY; ++i) {
//         for (int j = 0; j < NX; ++j) {
//             maxU = std::max(maxU, std::abs(U_curr[i][j]));
//         }
//     }
//     return maxU;
// }

__inline __attribute__((always_inline)) void WaveSolver::updateWaveField(const int n) {
    const uint32_t gridStride = NX * NY;

    const double tau2 = tau*tau;

    for (int i = 1; i < NY-1; ++i) {
        for (int j = 1; j < NX-1; ++j) {
            // Коэффициенты для производных по x
            const double px1 = (P[access(j,i-1)] + P[access(j,i)]) * inv2hx2;// / (2 * hx * hx);
            const double px2 = (P[access(j-1,i-1)] + P[access(j-1,i)]) * inv2hx2;// / (2 * hx * hx);

            // Коэффициенты для производных по y
            const double py1 = (P[access(j-1,i)] + P[access(j,i)]) * inv2hy2;// / (2 * hy * hy);
            const double py2 = (P[access(j-1,i-1)] + P[access(j,i-1)]) * inv2hy2;// / (2 * hy * hy);

//            double source;
//            if (i == SY && j == SX) {
//                source = std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
//            } else source = 0.0;

            // Вычисление следующего временного слоя
            data[gridStride*nextGridIndex + access(j,i)] = 2 * data[gridStride*currentGridIndex + access(j,i)] - data[gridStride*prevGridIndex + access(j,i)] +
                          tau2 * (
//                              source +
                              //calculateSource(n, i, j) +
                              px1 * (data[gridStride*currentGridIndex + access(j+1,i)] - data[gridStride*currentGridIndex + access(j,i)]) +
                              px2 * (data[gridStride*currentGridIndex + access(j-1,i)] - data[gridStride*currentGridIndex + access(j,i)]) +
                              py1 * (data[gridStride*currentGridIndex + access(j,i+1)] - data[gridStride*currentGridIndex + access(j,i)]) +
                              py2 * (data[gridStride*currentGridIndex + access(j,i-1)] - data[gridStride*currentGridIndex + access(j,i)])
                          );
        }
    }
    const double t = n * tau;
    const double arg = 2 * std::numbers::pi * f0 * (t - t0);
    data[gridStride*nextGridIndex + access(SY, SX)] += tau2*std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
}

void WaveSolver::solve() {
    const auto start = std::chrono::steady_clock::now();

    for (int n = 0; n < NT; ++n) {
        updateWaveField(n);
        // const double t = n * tau;
        // const double arg = 2 * std::numbers::pi * f0 * (t - t0);
        // for (int i = 1; i < NY-1; ++i) {
        //     for (int j = 1; j < NX-1; ++j) {
        //         // Коэффициенты для производных по x
        //         const double px1 = (P[access(j,i-1)] + P[access(j,i)]) / (2 * hx * hx);
        //         const double px2 = (P[access(j-1,i-1)] + P[access(j-1,i)]) / (2 * hx * hx);
        //
        //         // Коэффициенты для производных по y
        //         const double py1 = (P[access(j-1,i)] + P[access(j,i)]) / (2 * hy * hy);
        //         const double py2 = (P[access(j-1,i-1)] + P[access(j,i-1)]) / (2 * hy * hy);
        //
        //         double source;
        //         if (i == SY && j == SX) {
        //             source = std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
        //         } else source = 0.0;
        //
        //         // Вычисление следующего временного слоя
        //         U_next[access(j,i)] = 2 * U_curr[access(j,i)] - U_prev[access(j,i)] +
        //                       tau * tau * (
        //                           source +
        //                           //calculateSource(n, i, j) +
        //                           px1 * (U_curr[access(j+1,i)] - U_curr[access(j,i)]) +
        //                           px2 * (U_curr[access(j-1,i)] - U_curr[access(j,i)]) +
        //                           py1 * (U_curr[access(j,i+1)] - U_curr[access(j,i)]) +
        //                           py2 * (U_curr[access(j,i-1)] - U_curr[access(j,i)])
        //                       );
        //     }
        // }

        // Обновление слоёв
        // U_prev.swap(U_curr);
        // U_curr.swap(U_next);
        prevGridIndex = ++prevGridIndex % 3;
        currentGridIndex = ++currentGridIndex % 3;
        nextGridIndex = ++nextGridIndex % 3;

#ifdef PLOT
        std::string filename = "double" + std::string(5 - std::to_string(n).length(), '0') + std::to_string(n);
        plotter.updatePlot(&(data[NX*NY * currentGridIndex]), filename, false, NX*NY);
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
#endif
        // currentMaxU = findMaxAbsU();
        // // Вывод прогресса каждые 10% итераций
        if (n % (NT/10) == 0) {
            std::cout << "Progress: " << (n * 100.0 / NT) << std::endl;//<< "%, Max U: " << currentMaxU << std::endl;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    totalTime = diff.count();
}