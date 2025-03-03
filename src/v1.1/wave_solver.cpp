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
    initializeArrays();
}

void WaveSolver::initializeArrays() {
    // Инициализация массивов нулями
    U_curr.resize(NY, std::vector<double>(NX, 0.0));
    U_prev.resize(NY, std::vector<double>(NX, 0.0));
    U_next.resize(NY, std::vector<double>(NX, 0.0));
    P.resize(NY, std::vector<double>(NX));

    // Инициализация фазовой скорости
    for (int i = 0; i < NY; ++i) {
        for (int j = 0; j < NX; ++j) {
            P[i][j] = (j < NX / 2) ? 0.01 : 0.04;  // 0.1*0.1 или 0.2*0.2
        }
    }
}

void WaveSolver::saveToFile(const std::string& filename) const {
    std::ofstream out(filename, std::ios::binary);
    if (!out.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (const auto& row : U_curr) {
        out.write(reinterpret_cast<const char*>(row.data()), NX * sizeof(double));
        if (!out) {
            std::cerr << "Error: Failed to write data to file " << filename << std::endl;
            return;
        }
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
//
// double WaveSolver::findMaxAbsU() const {
//     double maxU = 0.0;
//     for (int i = 0; i < NY; ++i) {
//         for (int j = 0; j < NX; ++j) {
//             maxU = std::max(maxU, std::abs(U_curr[i][j]));
//         }
//     }
//     return maxU;
// }

/*__always_inline*/ void WaveSolver::updateWaveField(const int n) {
    const double t = n * tau;
    const double arg = 2 * std::numbers::pi * f0 * (t - t0);
    for (int i = 1; i < NY-1; ++i) {
        for (int j = 1; j < NX-1; ++j) {
            // Коэффициенты для производных по x
            const double px1 = (P[i-1][j] + P[i][j]) / (2 * hx * hx);
            const double px2 = (P[i-1][j-1] + P[i][j-1]) / (2 * hx * hx);

            // Коэффициенты для производных по y
            const double py1 = (P[i][j-1] + P[i][j]) / (2 * hy * hy);
            const double py2 = (P[i-1][j-1] + P[i-1][j]) / (2 * hy * hy);

            double source;
            if (i == SY && j == SX) {
                source = std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
            } else source = 0.0;

            // Вычисление следующего временного слоя
            U_next[i][j] = 2 * U_curr[i][j] - U_prev[i][j] +
                          tau * tau * (
                              source +
                            //   calculateSource(n, i, j) +
                              px1 * (U_curr[i][j+1] - U_curr[i][j]) +
                              px2 * (U_curr[i][j-1] - U_curr[i][j]) +
                              py1 * (U_curr[i+1][j] - U_curr[i][j]) +
                              py2 * (U_curr[i-1][j] - U_curr[i][j])
                          );
        }
    }
}

void WaveSolver::solve() {
    const auto start = std::chrono::steady_clock::now();

    for (int n = 0; n < NT; ++n) {
        // updateWaveField(n);

        const double t = n * tau;
        const double arg = 2 * std::numbers::pi * f0 * (t - t0);
        for (int i = 1; i < NY-1; ++i) {
            for (int j = 1; j < NX-1; ++j) {
                // Коэффициенты для производных по x
                const double px1 = (P[i-1][j] + P[i][j]) / (2 * hx * hx);
                const double px2 = (P[i-1][j-1] + P[i][j-1]) / (2 * hx * hx);

                // Коэффициенты для производных по y
                const double py1 = (P[i][j-1] + P[i][j]) / (2 * hy * hy);
                const double py2 = (P[i-1][j-1] + P[i-1][j]) / (2 * hy * hy);

                double source;
                if (i == SY && j == SX) {
                    source = std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
                } else source = 0.0;

                // Вычисление следующего временного слоя
                U_next[i][j] = 2 * U_curr[i][j] - U_prev[i][j] +
                              tau * tau * (
                                  source +
                                  //calculateSource(n, i, j) +
                                  px1 * (U_curr[i][j+1] - U_curr[i][j]) +
                                  px2 * (U_curr[i][j-1] - U_curr[i][j]) +
                                  py1 * (U_curr[i+1][j] - U_curr[i][j]) +
                                  py2 * (U_curr[i-1][j] - U_curr[i][j])
                              );
            }
        }

        // Обновление слоёв
        U_prev.swap(U_curr);
        U_curr.swap(U_next);

#ifdef PLOT
        std::string filename = "double" + std::string(5 - std::to_string(n).length(), '0') + std::to_string(n);
        plotter.updatePlot(U_curr, filename, false);
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