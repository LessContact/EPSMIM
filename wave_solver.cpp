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

double WaveSolver::calculateSource(const int n, const int i, const int j) const {
    if (i == SY && j == SX) {
        const double t = n * tau;
        const double arg = 2 * std::numbers::pi * f0 * ((n * t) - t0);
        return std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
    }
    return 0.0;
}

double WaveSolver::findMaxAbsU() const {
    double maxU = 0.0;
    for (int i = 0; i < NY; ++i) {
        for (int j = 0; j < NX; ++j) {
            maxU = std::max(maxU, std::abs(U_curr[i][j]));
        }
    }
    return maxU;
}

// void WaveSolver::updateWaveField(const int n) {
//     for (int i = 1; i < NY-1; ++i) {
//         for (int j = 1; j < NX-1; ++j) {
//             // Коэффициенты для производных по x
//             const double px1 = (P[j][i-1] + P[j][i]) / (2 * hx * hx);
//             const double px2 = (P[j-1][i-1] + P[j-1][i]) / (2 * hx * hx);
//
//             // Коэффициенты для производных по y
//             const double py1 = (P[j-1][i] + P[j][i]) / (2 * hy * hy);
//             const double py2 = (P[j-1][i-1] + P[j][i-1]) / (2 * hy * hy);
//
//             // Вычисление следующего временного слоя
//             U_next[j][i] = 2 * U_curr[j][i] - U_prev[j][i] +
//                           tau * tau * (
//                               calculateSource(n, i, j) +
//                               px1 * (U_curr[j+1][i] - U_curr[j][i]) +
//                               px2 * (U_curr[j-1][i] - U_curr[j][i]) +
//                               py1 * (U_curr[j][i+1] - U_curr[j][i]) +
//                               py2 * (U_curr[j][i-1] - U_curr[j][i])
//                           );
//         }
//     }
// }

void WaveSolver::updateWaveField(const int n) {
    for (int i = 1; i < NY-1; ++i) {
        for (int j = 1; j < NX-1; ++j) {
            // Коэффициенты для производных по x
            const double px1 = (P[i-1][j] + P[i][j]) / (2 * hx * hx);
            const double px2 = (P[i-1][j-1] + P[i][j-1]) / (2 * hx * hx);

            // Коэффициенты для производных по y
            const double py1 = (P[i][j-1] + P[i][j]) / (2 * hy * hy);
            const double py2 = (P[i-1][j-1] + P[i-1][j]) / (2 * hy * hy);

            // Вычисление следующего временного слоя
            U_next[i][j] = 2 * U_curr[i][j] - U_prev[i][j] +
                          tau * tau * (
                              calculateSource(n, i, j) +
                              px1 * (U_curr[i][j+1] - U_curr[i][j]) +
                              px2 * (U_curr[i][j-1] - U_curr[i][j]) +
                              py1 * (U_curr[i+1][j] - U_curr[i][j]) +
                              py2 * (U_curr[i-1][j] - U_curr[i][j])
                          );
        }
    }
}

// void WaveSolver::updateWaveField(const int n) {
//     const auto Flatten2D = [&](const uint32_t x, const uint32_t y) noexcept { return x + NX * y; };
//     for (uint32_t j{1}; j < NX - 1; ++j) {
//         for (uint32_t i{1}; i < NY - 1; ++i) {
//             const uint32_t currentIdx = Flatten2D(j, i);
//             const float impulse = static_cast<float>(std::tie(SX, SY) == std::tuple(j, i)) * calculateSource(n, i, j);
//             const float Px = 1 / (2 * hx * hx);
//             const float Py = 1 / (2 * hy * hy);
//             U_next[j][i] =
//                     2 * U_curr[j][i] - U_prev[j][i] +
//                     tau * tau *
//                         (impulse +
//                          ((U_curr[j + 1][i] - U_curr[j][i]) *
//                               (P[j][i-1] + P[j][i]) +
//                           (U_curr[j-1][i] - U_curr[j][i]) *
//                               (P[j-1][i-1] + P[j-1][i])) *
//                              Px +
//                          ((U_curr[j][i-1] - U_curr[j][i]) *
//                               (P[j-1][i] + P[j][i]) +
//                           (U_curr[j][i-1] - U_curr[j][i]) *
//                               (P[j-1][i-1] + P[j][i-1])) *
//                              Py);
//         }
//     }
//
//             // gridMax = std::max(gridMax, std::abs(U_next[currentIdx]));
// }

void WaveSolver::solve() {
    const auto start = std::chrono::steady_clock::now();

    for (int n = 0; n < NT; ++n) {
        updateWaveField(n);
        
        // Обновление слоёв
        U_prev = U_curr;
        U_curr = U_next;
        
        // Вычисление максимального значения
        currentMaxU = findMaxAbsU();

#ifdef PLOT
        std::string filename = "double" + std::string(5 - std::to_string(n).length(), '0') + std::to_string(n);
        plotter.updatePlot(U_curr, filename, false);
#endif
        // Вывод прогресса каждые 10% итераций
        if (n % (NT/10) == 0) {
            std::cout << "Progress: " << (n * 100.0 / NT) << "%, Max U: " << currentMaxU << std::endl;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> diff = end - start;
    totalTime = diff.count();
}