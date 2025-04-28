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

__inline __attribute__((always_inline)) void WaveSolver::updateWaveField(const int n, const int i) {
    const int gridStride = NX * NY;
    const double tau2 = tau*tau;

    double * __restrict__ current_data = &data[gridStride*currentGridIndex];
    double * __restrict__ next_data = &data[gridStride*nextGridIndex];
    const double * __restrict__ p_values = P.data();

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
        // currentMaxU = std::max(currentMaxU, std::abs(next_data[access(j,i)]));
    }
    const double t = n * tau;
    const double arg = 2 * std::numbers::pi * f0 * (t - t0);
    next_data[access(SY, SX)] += tau2*std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
}

void WaveSolver::solve() {
    const std::vector<int> steps{1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 24, 25, 30, 40, 60, 80, 90, 100};
    constexpr int rerunCount = 5;

    std::vector<std::pair<int, int>> stepResults(steps.size(), std::make_pair(0, 0));
    for (int rerunIndex = 0; rerunIndex < rerunCount; ++rerunIndex) {
        for (int stepIndex = 0; stepIndex < static_cast<int>(steps.size()); ++stepIndex) {
            const auto stepsPerPass = steps[stepIndex];

            printf("Solving [%u, %u], %u steps. Steps per pass: %u\n", NX, NY, NT, stepsPerPass);

            const auto start = std::chrono::steady_clock::now();

            for (int step = 0; step < NT; step += stepsPerPass) {

                // initialize
                for (int i{1}; i < stepsPerPass && i < NT; ++i) {
                    for (int j{0}; j < i; ++j) {

                        updateWaveField(step + j, i - j);

                        std::swap(currentGridIndex, nextGridIndex);
                        // currentGridIndex = nextGridIndex;
                        // nextGridIndex = (nextGridIndex + 1) % 2;
                    }

                    if (i % 2 == 1) {
                        std::swap(currentGridIndex, nextGridIndex);
                        // currentGridIndex = nextGridIndex;
                        // nextGridIndex = (nextGridIndex + 1) % 2;
                    }
                }

                // run
                for (int i{stepsPerPass}; i < NY - 1; ++i) {
                    for (int j{0}; j < stepsPerPass && j < NT; ++j) {

                        updateWaveField(step + j, i - j);

                        std::swap(currentGridIndex, nextGridIndex);
                        // currentGridIndex = nextGridIndex;
                        // nextGridIndex = (nextGridIndex + 1) % 2;
                    }

                    if (i % 2 == 1) {
                        std::swap(currentGridIndex, nextGridIndex);
                        // currentGridIndex = nextGridIndex;
                        // nextGridIndex = (nextGridIndex + 1) % 2;
                    }
                }

                std::swap(currentGridIndex, nextGridIndex);
                // currentGridIndex = nextGridIndex;
                // nextGridIndex = (nextGridIndex + 1) % 2;

                // finalize
                for (int j{1}; j < stepsPerPass && j < NT; ++j) {
                    for (int i{0}; i < j; ++i) {

                        updateWaveField(step + j, NY - i - 1);

                        std::swap(currentGridIndex, nextGridIndex);
                        // currentGridIndex = nextGridIndex;
                        // nextGridIndex = (nextGridIndex + 1) % 2;
                    }
                }
#ifdef PLOT
                std::string filename = "double" + std::string(5 - std::to_string(step).length(), '0') + std::to_string(step);
                plotter.updatePlot(&(data[NX*NY * currentGridIndex]), filename, false, NX*NY);
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
#endif
                // // Вывод прогресса каждые 10% итераций
                if (step % (NT/10) == 0) {
                    std::cout << "Progress: " << (step * 100.0 / NT) << std::endl;
                }
            }

            const auto end = std::chrono::steady_clock::now();
            const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            totalTime = diff.count();
            printf("Grid[%u, %u] took: %lu ms.\n", NX, NY, totalTime);

            stepResults[stepIndex].first = stepsPerPass;
            stepResults[stepIndex].second += totalTime;
        }
    }

    // Open file for writing (will overwrite existing file)
    std::ofstream outfile("performance.dat");

    if (outfile.is_open()) {
        for (const auto& point : stepResults) {
            outfile << point.first << "\t" << point.second / rerunCount << "\n";
        }
        outfile.close();
        std::cout << "Data successfully written to performance.dat" << std::endl;
    }
}