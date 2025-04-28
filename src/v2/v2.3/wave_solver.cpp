#include "wave_solver.h"

#include <fstream>
#include <iostream>

#include <immintrin.h>

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

    const __m256d _inv2hx2 = _mm256_set1_pd(inv2hx2);
    const __m256d _inv2hy2 = _mm256_set1_pd(inv2hy2);
    const __m256d _two = _mm256_set1_pd(2.0);
    const __m256d _tau2 = _mm256_set1_pd(tau2);
    __m256d vMax = _mm256_setzero_pd();
    __m256d maskForAbs = _mm256_set1_pd(-0.0f);

    for (int i = 1; i < NY-1; ++i) {
        // Process chunks of 4 grid points at a time
        for (int j = 1; j < NX-4; j += 4) {
            // Current grid offset
            const uint32_t curr_offset = gridStride*currentGridIndex;
            // Next grid offset
            const uint32_t next_offset = gridStride*nextGridIndex;

            // Load current values
            __m256d current = _mm256_loadu_pd(&data[curr_offset + access(j,i)]);
            __m256d next = _mm256_loadu_pd(&data[next_offset + access(j,i)]);

            // Load neighboring points for computation
            __m256d right = _mm256_loadu_pd(&data[curr_offset + access(j+1,i)]);
            __m256d left = _mm256_loadu_pd(&data[curr_offset + access(j-1,i)]);
            __m256d up = _mm256_loadu_pd(&data[curr_offset + access(j,i+1)]);
            __m256d down = _mm256_loadu_pd(&data[curr_offset + access(j,i-1)]);

            // Load P coefficients and calculate px1, px2, py1, py2 vectors
            __m256d p_j_im1 = _mm256_loadu_pd(&P[access(j,i-1)]);
            __m256d p_j_i = _mm256_loadu_pd(&P[access(j,i)]);
            __m256d p_jm1_im1 = _mm256_loadu_pd(&P[access(j-1,i-1)]);
            __m256d p_jm1_i = _mm256_loadu_pd(&P[access(j-1,i)]);

            // Calculate coefficient vectors
            __m256d px1 = _mm256_mul_pd(_mm256_add_pd(p_j_im1, p_j_i), _inv2hx2);
            __m256d px2 = _mm256_mul_pd(_mm256_add_pd(p_jm1_im1, p_jm1_i), _inv2hx2);
            __m256d py1 = _mm256_mul_pd(_mm256_add_pd(p_jm1_i, p_j_i), _inv2hy2);
            __m256d py2 = _mm256_mul_pd(_mm256_add_pd(p_jm1_im1, p_j_im1), _inv2hy2);

            // Compute differences
            __m256d diff_right = _mm256_sub_pd(right, current);
            __m256d diff_left = _mm256_sub_pd(left, current);
            __m256d diff_up = _mm256_sub_pd(up, current);
            __m256d diff_down = _mm256_sub_pd(down, current);

            const __m256d term1 = _mm256_mul_pd(diff_right, px1);
            const __m256d term2 = _mm256_mul_pd(diff_left, px2);
            const __m256d term3 = _mm256_mul_pd(diff_down, py1);
            const __m256d term4 = _mm256_mul_pd(diff_up, py2);

            const __m256d sum = _mm256_add_pd(_mm256_add_pd(term1, term2), _mm256_add_pd(term3, term4));

            __m256d result = _mm256_fmadd_pd(_two, current,
                _mm256_fmsub_pd(_tau2, sum, next));

            // Store the result
            _mm256_storeu_pd(&data[next_offset + access(j,i)], result);

            const __m256d absNewGrid = _mm256_andnot_pd(maskForAbs, result);
            vMax = _mm256_max_pd(vMax, absNewGrid);
        }

        // Handle remaining grid points with scalar code
        for (int j = NX - ((NX-1) % 4); j < NX-1; ++j) {
            // Scalar computation identical to the original code
            const double px1 = (P[access(j,i-1)] + P[access(j,i)]) * inv2hx2;
            const double px2 = (P[access(j-1,i-1)] + P[access(j-1,i)]) * inv2hx2;
            const double py1 = (P[access(j-1,i)] + P[access(j,i)]) * inv2hy2;
            const double py2 = (P[access(j-1,i-1)] + P[access(j,i-1)]) * inv2hy2;

            data[gridStride*nextGridIndex + access(j,i)] =
                2 * data[gridStride*currentGridIndex + access(j,i)] -
                data[gridStride*nextGridIndex + access(j,i)] +
                tau2 * (
                    px1 * (data[gridStride*currentGridIndex + access(j+1,i)] -
                          data[gridStride*currentGridIndex + access(j,i)]) +
                    px2 * (data[gridStride*currentGridIndex + access(j-1,i)] -
                          data[gridStride*currentGridIndex + access(j,i)]) +
                    py1 * (data[gridStride*currentGridIndex + access(j,i+1)] -
                          data[gridStride*currentGridIndex + access(j,i)]) +
                    py2 * (data[gridStride*currentGridIndex + access(j,i-1)] -
                          data[gridStride*currentGridIndex + access(j,i)])
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
        currentMaxU = std::numeric_limits<double>::lowest();
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
            std::cout << "Progress: " << (n * 100.0 / NT) << "%" << std::endl;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    totalTime = diff.count();
}