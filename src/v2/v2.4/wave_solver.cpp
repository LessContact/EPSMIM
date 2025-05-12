#include "wave_solver.h"

#include <iostream>
#include <immintrin.h>
#include <cstring>
#include <numbers>

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
    data_a.reset(static_cast<double*>(std::aligned_alloc(32, NY*NX*2*sizeof(double))));
    P_a.reset(static_cast<double*>(std::aligned_alloc(32, NY*NX*sizeof(double))));

    std::memset(data_a.get(), 0, NY*NX*2*sizeof(double));

    data = data_a.get();
    P = P_a.get();

    for (int i = 0; i < NY; ++i) {
        for (int j = 0; j < NX; ++j) {
            P[access(j,i)] = (j < NX / 2) ? 0.01 : 0.04;
        }
    }
}

// [a1, a2, a3, b0]
__inline __attribute__((always_inline)) static auto ShiftRight(const __m256d &a, const __m256d &b) {
    // step 1: t = [ a1, a2, a3, a0 ]
    __m256d t = _mm256_permute4x64_pd(a, 0b00111001);
    // step 2: b0 = [ b0, b0, b0, b0 ]
    __m256d b0 = _mm256_permute4x64_pd(b, 0b00000000);
    return _mm256_blend_pd(t, b0, 0b1000);
}
// [a3, b0, b1, b2]
__inline __attribute__((always_inline)) static auto ShiftLeft(const __m256d &a, const __m256d &b) {
    __m256d b_perm = _mm256_permute4x64_pd(b, 0b10010000);
    // 1. Indices 0, 0, 1, 2 -> [b0, b0, b1, b2]
    __m256d a_perm = _mm256_permute4x64_pd(a, 0b11111111);
    // 2. Index 3 replicated -> [a3, a3, a3, a3]
    return _mm256_blend_pd(b_perm, a_perm, 0x0001);
}

__inline __attribute__((always_inline)) void WaveSolver::updateWaveField(const int n) {
    const uint32_t gridStride = NX * NY;
    const double tau2 = tau*tau;

    const __m256d _inv2hx2 = _mm256_set1_pd(inv2hx2);
    const __m256d _inv2hy2 = _mm256_set1_pd(inv2hy2);
    const __m256d _two = _mm256_set1_pd(2.0);
    const __m256d _tau2 = _mm256_set1_pd(tau2);
    __m256d vMax = _mm256_setzero_pd();
    const __m256d maskForAbs = _mm256_set1_pd(-0.0);

    auto* GridNext = data + gridStride * nextGridIndex;
    auto* GridCurr = data + gridStride * currentGridIndex;

    const auto vSizeX = NX / 4;
    for (int i = 1; i < NY-1; ++i) {
        // first 4 elems
        int j = 1;
        for (; j < 4; ++j) {
            // Коэффициенты для производных по x
            const double px1 = (P[access(j,i-1)] + P[access(j,i)]) * inv2hx2;
            const double px2 = (P[access(j-1,i-1)] + P[access(j-1,i)]) * inv2hx2;
            // Коэффициенты для производных по y
            const double py1 = (P[access(j-1,i)] + P[access(j,i)]) * inv2hy2;
            const double py2 = (P[access(j-1,i-1)] + P[access(j,i-1)]) * inv2hy2;

            // Вычисление следующего временного слоя
            data[gridStride*nextGridIndex + access(j,i)] = 2 * data[gridStride*currentGridIndex + access(j,i)] - data[gridStride*nextGridIndex + access(j,i)] +
                          tau2 * (

                              px1 * (data[gridStride*currentGridIndex + access(j+1,i)] - data[gridStride*currentGridIndex + access(j,i)]) +
                              px2 * (data[gridStride*currentGridIndex + access(j-1,i)] - data[gridStride*currentGridIndex + access(j,i)]) +
                              py1 * (data[gridStride*currentGridIndex + access(j,i+1)] - data[gridStride*currentGridIndex + access(j,i)]) +
                              py2 * (data[gridStride*currentGridIndex + access(j,i-1)] - data[gridStride*currentGridIndex + access(j,i)])
                          );
            currentMaxU = std::max(currentMaxU, std::abs(data[gridStride*nextGridIndex + access(j,i)]));
        }

        // the big chunk
        __m256d gridPrev;
        __m256d gridCurr = _mm256_load_pd(&GridCurr[access(0u, i)]);
        __m256d gridNext = _mm256_load_pd(&GridCurr[access(j, i)]);

        __m256d phaseSpeedPrev;
        __m256d phaseSpeedBottomPrev;

        auto phaseSpeedCurr = _mm256_load_pd(&P[access(0u, i)]);
        auto phaseSpeedBottomCurr = _mm256_load_pd(&P[access(0u, i - 1)]);

        for (;j < NX - 4; j += 4) {
            gridPrev = gridCurr;
            gridCurr = gridNext;
            gridNext = _mm256_load_pd(&GridCurr[access(j + 4, i)]);

            const auto gridLeft = ShiftLeft(gridPrev, gridCurr);
            const auto gridBottom = _mm256_load_pd(&GridCurr[access(j, i - 1)]);
            const auto gridRight = ShiftRight(gridCurr, gridNext);
            const auto gridTop = _mm256_load_pd(&GridCurr[access(j, i + 1)]);

            const auto diffLeft = _mm256_sub_pd(gridLeft, gridCurr);
            const auto diffBottom = _mm256_sub_pd(gridBottom, gridCurr);
            const auto diffRight = _mm256_sub_pd(gridRight, gridCurr);
            const auto diffTop = _mm256_sub_pd(gridTop, gridCurr);

            phaseSpeedPrev = phaseSpeedCurr;
            phaseSpeedCurr = _mm256_load_pd(&P[access(j, i)]);
            const auto phaseSpeedLeft = ShiftLeft(phaseSpeedPrev, phaseSpeedCurr);

            phaseSpeedBottomPrev = phaseSpeedBottomCurr;
            phaseSpeedBottomCurr = _mm256_load_pd(&P[access(j, i - 1)]);
            const auto phaseSpeedBottomLeft = ShiftLeft(phaseSpeedBottomPrev, phaseSpeedBottomCurr);

            // Compute phase differences
            const __m256d vPx1 = _mm256_mul_pd(_mm256_add_pd(phaseSpeedBottomCurr, phaseSpeedCurr), _inv2hx2);
            const __m256d vPx2 = _mm256_mul_pd(_mm256_add_pd(phaseSpeedBottomLeft, phaseSpeedLeft), _inv2hx2);
            const __m256d vPy1 = _mm256_mul_pd(_mm256_add_pd(phaseSpeedLeft, phaseSpeedCurr), _inv2hy2);
            const __m256d vPy2 = _mm256_mul_pd(_mm256_add_pd(phaseSpeedBottomLeft, phaseSpeedBottomCurr), _inv2hy2);

            // Compute terms
            const __m256d term1 = _mm256_mul_pd(diffRight, vPx1);
            const __m256d term2 = _mm256_mul_pd(diffLeft, vPx2);
            const __m256d term3 = _mm256_mul_pd(diffBottom, vPy1);
            const __m256d term4 = _mm256_mul_pd(diffTop, vPy2);

            const __m256d vSum = _mm256_add_pd(_mm256_add_pd(term1, term2), _mm256_add_pd(term3, term4));

            const __m256d vCurrent = _mm256_load_pd(&GridNext[access(j, i)]);
            const __m256d vNewValue = _mm256_fmadd_pd(_two, gridCurr, _mm256_fmsub_pd(_tau2, vSum, vCurrent));
            _mm256_store_pd(&GridNext[access(j, i)], vNewValue);

            const __m256d absNewGrid = _mm256_andnot_pd(maskForAbs, vNewValue);
            vMax = _mm256_max_pd(vMax, absNewGrid);
        }

        // reduction horizontal (for doubles, AVX2)
        const __m128d vlow  = _mm256_castpd256_pd128(vMax);
        const __m128d vhigh = _mm256_extractf128_pd(vMax, 1);
        __m128d max128 = _mm_max_pd(vlow, vhigh);
        // shuffle high 64‐bit lane into low
        const __m128d hi = _mm_unpackhi_pd(max128, max128);
        // horizontal max
        max128 = _mm_max_sd(max128, hi);
        // extract final scalar
        currentMaxU = std::max(currentMaxU, _mm_cvtsd_f64(max128));

        // Handle remaining grid points with scalar code
        for (; j < NX - 1; ++j) {
            // Scalar computation identical to the original code
            const double px1 = (P[access(j,i-1)] + P[access(j,i)]) * inv2hx2;
            const double px2 = (P[access(j-1,i-1)] + P[access(j-1,i)]) * inv2hx2;
            const double py1 = (P[access(j-1,i)] + P[access(j,i)]) * inv2hy2;
            const double py2 = (P[access(j-1,i-1)] + P[access(j,i-1)]) * inv2hy2;

            data[gridStride*nextGridIndex + access(j,i)] =
                2 * data[gridStride*currentGridIndex + access(j,i)] -
                data[gridStride*nextGridIndex + access(j,i)] +
                tau2 * (
                    px1 * (data[gridStride*currentGridIndex + access(j+1,i)] - data[gridStride*currentGridIndex + access(j,i)]) +
                    px2 * (data[gridStride*currentGridIndex + access(j-1,i)] - data[gridStride*currentGridIndex + access(j,i)]) +
                    py1 * (data[gridStride*currentGridIndex + access(j,i+1)] - data[gridStride*currentGridIndex + access(j,i)]) +
                    py2 * (data[gridStride*currentGridIndex + access(j,i-1)] - data[gridStride*currentGridIndex + access(j,i)])
                );
            currentMaxU = std::max(currentMaxU, std::abs(data[gridStride*nextGridIndex + access(j,i)]));
        }
    }

    const double t = (n+1) * tau;
    const double arg = 2 * std::numbers::pi * f0 * (t - t0);
    data[gridStride*nextGridIndex + access(SY, SX)] += tau2*std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
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
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
#endif
        // Show progress every 10% iterations
        if (n % (NT/10) == 0) {
            std::cout << "Progress: " << (n * 100.0 / NT) << "%" << std::endl;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    totalTime = diff.count();
}
