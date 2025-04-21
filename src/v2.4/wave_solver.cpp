#include "wave_solver.h"

#include <iostream>
#include <immintrin.h>
#include <cstring>
#include <numbers>
#include <csignal> // For sigaction
#include <cstdlib> // For exit

// Global variables to store the current indices for debugging
static volatile int debug_j = -1, debug_i = -1;

// Signal handler for SIGSEGV
void handle_sigsegv(int signal) {
    std::cerr << "Segmentation fault caught! Debug info: access(j, i) = access("
              << debug_j << ", " << debug_i << ")\n" <<"Is aligned: " << ((debug_j + 600 * debug_i))%32 << std::endl;
    std::exit(EXIT_FAILURE);
}

void setup_signal_handler() {
    struct sigaction sa;
    sa.sa_handler = handle_sigsegv;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGSEGV, &sa, nullptr);
}

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
    // Allocate aligned memory
    data_a.reset(static_cast<double*>(std::aligned_alloc(32, NY*NX*2*sizeof(double))));
    P_a.reset(static_cast<double*>(std::aligned_alloc(32, NY*NX*sizeof(double))));

    // Initialize with zeros
    std::memset(data_a.get(), 0, NY*NX*2*sizeof(double));

    data = data_a.get();
    P = P_a.get();

    for (int i = 0; i < NY; ++i) {
        for (int j = 0; j < NX; ++j) {
            P[access(j,i)] = (j < NX / 2) ? 0.01 : 0.04;
        }
    }
}

__inline __attribute__((always_inline)) static auto ShiftRight(const auto& a, const auto& b) {
    // 0b0010_0001 = 0x21
    // bits [1:0] = 01 pick high 128 of a (a2,a3) goes into result[127:0]
    // bits [5:4] = 10 pick low 128 of b (b0,b1) goes into result[255:128]
    return _mm256_permute2f128_pd(a, b, 0b00100001);
    //a2,a3,b0,b1
};

__inline __attribute__((always_inline)) void WaveSolver::updateWaveField(const int n) {
    const uint32_t gridStride = NX * NY;
    const double tau2 = tau*tau;

    const __m256d _inv2hx2 = _mm256_set1_pd(inv2hx2);
    const __m256d _inv2hy2 = _mm256_set1_pd(inv2hy2);
    const __m256d _two = _mm256_set1_pd(2.0);
    const __m256d _tau2 = _mm256_set1_pd(tau2);
    __m256d vMax = _mm256_setzero_pd();
    const __m256d maskForAbs = _mm256_set1_pd(-0.0f);

    auto* vGridNext = reinterpret_cast<__m256d*>(data + gridStride * nextGridIndex);
    auto* vGridCurr = reinterpret_cast<__m256d*>(data + gridStride * currentGridIndex);
    const auto* vP = reinterpret_cast<const __m256d*>(P);
    // auto volatile var = (112 + 150*NX )%32;
    // std::cout << var << std::endl;
    // vGridNext[access(112,150)] = _mm256_set1_pd(0.0);

    const auto vSizeX = NX / 4;
    for (int i = 1; i < NY-1; ++i) {
        // first 4 elems
        for (int j = 1; j < 4; ++j) {
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

        __m256d va = vGridCurr[access(0,i)];
        __m256d vb = vGridCurr[access(1,i)];
        __m256d vl = ShiftRight(va,vb);

        __m256d pca = vP[access(0,i)];
        __m256d pcb = vP[access(1,i)];
        __m256d vPL = ShiftRight(pca,pcb);

        __m256d pba = vP[access(0,i-1)];
        __m256d pbb = vP[access(1,i-1)];
        __m256d vPBL = ShiftRight(pba,pbb);

        for (int j = 1; j < vSizeX - 1; ++j) {
            const __m256d vc = vGridCurr[access(j+1,i)];
            const __m256d vr = ShiftRight(vb,vc);
            const __m256d vTop = vGridCurr[access(j,i+1)];
            const __m256d vBottom = vGridCurr[access(j,i-1)];

            const __m256d pcc = vP[access(j+1,i)];
            const __m256d vPC = ShiftRight(pcb,pcc);
            const __m256d pbc = vP[access(j+1,i-1)];
            const __m256d vPB = ShiftRight(pbb,pbc);

            // get differentials
            const __m256d vPx1 = _mm256_mul_pd(_mm256_add_pd(vPB, vPC), _inv2hx2);
            const __m256d vPx2 = _mm256_mul_pd(_mm256_add_pd(vPBL, vPL), _inv2hx2);
            const __m256d vPy1 = _mm256_mul_pd(_mm256_add_pd(vPL, vPC), _inv2hy2);
            const __m256d vPy2 = _mm256_mul_pd(_mm256_add_pd(vPBL, vPB), _inv2hy2);

            //differences
            const __m256d vDiffRight = _mm256_sub_pd(vr, vb);
            const __m256d vDiffLeft = _mm256_sub_pd(vl, vb);
            const __m256d vDiffUp = _mm256_sub_pd(vTop, vb);
            const __m256d vDiffDown = _mm256_sub_pd(vBottom, vb);

            //terms
            const __m256d vTerm1 = _mm256_mul_pd(vDiffRight, vPx1);
            const __m256d vTerm2 = _mm256_mul_pd(vDiffLeft, vPx2);
            const __m256d vTerm3 = _mm256_mul_pd(vDiffDown, vPy1);
            const __m256d vTerm4 = _mm256_mul_pd(vDiffUp, vPy2);
            //sum
            const __m256d vSum = _mm256_add_pd(_mm256_add_pd(vTerm1, vTerm2), _mm256_add_pd(vTerm3, vTerm4));
            //result
            const __m256d vCur = vGridCurr[access(j,i)];
            const __m256d vResult = _mm256_fmadd_pd(_two, vb,
                _mm256_fmsub_pd(_tau2, vSum, vCur));
            vGridNext[access(j,i)] = vResult;

            const __m256d absNewGrid = _mm256_andnot_pd(maskForAbs, vResult);
            vMax = _mm256_max_pd(vMax, absNewGrid);

            //move right
            va = vb;
            vb = vc;
            vl = vr;

            pca = pcb;
            pcb = pcc;
            vPL = vPC;

            pba = pbb;
            pbb = pbc;
            vPBL = vPB;
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
            currentMaxU = std::max(currentMaxU, std::abs(data[gridStride*nextGridIndex + access(j,i)]));
        }
    }

    const double t = n * tau;
    const double arg = 2 * std::numbers::pi * f0 * (t - t0);
    data[gridStride*nextGridIndex + access(SY, SX)] += tau2*std::exp(-(arg * arg) / (gamma * gamma)) * std::sin(arg);
}

void WaveSolver::solve() {
    // setup_signal_handler(); // Set up the signal handler

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
        // Show progress every 10% iterations
        if (n % (NT/10) == 0) {
            std::cout << "Progress: " << (n * 100.0 / NT) << "%" << std::endl;
        }
    }

    const auto end = std::chrono::steady_clock::now();
    const auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    totalTime = diff.count();
}
