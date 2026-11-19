#include <omp.h>
#include <iomanip>
#include <string>
#include <algorithm>
#include <iostream>
/* assignment2.h: starter kit for ITCS 4145 Assignment 2 (Parts A and B).
*
* Include this header before any other header, since it sets the feature
* macro that makes clock_gettime visible. Build pthreads programs with
* gcc -O2 -pthread and OpenMP programs with gcc -O2 -fopenmp. The OpenMP
* programs time with omp_get_wtime, not with now_seconds.
*
* Contents:
* now_seconds() monotonic wall-clock time in seconds (Part A)
* generate_matrix() row-major dense matrix of small random integers
*
* All functions are static inline, so each .c file that includes this
* header gets its own copy, no separate object file is needed, and unused
* functions produce no warning under -Wall.
*/
#ifndef ASSIGNMENT2_H
#define ASSIGNMENT2_H
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
/* ---------- Timer -------------------------------------------------------- */
/* Seconds since an arbitrary fixed point; differences are elapsed wall time. */
static inline double now_seconds(void) {
struct timespec ts;
clock_gettime(CLOCK_MONOTONIC, &ts);
return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}
/* ---------- Random numbers ---------------------------------------------- */
/* xorshift64*: the same seed gives the same sequence on every machine. */
typedef struct { uint64_t state; } rng_t;
static inline void rng_seed(rng_t *r, uint64_t seed) {
r->state = seed * 0x9E3779B97F4A7C15ULL + 1;
}
static inline uint64_t rng_next(rng_t *r) {
r->state ^= r->state >> 12;
r->state ^= r->state << 25;
r->state ^= r->state >> 27;
return r->state * 2685821657736338717ULL;
}
/* Uniform double in [0, 1). */
static inline double rng_unit(rng_t *r) {
return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}
/* ---------- Dense matrices ------------------------------------------------ */
/* A rows x cols matrix in row-major order, element (i, j) at a[i*cols + j],
* filled with random integers from 1 to 9 stored as doubles. Products and
* sums of such entries are exact in double for the sizes of this
* assignment, so results can be compared for exact equality.
* Exits with a message if allocation fails. */
static inline double *generate_matrix(int rows, int cols, uint64_t seed) {
double *a = (double *)malloc((size_t)rows * (size_t)cols * sizeof(double));
if (a == NULL) {
fprintf(stderr, "generate_matrix: out of memory\n");
exit(1);
}
rng_t r;
rng_seed(&r, seed);
size_t count = (size_t)rows * (size_t)cols;
for (size_t i = 0; i < count; i++)
a[i] = (double)(1 + (int)(rng_next(&r) % 9));
return a;
}
#endif /* ASSIGNMENT2_H */

void loopKIJ(const double* a, const double* b, double* c, int M, int K, int N, int tileSize) {
    if (tileSize == 0) {
        for (int i = 0; i < M * N; i++) {
            c[i] = 0.0;
        }
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    c[i * N + j] = a[i * K + k] * b[k * N + j] + c[i * N + j];
                }
            }
        }

        return;
    }
    for (int i0 = 0; i0 < M; i0 += tileSize) {
        for (int j0 = 0; j0 < N; j0 += tileSize) {
            
            int endRow = i0 + tileSize;
            if (endRow > M){ 
                endRow = M;
            }

            int endCol = j0 + tileSize;
            if (endCol > N){ 
                endCol = N;
            }
            for (int i = i0; i < endRow; i++) {
                for (int j = j0; j < endCol; j++) {
                    c[i * N + j] = 0.0;
                }
            }

            for (int k0 = 0; k0 < K; k0 += tileSize) {
                int endK = k0 + tileSize;
                if (endK > K) endK = K;
                for (int k = k0; k < endK; k++) {
                    for (int i = i0; i < endRow; i++) {
                        for (int j = j0; j < endCol; j++) {
                            c[i * N + j] =
                                a[i * K + k] * b[k * N + j] + c[i * N + j];
                        }
                    }
                }
            }
        }
    }
}

void threadsLoopKIJ(const double* a, const double* b, double* c, int M, int K, int N, int tileSize, int startRow, int endRow, int startCol, int endCol) {
    if (tileSize == 0) {
        for (int i = startRow; i < endRow; i++) {
            for (int j = startCol; j < endCol; j++) {
                c[i * N + j] = 0.0;
            }
        }

        for (int k = 0; k < K; k++) {
            for (int i = startRow; i < endRow; i++) {
                for (int j = startCol; j < endCol; j++) {
                    c[i * N + j] = a[i * K + k] * b[k * N + j] + c[i * N + j];
                }
            }
        }
        return;
    }

    for (int i0 = startRow; i0 < endRow; i0 += tileSize) {
        for (int j0 = startCol; j0 < endCol; j0 += tileSize) {
            int tileEndRow = i0 + tileSize;
            if (tileEndRow > endRow) {
                tileEndRow = endRow;
            }

            int tileendCol = j0 + tileSize;
            if (tileendCol > endCol) {
                tileendCol = endCol;
            }

            for (int i = i0; i < tileEndRow; i++) {
                for (int j = j0; j < tileendCol; j++) {
                    c[i * N + j] = 0.0;
                }
            }

            for (int k0 = 0; k0 < K; k0 += tileSize) {
                int endK = k0 + tileSize;
                if (endK > K) {
                    endK = K;
                }

                for (int k = k0; k < endK; k++) {
                    for (int i = i0; i < tileEndRow; i++) {
                        for (int j = j0; j < tileendCol; j++) {
                            c[i * N + j] =
                                a[i * K + k] * b[k * N + j] + c[i * N + j];
                        }
                    }
                }
            }
        }
    }
}

void openmpSequential(const double* a, const double* b, double* c, int M, int K, int N, int tileSize) {
    loopKIJ(a, b, c, M, K, N, tileSize);
}

void openmpRun1D(const double* a, const double* b, double* c, int M, int K, int N, int tileSize, int p) {
    #pragma omp parallel
    {
        int rank = omp_get_thread_num();

        if (omp_get_num_threads() != p) {
            std::abort();
        }

        #pragma omp for schedule(static)
        for (int t = 0; t < p; t++) {
            // With p iterations and p threads, iteration t is handled by thread t.
            int rowsPerThread = M / p;
            int startRow = t * rowsPerThread;
            int endRow = startRow + rowsPerThread;

            if (t == p - 1) {
                endRow = M;
            }

            // Use rank to check the expected static assignment.
            if (rank == t) {
                threadsLoopKIJ(a, b, c, M, K, N, tileSize,
                               startRow, endRow, 0, N);
            }
        }
    }
}

void openmpRun2D(const double* a, const double* b, double* c, int M, int K, int N, int tileSize, int pr, int pc) {
    int p = pr * pc;
    #pragma omp parallel
    {
        int rank = omp_get_thread_num();

        if (omp_get_num_threads() != p) {
            std::abort();
        }

        #pragma omp for schedule(static)
        for (int t = 0; t < p; t++) {
            int bi = t / pc;
            int bj = t % pc;

            int rowsPerThread = M / pr;
            int startRow = bi * rowsPerThread;
            int endRow = startRow + rowsPerThread;

            if (bi == pr - 1) {
                endRow = M;
            }

            int colsPerThread = N / pc;
            int startCol = bj * colsPerThread;
            int endCol = startCol + colsPerThread;

            if (bj == pc - 1) {
                endCol = N;
            }

            if (rank == t) {
                threadsLoopKIJ(a, b, c, M, K, N, tileSize,
                               startRow, endRow, startCol, endCol);
            }
        }
    }
}

template <typename Function>
std::vector<double> measure(Function function) {
    std::vector<double> times;

    function(); 

    for (int i = 0; i < 10; i++) {
        double start = omp_get_wtime();
        function();
        double end = omp_get_wtime();
        times.push_back(end - start);
    }

    return times;
}

double printReport(std::vector<double> times, const std::string& name, int M, int K, int N) {
    double total = 0.0;
    for (int i = 0; i < static_cast<int>(times.size()); i++) {
        total += times[i];
    }
    double avg = total / times.size();

    std::sort(times.begin(), times.end());

    double med = (times[4] + times[5]) / 2.0;
    double minTime = times.front();
    double maxTime = times.back();
    double rate = (2.0 * M * K * N) / med / 1e9;

    std::cout << "\n" << name << "\n";
    std::cout << "Rate: " << rate << " GFLOP/s\n";
    std::cout << "Average: " << avg << " s"
              << " Median: " << med << " s"
              << " Min: " << minTime << " s"
              << " Max: " << maxTime << " s\n";

    return med;
}

bool checkCorrect(const double* actual, const double* expected, int M, int N) {
    for (int i = 0; i < M * N; i++) {
        if (actual[i] != expected[i]) {
            return false;
        }
    }
    return true;
}

void runSequentialBaselines(int M, int K, int N,
                            const double* a, const double* b, double* c,
                            double* expected, double& medianB0,
                            double& medianB64) {
    loopKIJ(a, b, expected, M, K, N, 0);

    std::vector<double> times = measure([&]() {
        openmpSequential(a, b, c, M, K, N, 0);
    });

    std::cout << (checkCorrect(c, expected, M, N)
                  ? "Correct sequential multiplication\n"
                  : "ERROR: sequential multiplication is incorrect\n");
    medianB0 = printReport(times, "Sequential b=0", M, K, N);

    times = measure([&]() {
        openmpSequential(a, b, c, M, K, N, 64);
    });

    std::cout << (checkCorrect(c, expected, M, N)
                  ? "Correct sequential multiplication\n"
                  : "ERROR: sequential multiplication is incorrect\n");
    medianB64 = printReport(times, "Sequential b=64", M, K, N);
}

void runOneDimensional(int M, int K, int N,
                       const double* a, const double* b, double* c,
                       const double* expected, double medianB0,
                       double medianB64) {
    int threadCounts[] = {1, 2, 4, 8, 12, 18, 36};
    int tileSizes[] = {0, 64};

    std::cout << "\n===== 1D partition: " << M << " x " << N << " =====\n";

    for (int tileIndex = 0; tileIndex < 2; tileIndex++) {
        int tileSize = tileSizes[tileIndex];
        double baseline = (tileSize == 0) ? medianB0 : medianB64;

        std::cout << "\nTile size b = " << tileSize << "\n";

        for (int i = 0; i < 7; i++) {
            int p = threadCounts[i];
            omp_set_num_threads(p); 

            std::vector<double> times = measure([&]() {
                openmpRun1D(a, b, c, M, K, N, tileSize, p);
            });

            bool correct = checkCorrect(c, expected, M, N);
            std::cout << (correct ? "Correct multiplication\n"
                                  : "ERROR: incorrect multiplication\n");

            std::string name = "1D b=" + std::to_string(tileSize)
                             + ", p=" + std::to_string(p);
            double med = printReport(times, name, M, K, N);
            double speedup = medianB0 / med;
            std::cout << "Speedup against sequential untiled: "
                      << speedup << "x\n";

            if (p == 1) {
                std::cout << "T(1)/Tseq (untiled): " << med / medianB0
                          << "\n";
            }
        }
    }
}

void runTwoDimensional(int M, int K, int N,
                       const double* a, const double* b, double* c,
                       const double* expected, double medianB0,
                       bool includeColumnGrid) {
    int gridRows[] = {1, 2, 2, 3, 3, 6};
    int gridCols[] = {2, 2, 4, 4, 6, 6};

    std::cout << "\n===== 2D partition: " << M << " x " << N
              << ", b=64 =====\n";

    for (int i = 0; i < 6; i++) {
        int pr = gridRows[i];
        int pc = gridCols[i];
        int p = pr * pc;
        omp_set_num_threads(p); 

        std::vector<double> times = measure([&]() {
            openmpRun2D(a, b, c, M, K, N, 64, pr, pc);
        });

        bool correct = checkCorrect(c, expected, M, N);
        std::cout << (correct ? "Correct multiplication\n"
                              : "incorrect multiplication\n");

        std::string name = "2D b=64, grid=" + std::to_string(pr)
                         + "x" + std::to_string(pc);
        double med = printReport(times, name, M, K, N);
        std::cout << "Speedup against sequential untiled: "
                  << medianB0 / med << "x\n";
    }

    if (includeColumnGrid) {
        omp_set_num_threads(36); 
        std::vector<double> times = measure([&]() {
            openmpRun2D(a, b, c, M, K, N, 64, 1, 36);
        });

        bool correct = checkCorrect(c, expected, M, N);
        std::cout << (correct ? "Correct multiplication\n"
                              : "incorrect multiplication\n");

        double med = printReport(times, "2D b=64, grid=1x36", M, K, N);
        std::cout << "Speedup against sequential untiled: "
                  << medianB0 / med << "x\n";
    }
}

int main() {
    std::cout << std::fixed << std::setprecision(6);

    int M = 1024;
    int K = 1024;
    int N = 1024;

    double* a = generate_matrix(M, K, 1);
    double* b = generate_matrix(K, N, 2);
    double* c = static_cast<double*>(std::malloc(
        static_cast<std::size_t>(M) * N * sizeof(double)));
    double* expected = static_cast<double*>(std::malloc(
        static_cast<std::size_t>(M) * N * sizeof(double)));

    if (c == nullptr || expected == nullptr) {
        std::fprintf(stderr, "Could not allocate matrix memory\n");
        return 1;
    }

    double median1024B0 = 0.0;
    double median1024B64 = 0.0;
    runSequentialBaselines(M, K, N, a, b, c, expected, median1024B0, median1024B64);

    runOneDimensional(M, K, N, a, b, c, expected, median1024B0, median1024B64);
    runTwoDimensional(M, K, N, a, b, c, expected, median1024B0, false);

    std::free(a);
    std::free(b);
    std::free(c);
    std::free(expected);

    M = 2048;
    K = 2048;
    N = 2048;

    a = generate_matrix(M, K, 1);
    b = generate_matrix(K, N, 2);
    c = static_cast<double*>(std::malloc(
        static_cast<std::size_t>(M) * N * sizeof(double)));
    expected = static_cast<double*>(std::malloc(
        static_cast<std::size_t>(M) * N * sizeof(double)));

    if (c == nullptr || expected == nullptr) {
        std::fprintf(stderr, "Could not allocate matrix memory\n");
        return 1;
    }

    double median2048B0 = 0.0;
    double median2048B64 = 0.0;
    runSequentialBaselines(M, K, N, a, b, c, expected, median2048B0, median2048B64);

    runOneDimensional(M, K, N, a, b, c, expected, median2048B0, median2048B64);
    runTwoDimensional(M, K, N, a, b, c, expected, median2048B0, true);

    std::free(a);
    std::free(b);
    std::free(c);
    std::free(expected);

    // Q4: non-divisible dimensions, 6 x 6 grid, b=64.
    M = 1000;
    K = 1024;
    N = 1536;

    a = generate_matrix(M, K, 1);
    b = generate_matrix(K, N, 2);
    c = static_cast<double*>(std::malloc(static_cast<std::size_t>(M) * N * sizeof(double)));
    expected = static_cast<double*>(std::malloc(static_cast<std::size_t>(M) * N * sizeof(double)));

    if (c == nullptr || expected == nullptr) {
        std::fprintf(stderr, "Could not allocate matrix memory\n");
        return 1;
    }

    loopKIJ(a, b, expected, M, K, N, 0);
    omp_set_num_threads(36); 

    std::vector<double> q4Times = measure([&]() {
        openmpRun2D(a, b, c, M, K, N, 64, 6, 6);
    });

    std::cout << (checkCorrect(c, expected, M, N)
                  ? "Q4: Correct multiplication\n"
                  : "Q4: incorrect multiplication\n");
    printReport(q4Times, "Q4 2D b=64, grid=6x6", M, K, N);

    std::free(a);
    std::free(b);
    std::free(c);
    std::free(expected);

    return 0;
}
