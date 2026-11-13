
/* assignment2.hpp: starter kit for ITCS 4145 Assignment 2 (Parts A and B).
 *
 * Include this header before any other header.
 * Build pthreads programs with g++ -O2 -pthread.
 * Build OpenMP programs with g++ -O2 -fopenmp.
 *
 * Contents:
 * now_seconds()       monotonic wall-clock time in seconds (Part A)
 * generate_matrix()   row-major dense matrix of small random integers
 *
 * All functions are static inline.
 */

#ifndef ASSIGNMENT2_HPP
#define ASSIGNMENT2_HPP

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>





// ---------- Timer --------------------------------------------------------

// Seconds since an arbitrary fixed point;
// differences are elapsed wall time.
static inline double now_seconds() {
    timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);

    return static_cast<double>(ts.tv_sec)
         + 1e-9 * static_cast<double>(ts.tv_nsec);
}

// ---------- Random numbers ----------------------------------------------

// xorshift64*: the same seed gives the same sequence on every machine.
struct rng_t {
    std::uint64_t state;
};

static inline void rng_seed(rng_t *r, std::uint64_t seed) {
    r->state = seed * 0x9E3779B97F4A7C15ULL + 1;
}

static inline std::uint64_t rng_next(rng_t *r) {
    r->state ^= r->state >> 12;
    r->state ^= r->state << 25;
    r->state ^= r->state >> 27;

    return r->state * 2685821657736338717ULL;
}

// Uniform double in [0, 1).
static inline double rng_unit(rng_t *r) {
    return static_cast<double>(rng_next(r) >> 11)
         * (1.0 / 9007199254740992.0);
}

// ---------- Dense matrices -----------------------------------------------

// A rows x cols matrix in row-major order.
// Element (i, j) is stored at a[i * cols + j].
// Filled with random integers from 1 to 9 stored as doubles.
static inline double *generate_matrix(int rows, int cols, std::uint64_t seed) {
    double *a = static_cast<double *>(
        std::malloc(
            static_cast<std::size_t>(rows)
            * static_cast<std::size_t>(cols)
            * sizeof(double)
        )
    );

    if (a == nullptr) {
        std::fprintf(stderr, "generate_matrix: out of memory\n");
        std::exit(1);
    }

    rng_t r;
    rng_seed(&r, seed);

    std::size_t count =
        static_cast<std::size_t>(rows)
        * static_cast<std::size_t>(cols);

    for (std::size_t i = 0; i < count; i++) {
        a[i] = static_cast<double>(
            1 + static_cast<int>(rng_next(&r) % 9)
        );
    }

    return a;
}


//start of assignment code


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


std::vector<double> timer(const double* a, const double* b, double* c, int M, int K, int N, int tileSize){
    std::vector<double> timeList; 
    loopKIJ(a, b, c, M, K, N, tileSize);

    for (int i = 0; i <10; i++){
        double start = now_seconds();
        loopKIJ(a, b, c, M, K, N, tileSize);
        double end = now_seconds();
        double elapsed = end - start;
        timeList.push_back(elapsed);
    }
    std::cout << "Result: " << c[0] << std::endl;
    return timeList;

}

void updateCSV(std::string version, std::string MKN, double rate){
    std::ofstream fin;
    fin.open("../data/matrixMult.csv", std::ios_base::app);
    fin << version<< "," <<MKN << "," << rate << std::endl;
}
double average(std::vector<double>& a){
    double avg = 0.0; 
    for(int i = 0; i<a.size(); i++){
        avg = a[i]+avg;
    }
    avg = avg/a.size();
    return avg;
}
double median(std::vector<double> a){
    return (a[4] + a[5]) / 2.0;
}
void report(std::vector<double>& timeList, std::string loop, int M, int N, int K){
    double avg = average(timeList);

    std::sort(timeList.begin(), timeList.end());

    double med = median(timeList);
    double max = *std::max_element(timeList.begin(), timeList.end());
    double min = *std::min_element(timeList.begin(), timeList.end());

    double rate = (2.0*M*K*N) / med / 1e9;

    std::cout << "Report of "  << loop <<" loop "<< ", matrix size = "<<M << " * " <<N << std::endl;
    std::cout << "Rate: " << rate << " GFLOP/s" << std::endl;
    std::cout << "Average: " << avg << " s"
              << " Median: " << med << " s"
              << " Min: " << min << " s"
              << " Max: " << max << " s"
              << std::endl;

    std::string MKN = std::to_string(M)+"x"+std::to_string(K)+"x"+std::to_string(N);
    updateCSV(loop, MKN, rate);
}

int main(){
    int M = 2048;
    int K = 2048;
    int N = 2048;

    double* a = generate_matrix(M, K, 1);
    double* b = generate_matrix(K, N, 2);
    double* c = (double*)malloc(M * N * sizeof(double));
    int tileSize[] = {0, 16, 32, 64, 128, 256};

    double* expected = (double*)malloc(M * N * sizeof(double));
    
    for (int i = 0; i < 6; i++) {
        int currentSize = tileSize[i];
        loopKIJ(a, b, expected, M, K, N, 0);
        
        std::vector<double> timeList = timer(a, b, c, M, K, N, currentSize);
        bool correct = true;

        for (int j = 0; j < M * N; j++) {
            if (c[j] != expected[j]) {
                correct = false;
                break;
            }
        }
        std::cout << (correct ? "Correct multiplication" : "incorrect multiplication") << std::endl;
        std::string name = "b=" + std::to_string(currentSize);
        report(timeList, name, M, N, K);
    }


    free(expected);
    free(a);
    free(b);
    free(c);

    M = K = N = 1024;

    a = generate_matrix(M, K, 1);
    b = generate_matrix(K, N, 2);
    c = (double*)malloc(M * N * sizeof(double));

    std::vector<double> timeList = timer(a,b,c,M,K,N, 0);
    report(timeList, "b=0", M, N, K);

    free(a);
    free(b);
    free(c);

}

#endif 
