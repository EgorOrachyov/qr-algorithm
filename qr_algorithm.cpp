#include <chrono>
#include <numeric>
#include <vector>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <iomanip>

template<typename T = double>
struct Matrix {
    explicit Matrix(int n, T init = 0) : n(n), v(n * n, init) {}

    int n;
    std::vector<T> v;
};

template<typename T>
Matrix<T> load_matrix(const std::string &path) {
    std::ifstream file(path);
    if (!file.is_open())
        throw std::exception("Failed to open file");

    int m, n, n_values;
    file >> m >> n >> n_values;
    assert(m == n);

    Matrix<T> M(n);

    for (int k = 0; k < n_values; k++) {
        int i, j;
        T v;
        file >> i >> j >> v;
        i -= 1;
        j -= 1;
        assert(0 <= i && i < n);
        assert(0 <= j && j < n);
        M.v[i * n + j] = v;
    }

    return std::move(M);
}

template<typename T>
void print_matrix(const std::string &name, const Matrix<T> &M, int width = 8, int precision = 5) {
    std::cout << std::setprecision(precision) << "Matrix " << name << std::endl;
    for (int i = 0; i < M.n; i++) {
        for (int j = 0; j < M.n; j++) {
            std::cout << std::setw(width) << M.v[i * M.n + j] << " ";
        }
        std::cout << std::endl;
    }
}

template<typename T>
void qr_decomposition(const Matrix<T> &A /* Col-major */,
                      Matrix<T> &Q /* Col-major */,
                      Matrix<T> &R /* Row-major */) {
    assert(A.n == Q.n);
    auto n = A.n;

    Q = A;

#pragma omp parallel for
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
            R.v[i * n + k] = 0;

    for (int i = 0; i < n; i++)
        R.v[i * n + i] = 1;

    for (int i = 0; i < n; i++) {
        // Normalize q_i = q_i / || q_i ||
        T norm = 0;
        for (int k = 0; k < n; k++)
            norm += Q.v[i * n + k] * Q.v[i * n + k];
        norm = std::sqrt(norm);
        for (int k = 0; k < n; k++)
            Q.v[i * n + k] /= norm;
        R.v[i * n + i] *= norm;

        // For each next q in [i + 1, n]
        // q_j = q_j - <q_j, q_i> q_i

#pragma omp parallel for
        for (int j = i + 1; j < n; j++) {
            T q_i_dot_q_j = 0;
            for (int k = 0; k < n; k++)
                q_i_dot_q_j += Q.v[i * n + k] * Q.v[j * n + k];
            for (int k = 0; k < n; k++)
                Q.v[j * n + k] -= q_i_dot_q_j * Q.v[i * n + k];
            R.v[j * n + i] += q_i_dot_q_j;
        }
    }
}

template<typename T>
void multiply(const Matrix<T> &A /* Row-major */,
              const Matrix<T> &B /* Col-major */,
              Matrix<T> &C /* Col-major */) {
    assert(A.n == B.n);
    assert(A.n == C.n);
    auto n = A.n;

#pragma omp parallel for
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            T sum = 0;
            for (int k = 0; k < n; k++) {
                sum += A.v[i * n + k] * B.v[j * n + k];
            }
            // Store as column major
            C.v[j * n + i] = sum;
        }
    }
}

int main(int argc, const char *const *argv) {
    using clock = std::chrono::steady_clock;
    using ns = std::chrono::nanoseconds;

    int samples = 4;
    int threads_count = 8;
    int iterations = 10;
    double average = 0.0;
    double sd = 0.0;
    std::string matrix_path = "../a_2500.mtx";
    std::vector<double> times;

    omp_set_num_threads(threads_count);

    Matrix<double> A = load_matrix<double>(matrix_path); // Column major
    Matrix<double> A_next(A.n); // Column major
    Matrix<double> Q(A.n); // Column major
    Matrix<double> R(A.n); // Row major

    for (int i = 0; i < samples; i++) {
        auto time_point = clock::now();

        for (int k = 0; k < iterations; k++) {
            qr_decomposition(A, Q, R);
            multiply(R, Q, A_next);
            std::swap(A, A_next);
        }

        times.push_back(static_cast<double>(std::chrono::duration_cast<ns>(clock::now() - time_point).count()) / 1e9f);
    }

    average = std::reduce(times.begin(), times.end(), 0.0) / static_cast<double>(samples);
    sd = std::transform_reduce(times.begin(), times.end(), 0.0, std::plus<>(),
                               [=](auto x) { return (x - average) * (x - average); })
         / static_cast<double>(samples - 1);
    std::cout << average << " sec " << sd << " sec" << std::endl;

    return 0;
}