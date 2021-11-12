#ifndef QR_DECOMPOSITION_MATH_HPP
#define QR_DECOMPOSITION_MATH_HPP

#include <vector>
#include <unordered_set>
#include <random>
#include <iostream>
#include <iomanip>

namespace qr {

    static const std::size_t PRECISION = 5;
    static const std::size_t WIDTH = 10;

    struct PairHash {
    public:
        template<typename T, typename U>
        std::size_t operator()(const std::pair<T, U> &x) const {
            return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
        }
    };

    template<typename T>
    class Vec {
    public:
        using Type = T;
        using Index = unsigned int;

        Vec() = default;

        Vec(std::size_t n, std::vector<Index> indices, std::vector<T> values)
                : n(n), indices(std::move(indices)), values(std::move(values)) {

        }

        T dot(const Vec &v) const {
            std::size_t i = 0;
            std::size_t j = 0;
            std::size_t i_end = nvals();
            std::size_t j_end = v.nvals();

            T result{};

            while (i != i_end && j != j_end) {
                if (indices[i] == v.indices[j]) {
                    result += values[i] * values[j];
                    i++;
                    j++;
                } else if (indices[i] < v.indices[j])
                    i++;
                else
                    j++;
            }

            return result;
        }

        T norm() const {
            T result{};

            for (std::size_t i = 0; i < nvals(); i++) {
                result += values[i] * values[i];
            }

            return result;
        }

        void divide(T value) {
            for (auto &v: values) {
                v /= value;
            }
        }

        Vec project(T alpha, const Vec &q) {
            std::size_t i = 0;
            std::size_t j = 0;
            std::size_t i_end = nvals();
            std::size_t j_end = q.nvals();

            std::vector<Index> res_indices;
            std::vector<T> res_values;

            while (i != i_end && j != j_end) {
                if (indices[i] == q.indices[j]) {
                    res_indices.push_back(i);
                    res_values.push_back(values[i] - alpha * q.values[j]);
                    i++;
                    j++;
                } else if (indices[i] < q.indices[j]) {
                    res_indices.push_back(i);
                    res_values.push_back(values[i]);
                    i++;
                } else {
                    res_indices.push_back(j);
                    res_values.push_back(-alpha * q.values[j]);
                    j++;
                }
            }

            return Vec(dim(), std::move(res_indices), std::move(res_values));
        }

        std::size_t dim() const {
            return n;
        }

        std::size_t nvals() const {
            return indices.size();
        }

    public:
        std::size_t n;
        std::vector<Index> indices;
        std::vector<Type> values;
    };

    template<typename Stream, typename T>
    Stream &operator<<(Stream &stream, const Vec<T> &v) {
        std::size_t i = 0;
        std::size_t j = 0;
        std::size_t i_end = v.dim();
        std::size_t j_end = v.nvals();

        while (i < i_end && j < j_end) {
            if (i == v.indices[j]) {
                stream << "[" << i << "]=" << std::setw(WIDTH) << v.values[j] << " ";
                i++;
                j++;
            } else if (i < v.indices[j]) {
                stream << "[" << i << "]=" << std::setw(WIDTH) << static_cast<T>(0) << " ";
                i++;
            }
        }

        while (i < i_end) {
            stream << "[" << i << "]=" << std::setw(WIDTH) << static_cast<T>(0) << " ";
            i++;
        }

        return stream;
    }

    template<typename T>
    class Mat {
    public:

        Mat() = default;

        Mat(std::size_t n, std::vector<Vec<T>> cols)
                : n(n), cols(std::move(cols)) {

        }

        static Mat generate_symmetric(std::size_t n, double fill_factor = 0.2, int seed = 0) {
            using p = std::pair<unsigned int, double>;
            auto max_values = n * n;
            auto values_to_gen = static_cast<std::size_t>(static_cast<double>(max_values) * fill_factor);
            std::default_random_engine engine(seed);
            std::uniform_int_distribution<unsigned int> indices(0, n - 1);
            std::uniform_real_distribution<T> values;
            std::vector<std::unordered_set<p, PairHash>> columns(n);

            for (std::size_t k = 0; k < values_to_gen; k++) {
                auto i = indices(engine);
                auto j = indices(engine);
                auto v = values(engine);

                columns[i].emplace(j, v);
                columns[j].emplace(i, v);
            }

            std::vector<Vec<T>> res_cols;
            res_cols.reserve(n);

            for (const auto &column: columns) {
                std::vector<p> tmp(column.size());
                std::copy(column.begin(), column.end(), tmp.begin());
                std::sort(tmp.begin(), tmp.end(), [](const p &a, const p &b) { return a.first < b.first; });

                std::vector<unsigned int> res_indices;
                std::vector<T> res_values;
                res_indices.reserve(tmp.size());
                res_values.reserve(tmp.size());

                for (const auto &entry: tmp) {
                    res_indices.push_back(entry.first);
                    res_values.push_back(entry.second);
                }

                res_cols.emplace_back(n, std::move(res_indices), std::move(res_values));
            }

            return Mat(n, std::move(res_cols));
        }

    public:
        std::size_t n = 0;
        std::vector<Vec<T>> cols;
    };

    template<typename Stream, typename T>
    Stream &operator<<(Stream &stream, const Mat<T> &m) {
        std::size_t n = m.n;
        stream << std::setprecision(PRECISION);
        for (std::size_t i = 0; i < n; i++) {
            stream << "[ col " << i << "]: ";
            stream << m.cols[i] << std::endl;
        }
        return stream;
    }

}

#endif //QR_DECOMPOSITION_MATH_HPP
