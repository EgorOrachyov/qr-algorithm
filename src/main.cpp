#include "math.hpp"
#include "qr_omp.hpp"
#include "qr_seq.hpp"

int main(int argc, const char* const *argv) {
    std::size_t n = 4;
    auto m = qr::Mat<double>::generate_symmetric(n);

    std::cout << m;

    return 0;
}