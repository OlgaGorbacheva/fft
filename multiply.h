#ifndef MULTIPLY_H
#define MULTIPLY_H

#include <complex>
#include <vector>
#include <cmath>
#include <boost/thread.hpp>
#include <iostream>
#include <functional>

namespace my {
class pol_multiply;
}

class my::pol_multiply {
public:
    typedef std::complex<double> c_num;

private:
    std::vector<c_num> bit_reverse_copy(std::vector<c_num> const &a);
    int rev(int k, int n);
    void find_roots(int m, int n, bool mod, std::vector<c_num> &roots);

    void iterative_fft(std::vector<c_num> &a, bool m); // 1 - прямой, 0 - обратный
    void non_parallel_iterative_fft(std::vector<c_num> &a, bool m); // 1 - прямой, 0 - обратный
public:
    pol_multiply();

    std::vector<c_num> fft(std::vector<c_num> const &a, bool par = true);
    std::vector<c_num> rev_fft(std::vector<c_num> const &a,  bool par = true);

    std::vector<c_num> operator()(std::vector<c_num> const &a, std::vector<c_num> const &b); //собственно умножение
};

#endif // MULTIPLY_H
