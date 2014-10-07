//#define BOOST_TEST_MODULE fft_test
#include "multiply.h"
#include <vector>
#include <iostream>
//#include <boost/test/included/unit_test.hpp>
#include <ctime>

using namespace std;
using namespace my;

const double eps = 1e-6;

int main() {
    pol_multiply mul;
    vector<pol_multiply::c_num> a = {4, 3, 2, 1, 0, 0, 0, 0, 4, 3, 2, 1, 0, 0, 0, 0};
//    vector<pol_multiply::c_num> b = {1, 1};
    vector<pol_multiply::c_num> c = mul.fft(a);//mul(a, b);
    for (unsigned int i = 0; i < c.size(); ++i) {
        cout << c[i] << endl;
    }
    return 0;
}

//vector<pol_multiply::c_num> plain_multiply(vector<pol_multiply::c_num> const &a, vector<pol_multiply::c_num> const &b) {
//    size_t n = 2 * max(a.size(), b.size());
//    vector<pol_multiply::c_num> res(n);
//    for (unsigned int i = 0; i < a.size(); ++i) {
//        for (unsigned int j = 0; j < b.size(); ++j) {
//            res[i + j] += a[i] * b[j];
//        }
//    }
//    return res;
//}

//void make_rand(vector<pol_multiply::c_num> &v) {
//    srand (time(NULL));
//    for (int i = 0; i < 100; i++) {
//        v.push_back(rand() % 100);
//    }
//}

//BOOST_AUTO_TEST_SUITE (fft_test)

//BOOST_AUTO_TEST_CASE (multiply) {
//    pol_multiply mul;
//    vector<pol_multiply::c_num> a, b, c, d;
//    make_rand(a);
//    make_rand(b);
//    c = plain_multiply(a, b);
//    cout << "i'm hear" << endl;
//    d = mul(a, b);
//    int n = max(c.size(), d.size());
//    cout << "i'm hear" << endl;
//    c.resize(n, 0);
//    d.resize(n, 0);
//    for (int i = 0; i < n; ++i) {
//        BOOST_CHECK_LT(fabs(c[i].real() - d[i].real()), eps);
//    }
//}

//BOOST_AUTO_TEST_CASE (parallel) {
//    pol_multiply mul;
//    vector<pol_multiply::c_num> a, fft, p_fft;
//    make_rand(a);
//    a.resize(128, 0);
//    clock_t par, non_par;
//    non_par = clock();
//    fft = mul.fft(a, 1);
//    cout << "i'm hear" << endl;
//    non_par = clock() - non_par;
//    par = clock();
//    p_fft = mul.fft(a, 0);
//    cout << "i'm hear" << endl;
//    par = clock() - par;
//    cout << (long double) (non_par - par) / CLOCKS_PER_SEC  << endl;
//    int n = fft.size();
//    for (int i = 0; i < n; ++i) {
//        BOOST_CHECK_LT(fabs(fft[i].real() - p_fft[i].real()), eps);
//        BOOST_CHECK_LT(fabs(fft[i].imag() - p_fft[i].imag()), eps);
//    }
//}

//BOOST_AUTO_TEST_SUITE_END( )
