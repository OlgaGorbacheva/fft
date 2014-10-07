#ifndef MULTIPLY_CPP
#define MULTIPLY_CPP

#include "multiply.h"
#include "threadpool.h"
#include <iostream>
#include <future>

#define PI 3.14159265358979323846

using namespace std;
using namespace my;

pol_multiply::pol_multiply() {
}
#include <iostream>
vector<pol_multiply::c_num> pol_multiply::fft(vector<pol_multiply::c_num> const &a, bool par) {
    if(a.size() == 0) {
        return vector<c_num>();
    }
    vector<c_num> answer = bit_reverse_copy(a);
    if (par)
        iterative_fft(answer, 1); // 1 - прямой, 0 - обратный
    else non_parallel_iterative_fft(answer, 1);
    return answer;
}

vector<pol_multiply::c_num> pol_multiply::rev_fft(vector<pol_multiply::c_num> const &a, bool par) {
    if(a.size() == 0) {
        return vector<c_num>();
    }
    vector<c_num> answer = bit_reverse_copy(a);
    if (par)
        iterative_fft(answer, 0); // 1 - прямой, 0 - обратный
    else non_parallel_iterative_fft(answer, 0);
    for (auto i = answer.begin(); i != answer.end(); ++i) {
        (*i) /= answer.size();
    }
    return answer;
}


int pol_multiply::rev(int k, int n) {
    int rev = 0;
    n = ceil(log(n)/log(2));
    for (int i = 0; i < n; i++) {
        int mask = 1 << i;
        mask = mask & k; //нашли i-й бит с конца
        if (i < n / 2)
            mask = mask << (n - 2 * i - 1); //сдвинули на место
        else {
            mask = mask >> i;
            mask = mask << (n - i - 1);
        }
        rev = rev | mask;
    }
    return rev;
}

vector<pol_multiply::c_num> pol_multiply::bit_reverse_copy(vector<c_num> const &_a) {
    vector<c_num> a(_a.size());
    int n = a.size();
    for (int i = 0; i < n; ++i) {
        a[rev(i, n)] = _a[i];
    }
    return a;
}

void pol_multiply::find_roots(int m, int n, bool mod, vector<c_num> &roots) {
    double alpha = 2 * ((mod) ? 1. : -1.) * PI / m;
    const c_num angle (cos(alpha), sin(alpha));
    roots.resize(n);
    roots[0] = 1;
    for (int i = 1; i <= m / 2; ++i) {
        roots[i] = roots[i - 1] * angle;
    }
}

#include <iostream>

void pol_multiply::iterative_fft(vector<pol_multiply::c_num> &a, bool mod) {
    int n = a.size();
    double lg_n = ceil(log(n) / log(2));
    threadpool pool;
    cout << "Размер массива " << n << ", логарифм от него " << lg_n << endl;
    for (int s = 1; s < lg_n + 1; ++s) {
        int m = 1 << s;
        cout << "Внешний цикл, s = " << s << ", m = " << m << endl;
        vector<c_num> roots;
        find_roots(m, n, mod, roots);
        std::function<void(int, int)> f = [&](int k1, int k2) {
            for (int k = k1; k < k2; k += m) {
                for (int j = 0; j < m / 2; ++j) {
                    c_num t = roots[j] * a[k + j + m / 2];
                    c_num u = a[k + j];
                    a[k + j] = u + t;
                    a[k + j + m / 2] = u - t;
                }
            }
        };
        int j = max(n / 8, m);
        cout << "Начинаю херачить в пул, размер участков " << j << endl;
        for (int i = 0; i < n; i += j) {
            pool.add<void>(f, i, min(i + j, n));
            cout << i << endl;
        }
        cout << "Все команды в пуле, ждем-с" << endl;
        pool.wait();
        cout << "О чудо! О великий бог многопоточности!" << endl;
    }
}

void pol_multiply::non_parallel_iterative_fft(vector<pol_multiply::c_num> &a, bool mod) {
    int n = a.size();
    double lg_n = ceil(log(n) / log(2));
    for (int s = 1; s < lg_n + 1; ++s) {
        int m = 1 << s;
        vector<c_num> roots;
        find_roots(m, n, mod, roots);
        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < m / 2; ++j) {
                c_num t = roots[j] * a[k + j + m / 2];
                c_num u = a[k + j];
                a[k + j] = u + t;
                a[k + j + m / 2] = u - t;
            }
        }
    }
}

vector<pol_multiply::c_num> pol_multiply::operator()(vector<pol_multiply::c_num> const &_a, vector<pol_multiply::c_num> const &_b) {
    vector<c_num> a(_a), b(_b);
    int n = max(a.size(), b.size());
    int n_log = ceil(log(n) / log(2));
    n = 1 << n_log;
    a.resize(2 * n, 0);
    b.resize(2 * n, 0);
    vector<c_num> a_vals = fft(a);
    vector<c_num> b_vals = fft(b);
    int m = a_vals.size();
    vector<c_num> c_vals(m);

    for (int i = 0; i < m; ++i)
      c_vals[i] = a_vals[i] * b_vals[i];

    vector<c_num> c = rev_fft(c_vals, 0);
    return c;
}

#endif //MULTIPLY_CPP
