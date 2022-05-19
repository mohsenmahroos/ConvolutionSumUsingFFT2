// Accepted solution of FOUR ARRAYS problem 
// May 2022 Long-time Competitive Programming Contest 1
// at CodeChef.com

#include "FFT.hpp"
#include <iostream>
#include <numeric>

int main() {
    using namespace std;
    using namespace FFT;
    using i_limits = numeric_limits<size_t>;
    using FFT_engine = FFT::engine<18,float>;
    using ivector = vector<size_t>;
    using cvector = FFT_engine::cvector;
    using cscalar = FFT_engine::cscalar;
    using ull = unsigned long long;
    const auto rounded_real = [&](const cscalar& z) { return round(real(z)); };
    const FFT_engine fft;
    const size_t M = 4, N = fft.N, u_max = i_limits::max(), u_min = i_limits::min();
    ivector array_size(M), f_min(M,u_max), f_max(M,u_min), g(N), h(N);
    vector<ivector> f(M,ivector(N));
    vector<cvector> F(M,cvector(N));
    ull K;
    const auto y_index = [&](ull m, ull x) {
        const ull y = N-1;
        return x == 0 ? y: min(y,m/x); };
    const auto p_count = [&](ull m) {
        ull n = 0;
        for (size_t x = 0; x < N and n < K; ++x)
            if (g[x] > 0)
                n += 1ll*g[x]*h[y_index(m,x)];
        return n; };
    const auto binary_search = [&](ull l, ull r) {
        ull ans = r;
        while (l <= r) {
            const ull m = (l+r)/2;
            if (p_count(m) < K)
                l = m+1;
            else
                r = m-1, ans = min(ans,m); }
        return ans; };
    const auto product_of_sums = [&](const ivector &x) {
        return 1ull*(x[0]+x[1])*(x[2]+x[3]); };
    cin.tie(nullptr)->sync_with_stdio(false);
    for (size_t i = 0; i < M; ++i)
        cin >> array_size[i];
    cin >> K;
    for (size_t i = 0; i < M; ++i)
        for (size_t elem; array_size[i] > 0; --array_size[i], ++f[i][elem])
            cin >> elem,
            f_min[i] = min(f_min[i],elem),
            f_max[i] = max(f_max[i],elem);
    for (size_t i = 0; i < M; ++i)
        for (size_t j = 0; j < N; ++j)
            F[i][j] = cscalar(f[i][j],0);
    const auto G = fft.convolution_sum(F[0],F[1]);
    const auto H = fft.convolution_sum(F[2],F[3]);
    for (size_t i = 0; i < N; ++i)
        g[i] = rounded_real(G[i]),
        h[i] = rounded_real(H[i]);
    for (size_t y = 1; y < N; ++y)
        h[y] += h[y-1];
    cout << binary_search(product_of_sums(f_min),product_of_sums(f_max)); }
