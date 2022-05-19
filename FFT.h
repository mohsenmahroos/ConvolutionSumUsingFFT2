#ifndef __FFT__
#define __FFT__
#include <cmath>
#include <complex>
#include <vector>

namespace FFT {
    template<size_t size, class number_t = float>
    struct FFT_engine {
        enum type {FORWARD,INVERSE};
        using cscalar = std::complex<number_t>;
        using cvector = std::vector<cscalar>;
        using cmatrix = std::vector<cvector>;
        using ctensor = std::vector<cmatrix>;
        const size_t N = 1<<size;
    private:
        struct unity_roots: ctensor {
            inline unity_roots() : ctensor(2,cmatrix(size)) {
                number_t phi = 4.0*atan(1);
                for (size_t s = 0, k = 1; s < size; ++s, k <<= 1, phi *= 0.5)
                    for (size_t j = 0; j < k; ++j) {
                        number_t theta = j*phi, x = cos(theta), y = sin(theta);
                        for (size_t u = 0; u < 2; ++u, y = -y)
                            this->at(u).at(s).emplace_back(x,y); } } };
        const unity_roots WP;
    public:
        inline FFT_engine() : WP() {}
        inline void transform(cvector &F, type t) const {
            for (size_t l = N>>1, m = N-1, j = 0, i = 1; i < m; ++i) {
                size_t k = l;
                while (j ^= k, j < k)
                    k >>= 1;
                if (i < j)
                    swap(F[i],F[j]); }
            for (size_t u = 0, k = 1, w = 2; u < size; ++u, k = w, w <<= 1)
                for (size_t i = 0; i < N; i += w)
                    for (size_t v = 0; v < k; ++v) {
                        const auto l = i+v, r = l+k;
                        const auto a = F[l], b = F[r]*WP[t][u][v];
                        F[l] = a+b, F[r] = a-b; }
            if (t == INVERSE)
                for (auto &value: F)
                    value /= N; }
        inline cvector convolution_sum(cvector F, cvector G) const {
            cvector H(N);
            transform(F,FORWARD), transform(G,FORWARD);
            for (size_t i = 0; i < N; ++i)
                H[i] = F[i]*G[i];
            transform(H,INVERSE);
            return H; } }; }
#endif // __FFT__
