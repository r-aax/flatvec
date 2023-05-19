#include "stdafx.h"

#include "cases.h"

namespace fv
{
    void scase_arith_f32(int n,
                         float* a,
                         float* b,
                         float* c)
    {
        for (int i = 0; i < n; i++)
        {
            c[i] = (a[i] + b[i]) + (a[i] - b[i])
                + (a[i] * b[i]) + (a[i] / b[i])
                + std::min(a[i], b[i]) + std::max(a[i], b[i]);
        }
    }

    void vcase_arith_f32(int n,
                         float* a,
                         float* b,
                         float* c)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM va = _mm512_load_ps(a + sh);
            ZMM vb = _mm512_load_ps(b + sh);

            ZMM vadd = _mm512_add_ps(va, vb);
            ZMM vsub = _mm512_sub_ps(va, vb);
            ZMM vmul = _mm512_mul_ps(va, vb);
            ZMM vdiv = _mm512_div_ps(va, vb);
            ZMM vmin = _mm512_min_ps(va, vb);
            ZMM vmax = _mm512_max_ps(va, vb);

            ZMM vc = _mm512_add_ps(vadd, vsub);
            vc = _mm512_add_ps(vc, vmul);
            vc = _mm512_add_ps(vc, vdiv);
            vc = _mm512_add_ps(vc, vmin);
            vc = _mm512_add_ps(vc, vmax);

            _mm512_store_ps(c + sh, vc);
        }
    }
}
