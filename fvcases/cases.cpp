#include "stdafx.h"

#include "cases.h"

namespace fv
{
    void scase_add_f32(int n,
                       float* a,
                       float* b,
                       float* c)
    {
        for (int i = 0; i < n; i++)
        {
            c[i] = a[i] + b[i];
        }
    }

    void vcase_add_f32(int n,
                       float* a,
                       float* b,
                       float* c)
    {
        assert(n % ZMM::f32_count == 0);
        int vn = n / ZMM::f32_count;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::f32_count;

            ZMM va = _mm512_load_ps(a + sh);
            ZMM vb = _mm512_load_ps(b + sh);
            ZMM vc = _mm512_add_ps(va, vb);

            _mm512_store_ps(c + sh, vc);
        }
    }
}
