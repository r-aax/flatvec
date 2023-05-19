#include "stdafx.h"

#include "intrinsics.h"

namespace fv
{
    ZMM _mm512_load_ps(void const* mem_addr)
    {
        ZMM r;
        float const* faddr = static_cast<float const*>(mem_addr);

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, faddr[i]);
        }

        return r;
    }

    void _mm512_store_ps(void* mem_addr,
                         ZMM a)
    {
        float* faddr = static_cast<float*>(mem_addr);

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            faddr[i] = a.getF32(i);
        }
    }

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, a.getF32(i) + b.getF32(i));
        }

        return r;
    }

    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, a.getF32(i) - b.getF32(i));
        }

        return r;
    }

    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, a.getF32(i) * b.getF32(i));
        }

        return r;
    }

    ZMM _mm512_div_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, a.getF32(i) / b.getF32(i));
        }

        return r;
    }

    ZMM _mm512_min_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, std::min(a.getF32(i), b.getF32(i)));
        }

        return r;
    }

    ZMM _mm512_max_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::f32_count; i++)
        {
            r.setF32(i, std::max(a.getF32(i), b.getF32(i)));
        }

        return r;
    }
}
