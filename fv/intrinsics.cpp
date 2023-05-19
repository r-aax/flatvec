#include "stdafx.h"

#include "intrinsics.h"

namespace fv
{
    // Memory access.

    ZMM _mm512_load_ps(void const* mem_addr)
    {
        ZMM r;
        float const* faddr = static_cast<float const*>(mem_addr);

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, faddr[i]);
        }

        return r;
    }

    void _mm512_store_ps(void* mem_addr,
                         ZMM a)
    {
        float* faddr = static_cast<float*>(mem_addr);

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            faddr[i] = a.get<float>(i);
        }
    }

    // Arithmetic operations with 2 arguments.

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return x + y; });
    }

    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return x - y; });
    }

    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return x * y; });
    }

    ZMM _mm512_div_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return x / y; });
    }

    ZMM _mm512_min_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return std::min(x, y); });
    }

    ZMM _mm512_max_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return std::max(x, y); });
    }

    // Compare operations.

    Mask _mm512_cmpeq_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x == y; });
    }

    Mask _mm512_cmple_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x <= y; });
    }

    Mask _mm512_cmplt_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x < y; });
    }

    Mask _mm512_cmpneq_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x != y; });
    }

    Mask _mm512_cmpnle_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !(x <= y); });
    }

    Mask _mm512_cmpnlt_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !(x < y); });
    }

    Mask _mm512_cmpord_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !std::isnan(x) && !std::isnan(y); });
    }

    Mask _mm512_cmpunord_ps_mask(ZMM a,
                                 ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return std::isnan(x) || std::isnan(y); });
    }

    // Blend operations.

    ZMM _mm512_mask_blend_ps(Mask k,
                             ZMM a,
                             ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, k.is_set(i) ? b.get<float>(i) : a.get<float>(i));
        }

        return r;
    }
}
