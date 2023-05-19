#ifndef INTRINSICS_H
#define INTRINSICS_H

#include "stdafx.h"

#include "zmm.h"
#include "mask.h"

namespace fv
{
    // Memory access.

    ZMM _mm512_load_ps(void const* mem_addr);

    void _mm512_store_ps(void* mem_addr,
                         ZMM a);

    // Arithmetic operations with 2 arguments.

    template <typename T>
    ZMM arith2(ZMM a,
               ZMM b,
               std::function<T(T, T)> op)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<T>(); i++)
        {
            r.set<T>(i, op(a.get<T>(i), b.get<T>(i)));
        }

        return r;
    }

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_div_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_min_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_max_ps(ZMM a,
                      ZMM b);

    // Compare operations.

    template <typename T>
    Mask<16> compare(ZMM a,
                     ZMM b,
                     std::function<bool(T, T)> op)
    {
        Mask<16> m;

        for (int i = 0; i < ZMM::count<T>(); i++)
        {
            m.set(i, op(a.get<T>(i), b.get<T>(i)));
        }

        return m;
    }

    Mask<16> _mm512_cmpeq_ps_mask(ZMM a,
                                  ZMM b);

    Mask<16> _mm512_cmple_ps_mask(ZMM a,
                                  ZMM b);

    Mask<16> _mm512_cmplt_ps_mask(ZMM a,
                                  ZMM b);

    Mask<16> _mm512_cmpneq_ps_mask(ZMM a,
                                   ZMM b);

    Mask<16> _mm512_cmpnle_ps_mask(ZMM a,
                                   ZMM b);

    Mask<16> _mm512_cmpnlt_ps_mask(ZMM a,
                                   ZMM b);

    Mask<16> _mm512_cmpord_ps_mask(ZMM a,
                                   ZMM b);

    Mask<16> _mm512_cmpunord_ps_mask(ZMM a,
                                     ZMM b);

    // Blend operations.

    ZMM _mm512_mask_blend_ps(Mask<16> m,
                             ZMM a,
                             ZMM b);
}

#endif
