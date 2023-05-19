#ifndef INTRINSICS_H
#define INTRINSICS_H

#include "stdafx.h"

#include "zmm.h"

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
               std::function<T(T, T)> op2)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<T>(); i++)
        {
            r.set<T>(i, op2(a.get<T>(i), b.get<T>(i)));
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
}

#endif
