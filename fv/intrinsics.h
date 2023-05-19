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
