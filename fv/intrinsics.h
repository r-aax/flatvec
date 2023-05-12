#ifndef INTRINSICS_H
#define INTRINSICS_H

#include "stdafx.h"

#include "zmm.h"

namespace fv
{
    ZMM _mm512_load_ps(void const* mem_addr);

    void _mm512_store_ps(void* mem_addr,
                         ZMM a);

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b);
}

#endif
