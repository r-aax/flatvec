#ifndef INTRINSICS_EXT_H
#define INTRINSICS_EXT_H

#include "stdafx.h"

#include "zmm.h"
#include "mask.h"

namespace fv
{
    // Compare operations.

    Mask<16> _mm512_cmpge_ps_mask(ZMM a,
                                  ZMM b);

    Mask<16> _mm512_cmpgt_ps_mask(ZMM a,
                                  ZMM b);
}

#endif
