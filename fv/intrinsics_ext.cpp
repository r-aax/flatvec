#include "stdafx.h"

#include "intrinsics_ext.h"

#include "intrinsics.h"

namespace fv
{
    // Compare operations.

    Mask<16> _mm512_cmpge_ps_mask(ZMM a,
                                  ZMM b)
    {
        return _mm512_cmpnlt_ps_mask(a, b);
    }

    Mask<16> _mm512_cmpgt_ps_mask(ZMM a,
                                  ZMM b)
    {
        return _mm512_cmpnle_ps_mask(a, b);
    }
}
