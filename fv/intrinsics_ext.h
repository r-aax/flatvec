#include "stdafx.h"

#ifndef INTRINSICS_EXT_H
#define INTRINSICS_EXT_H

#include "zmm.h"
#include "mask.h"

namespace fv
{
    // Compare operations.

    Mask _mm512_cmpge_ps_mask(ZMM& a,
                              ZMM& b);

    Mask _mm512_mask_cmpge_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b);

    //

    Mask _mm512_cmpgt_ps_mask(ZMM& a,
                              ZMM& b);

    Mask _mm512_mask_cmpgt_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b);
}

#endif
