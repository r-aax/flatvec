#include "stdafx.h"

#include "intrinsics.h"

namespace fv
{
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
}
