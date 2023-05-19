#include "stdafx.h"

#include "zmm.h"

namespace fv
{
    // Constructors.

    ZMM::ZMM()
    {
        clear();
    }

    ZMM::ZMM(const ZMM& z)
    {
        for (int i = 0; i < count<int32_t>(); i++)
        {
            set<int32_t>(i, z.get<int32_t>(i));
        }
    }

    // Clear.
    void ZMM::clear()
    {
        for (int i = 0; i < count<int32_t>(); i++)
        {
            set<int32_t>(i, 0);
        }
    }
}
