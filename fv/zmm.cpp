#include "stdafx.h"

#include "zmm.h"

namespace fv
{
    // Constructors.

    /// <summary>
    /// Constructor.
    /// </summary>
    ZMM::ZMM()
    {
        clear();

        GS.create_vector();
    }

    /// <summary>
    /// Copy constructor.
    /// </summary>
    /// <param name="z">Another ZMM register.</param>
    ZMM::ZMM(const ZMM& z)
    {
        for (int i = 0; i < count<int32_t>(); i++)
        {
            set<int32_t>(i, z.get<int32_t>(i));
        }

        GS.copy_vector();
    }

    /// <summary>
    /// Clear.
    /// </summary>
    void ZMM::clear()
    {
        for (int i = 0; i < count<int32_t>(); i++)
        {
            set<int32_t>(i, 0);
        }
    }
}
