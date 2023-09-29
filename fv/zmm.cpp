#include "stdafx.h"

#include "zmm.h"

#include "global_stat.h"
#include "control_graph.h"

namespace fv
{
    // Constructors.

    /// <summary>
    /// Constructor.
    /// </summary>
    ZMM::ZMM()
    {
        id = CG.zmm_id();
        data = new int8_t[bytes];

        clear();

        CG.reglink("new z", id);
        GS.create_vector();
    }

    /// <summary>
    /// Copy constructor.
    /// </summary>
    /// <param name="z">Another ZMM register.</param>
    ZMM::ZMM(const ZMM& z)
    {
        id = CG.zmm_id();
        data = new int8_t[bytes];

        for (int i = 0; i < count<int32_t>(); ++i)
        {
            set<int32_t>(i, z.get<int32_t>(i));
        }

        GS.copy_vector();
        CG.reglink("copy z", z.get_id(), id);
    }

    /// <summary>
    /// Assign operator.
    /// </summary>
    /// <param name="z">Another ZMM register.</param>
    /// <returns>ZMM register.</returns>
    ZMM& ZMM::operator=(const ZMM& z)
    {
        if (this != &z)
        {
            // We do not need new identifier.
            // We do not need new data allocation.

            for (int i = 0; i < count<int32_t>(); ++i)
            {
                set<int32_t>(i, z.get<int32_t>(i));
            }

            GS.copy_vector();
            CG.reglink("rewrite z copy from " + std::to_string(z.get_id()),
                       z.get_id(), id);
        }

        return *this;
    }

    /// <summary>
    /// Move constructor.
    /// </summary>
    /// <param name="z">ZMM register.</param>
    ZMM::ZMM(ZMM&& z)
    {
        // We do not need new number for move constructor.
        // It is single use of z register.
        id = z.get_id();
        data = z.data;

        z.data = nullptr;

        // We cannon register id with "move z"
        // because it is source and destination at a time.
        // CG.reglink("move z", id);

        GS.move_vector();
    }

    /// <summary>
    /// Move assignment.
    /// </summary>
    /// <param name="z">ZMM register.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM& ZMM::operator=(ZMM&& z)
    {
        if (this != &z)
        {
            // Before this operation we had some register,
            // and now we rewrite it.
            CG.reglink("rewrite z move from " + std::to_string(z.get_id()),
                       id);

            // We do not need new number for move assignment.
            // It is single use of z register.
            id = z.get_id();

            delete[] data;
            data = z.data;
            z.data = nullptr;

            GS.move_vector();
        }

        return *this;
    }

    /// <summary>
    /// Destructor.
    /// </summary>
    ZMM::~ZMM()
    {
        delete[] data;
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
