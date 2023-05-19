#include "stdafx.h"

#include "mask.h"

namespace fv
{
    // Access to data.

    bool Mask::get(int i)
    {
        return data[i];
    }

    bool Mask::is_set(int i)
    {
        return get(i);
    }

    void Mask::set(int i,
                   bool v)
    {
        data[i] = v;
    }

    // Generate full mask.

    void Mask::set_full()
    {
        for (int i = 0; i < bits; i++)
        {
            data[i] = true;
        }
    }

    Mask Mask::full()
    {
        Mask k;

        k.set_full();

        return k;
    }
}
