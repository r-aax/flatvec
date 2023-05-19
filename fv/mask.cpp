#include "stdafx.h"

#include "mask.h"

namespace fv
{
    // Access to data.

    bool Mask::get(int i) const
    {
        return data[i];
    }

    bool Mask::is_set(int i) const
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

    // Check mask.

    bool Mask::is_empty() const
    {
        for (int i = 0; i < bits; i++)
        {
            if (is_set(i))
            {
                return false;
            }
        }

        return true;
    }
}
