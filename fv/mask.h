#ifndef MASK_H
#define MASK_H

#include "stdafx.h"

namespace fv
{
    // Mask.
    template <int Size>
    class Mask
    {

    private:

        // Mask data.
        bool data[Size];

    public:

        // Access to data.

        bool get(int i)
        {
            return data[i];
        }

        bool is_set(int i)
        {
            return get(i);
        }

        void set(int i,
                 bool v)
        {
            data[i] = v;
        }
    };
}

#endif
