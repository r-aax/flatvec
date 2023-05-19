#ifndef MASK_H
#define MASK_H

#include "stdafx.h"

namespace fv
{
    // Mask.
    class Mask
    {

    public:

        // Bits count.

        static const int bits = 64;

    private:

        // Mask data.
        bool data[bits];

    public:

        // Access to data.

        bool get(int i) const;

        bool is_set(int i) const;

        void set(int i,
                 bool v);

        // Generate full mask.

        void set_full();

        static Mask full();

        // Check mask.

        bool is_empty() const;
    };
}

#endif
