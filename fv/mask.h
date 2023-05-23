#include "stdafx.h"

#ifndef MASK_H
#define MASK_H

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

        // Constructor.

        Mask();

        // Access to data.

        bool get(int i) const;

        bool is_set(int i) const;

        void set(int i,
                 bool v);

        // Generate full mask.

        void clear();

        void set_full();

        static Mask full();

        // Check mask.

        bool is_empty() const;
    };
}

#endif
