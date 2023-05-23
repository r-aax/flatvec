#include "stdafx.h"

#ifndef MASK_H
#define MASK_H

namespace fv
{
    /// <summary>
    /// Mask.
    /// </summary>
    class Mask
    {

    public:

        // Bits count.

        /// <summary>
        /// Bits count.
        /// </summary>
        static const int bits = 64;

    private:

        /// <summary>
        /// Data.
        /// </summary>
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
