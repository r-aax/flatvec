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

        bool is_true(int i) const;

        void set(int i,
                 bool v);

        // Generate masks.

        void set_tail(int n,
                      bool v);

        /// <summary>
        /// Clear tail.
        /// </summary>
        /// <param name="n">Tail size.</param>
        inline void clear_tail(int n)
        {
            set_tail(n, false);
        }

        /// <summary>
        /// Clear all bits.
        /// </summary>
        inline void clear()
        {
            clear_tail(bits);
        }

        /// <summary>
        /// Fill tail.
        /// </summary>
        /// <param name="n">Tail size.</param>
        inline void fill_tail(int n)
        {
            set_tail(n, true);
        }

        /// <summary>
        /// Fill all bits.
        /// </summary>
        void fill()
        {
            fill_tail(bits);
        }

        static Mask full();

        // Check mask.

        bool is_empty() const;
    };
}

#endif
