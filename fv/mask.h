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

        static Mask full_tail(int n);

        /// <summary>
        /// Mask with full tail.
        /// </summary>
        /// <returns>Mask with full tail.</returns>
        inline static Mask full()
        {
            return full_tail(bits);
        }

        // Check mask.

        bool is_empty_tail(int n) const;

        /// <summary>
        /// Check if tail is empty.
        /// </summary>
        /// <returns>
        /// true - if mask is empty,
        /// false - if mask is not empty.
        /// </returns>
        inline bool is_empty() const
        {
            return is_empty_tail(bits);
        }
    };
}

#endif
