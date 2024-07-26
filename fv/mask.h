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
        /// Use here the most possible length of masks.
        /// </summary>
        static const int bits = 64;

    private:

        /// <summary>
        /// Identifier.
        /// </summary>
        int id;

        /// <summary>
        /// Data.
        /// </summary>
        bool data[bits];

    public:

        // Constructors.

        Mask();

        Mask(uint64_t bm);

        Mask(const Mask& m);

        Mask& operator=(const Mask& m);

        Mask(Mask&& m);

        Mask& operator=(Mask&& m);

        // Access to data.

        /// <summary>
        /// Get identifier.
        /// </summary>
        /// <returns>Identifier.</returns>
        inline int get_id() const
        {
            return id;
        }

        bool get(int i) const;

        bool is_true(int i) const;

        void set(int i,
                 bool v);

        uint64_t binary_mask() const;

        operator uint64_t() const;

        /// <summary>
        /// Operator |.
        /// </summary>
        /// <param name="m">Second mask.</param>
        /// <returns>New mask.</returns>
        Mask operator|(const Mask& m) const
        {
            return Mask(binary_mask() | m.binary_mask());
        }

        std::string view() const;

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

        // Get information about mask density.

        /// <summary>
        /// Count of 1 bits in the tail.
        /// </summary>
        /// <param name="tail_len">Tail length.</param>
        /// <returns>Count of 1 bits count in the tail.</returns>
        int popcnt_tail(int tail_len);

        /// <summary>
        /// Mask density.
        /// </summary>
        /// <param name="tail_len">Length of tail.</param>
        /// <returns>Density.</returns>
        double density(int tail_len);
    };
}

#endif
