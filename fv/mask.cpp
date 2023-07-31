#include "stdafx.h"

#include "mask.h"

namespace fv
{
    // Constructor.

    /// <summary>
    /// Default constructor.
    /// </summary>
    Mask::Mask()
    {
        clear();
    }

    // Access to data.

    /// <summary>
    /// Get mask element.
    /// </summary>
    /// <param name="i">Index.</param>
    /// <returns>Mask element.</returns>
    bool Mask::get(int i) const
    {
        return data[i];
    }

    /// <summary>
    /// Check if mask element is set.
    /// </summary>
    /// <param name="i">Index.</param>
    /// <returns>Mask element value.</returns>
    bool Mask::is_true(int i) const
    {
        return get(i);
    }

    /// <summary>
    /// Set mask element.
    /// </summary>
    /// <param name="i">Index.</param>
    /// <param name="v">Value.</param>
    void Mask::set(int i,
                   bool v)
    {
        data[i] = v;
    }

    // Generate masks.

    /// <summary>
    /// Set tail.
    /// </summary>
    /// <param name="n">Tail size.</param>
    /// <param name="v">Value.</param>
    void Mask::set_tail(int n,
                        bool v)
    {
        for (int i = 0; i < n; i++)
        {
            data[i] = v;
        }
    }

    /// <summary>
    /// Generate full mask.
    /// </summary>
    /// <param name="n">Tail size.</param>
    /// <returns>Full mask.</returns>
    Mask Mask::full_tail(int n)
    {
        Mask k;

        k.fill_tail(n);

        return k;
    }

    // Check mask.

    /// <summary>
    /// Check if mask is empty.
    /// </summary>
    /// <param name="n">Tail size.</param>
    /// <returns>
    /// true - if mask is empty,
    /// false - if mask is not empty.
    /// </returns>
    bool Mask::is_empty_tail(int n) const
    {
        for (int i = 0; i < n; i++)
        {
            if (is_true(i))
            {
                return false;
            }
        }

        return true;
    }

    // Get information about mask density.

    /// <summary>
    /// Count of 1 bits in the tail.
    /// </summary>
    /// <param name="tail_len">Tail length.</param>
    /// <returns>Count of 1 bits count in the tail.</returns>
    int Mask::tail_1_bits_count(int tail_len)
    {
        int cnt = 0;

        for (int i = 0; i < tail_len; i++)
        {
            if (data[i])
            {
                cnt++;
            }
        }

        return cnt;
    }

    /// <summary>
    /// Mask density.
    /// </summary>
    /// <param name="tail_len">Length of tail.</param>
    /// <returns>Density.</returns>
    double Mask::density(int tail_len)
    {
        double b1 = static_cast<double>(tail_1_bits_count(tail_len));
        double b = static_cast<double>(tail_len);

        return b1 / b;
    }
}
