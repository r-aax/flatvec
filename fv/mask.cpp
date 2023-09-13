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

    /// <summary>
    /// Constructor from int.
    /// </summary>
    /// <param name="bm">Binary representation.</param>
    Mask::Mask(int bm)
        : Mask()
    {
        for (int i = 0; i < Mask::bits; i++)
        {
            data[i] = ((1 << i) != 0x0);
        }
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

    /// <summary>
    /// Mask in binary representation.
    /// </summary>
    /// <returns>Mask binary representation.</returns>
    int Mask::binary_mask() const
    {
        int bm = 0x0;

        for (int i = 0; i < Mask::bits; i++)
        {
            bm |= (data[i] << i);
        }

        return bm;
    }

    /// <summary>
    /// String view.
    /// </summary>
    /// <returns>String.</returns>
    std::string Mask::view() const
    {
        std::string s = "";

        for (int i = 0; i < Mask::bits; i++)
        {
            s = (data[i] ? "1" : "0") + s;
        }

        return "[" + s + "]";
    }

    /// <summary>
    /// Check equality.
    /// </summary>
    /// <param name="bm">Mask binary representation.</param>
    /// <returns>
    /// true - if mask matches given binary representation,
    /// false - otherwise.
    /// </returns>
    bool Mask::operator==(int bm) const
    {
        return binary_mask() == bm;
    }

    /// <summary>
    /// Check inequality.
    /// </summary>
    /// <param name="bm">Mask binary representation.</param>
    /// <returns>
    /// true - if mask does not match binary representation,
    /// false - otherwise.
    /// </returns>
    bool Mask::operator!=(int bm) const
    {
        return !(*this == bm);
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
    int Mask::popcnt_tail(int tail_len)
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
        double b1 = static_cast<double>(popcnt_tail(tail_len));
        double b = static_cast<double>(tail_len);

        return b1 / b;
    }
}
