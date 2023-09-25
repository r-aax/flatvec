#include "stdafx.h"

#include "mask.h"

#include "global_stat.h"
#include "control_graph.h"

namespace fv
{
    // Constructor.

    /// <summary>
    /// Default constructor.
    /// </summary>
    Mask::Mask()
    {
        id = CG.mask_id();

        clear();

        CG.reglink("new m", id);
    }

    /// <summary>
    /// Constructor from int.
    /// </summary>
    /// <param name="bm">Binary representation.</param>
    Mask::Mask(uint64_t bm)
    {
        id = CG.mask_id();

        for (int i = 0; i < Mask::bits; ++i)
        {
            data[i] = ((bm & (static_cast<uint64_t>(1) << i)) != 0x0);
        }

        CG.reglink("new m from const", id);
    }

    /// <summary>
    /// Copy constructor.
    /// </summary>
    /// <param name="m">Another mask.</param>
    Mask::Mask(const Mask& m)
    {
        id = m.get_id();

        for (int i = 0; i < Mask::bits; ++i)
        {
            data[i] = m.data[i];
        }

        CG.reglink("copy m", m.get_id(), id);
    }

    /// <summary>
    /// Assigmnent operator.
    /// </summary>
    /// <param name="m">Another mask.</param>
    /// <returns>Result mask.</returns>
    Mask& Mask::operator=(const Mask& m)
    {
        if (this != &m)
        {
            for (int i = 0; i < Mask::bits; ++i)
            {
                data[i] = m.data[i];
            }

            CG.reglink("rewrite m copy from " + std::to_string(m.get_id()), m.get_id(), id);
        }

        return *this;
    }

    /// <summary>
    /// Move constructor.
    /// </summary>
    /// <param name="m">Source mask.</param>
    Mask::Mask(Mask&& m)
    {
        // We do not need new identifier.
        id = m.get_id();

        for (int i = 0; i < Mask::bits; ++i)
        {
            data[i] = m.data[i];
        }
    }

    /// <summary>
    /// Move assignment.
    /// </summary>
    /// <param name="m">Another mask.</param>
    /// <returns>Result mask.</returns>
    Mask& Mask::operator=(Mask&& m)
    {
        if (this != &m)
        {
            // Before this operation we have some mask and now we rewrite it.
            CG.reglink("rewrite m move from " + std::to_string(m.get_id()), id);

            // We do not need new number for move assihnment.
            id = m.get_id();

            for (int i = 0; i < Mask::bits; ++i)
            {
                data[i] = m.data[i];
            }
        }

        return *this;
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
    uint64_t Mask::binary_mask() const
    {
        uint64_t bm = 0x0;

        for (int i = 0; i < Mask::bits; i++)
        {
            bm |= (static_cast<uint64_t>(data[i]) << i);
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
    bool Mask::operator==(uint64_t bm) const
    {
        // Maybe it is the last use.
        CG.reglink("control", get_id());

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
    bool Mask::operator!=(uint64_t bm) const
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
