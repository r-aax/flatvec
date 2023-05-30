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
    /// <returns>Full mask.</returns>
    Mask Mask::full()
    {
        Mask k;

        k.fill();

        return k;
    }

    // Check mask.

    /// <summary>
    /// Check if mask is empty.
    /// </summary>
    /// <returns>
    /// true - if mask is empty,
    /// false - if mask is not empty.
    /// </returns>
    bool Mask::is_empty() const
    {
        for (int i = 0; i < bits; i++)
        {
            if (is_true(i))
            {
                return false;
            }
        }

        return true;
    }
}
