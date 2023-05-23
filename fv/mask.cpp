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
    bool Mask::is_set(int i) const
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

    // Generate full mask.

    /// <summary>
    /// Clear mask.
    /// </summary>
    void Mask::clear()
    {
        for (int i = 0; i < bits; i++)
        {
            data[i] = false;
        }
    }

    /// <summary>
    /// Set full mask.
    /// </summary>
    void Mask::set_full()
    {
        for (int i = 0; i < bits; i++)
        {
            data[i] = true;
        }
    }

    /// <summary>
    /// Generate full mask.
    /// </summary>
    /// <returns>Full mask.</returns>
    Mask Mask::full()
    {
        Mask k;

        k.set_full();

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
            if (is_set(i))
            {
                return false;
            }
        }

        return true;
    }
}
