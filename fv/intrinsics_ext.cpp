#include "stdafx.h"

#include "intrinsics_ext.h"

#include "intrinsics.h"

namespace fv
{
    // Compare operations.

    /// <summary>
    /// Semantic for cmpge.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpge_ps_mask(ZMM& a,
                              ZMM& b)
    {
        return _mm512_cmpnlt_ps_mask(a, b);
    }

    /// <summary>
    /// Semantic for cmpge.
    /// </summary>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_mask_cmpge_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b)
    {
        return _mm512_mask_cmpnlt_ps_mask(k, a, b);
    }

    //

    /// <summary>
    /// Semantic for cmpgt.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpgt_ps_mask(ZMM& a,
                              ZMM& b)
    {
        return _mm512_cmpnle_ps_mask(a, b);
    }

    /// <summary>
    /// Semantic for cmpgt.
    /// </summary>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_mask_cmpgt_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b)
    {
        return _mm512_mask_cmpnle_ps_mask(k, a, b);
    }
}
