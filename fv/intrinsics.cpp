#include "stdafx.h"

#include "intrinsics.h"

namespace fv
{
    // Init operations.

    /// <summary>
    /// Semantic for set1.
    /// </summary>
    /// <param name="a">Value.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_set1_ps(float a)
    {
        ZMM dst;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            dst.set<float>(i, a);
        }

        return dst;
    }

    // Memory access.

    /// <summary>
    /// Semantic for load.
    /// </summary>
    /// <param name="mem_addr">Address.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_load_ps(void const* mem_addr)
    {
        ZMM r;
        float const* faddr = static_cast<float const*>(mem_addr);

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, faddr[i]);
        }

        return r;
    }

    /// <summary>
    /// Semantic for store.
    /// </summary>
    /// <param name="mem_addr">Address.</param>
    /// <param name="a">Stored register.</param>
    void _mm512_store_ps(void* mem_addr,
                         ZMM a)
    {
        float* faddr = static_cast<float*>(mem_addr);

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            faddr[i] = a.get<float>(i);
        }
    }

    // Arithmetic operations with 1 argument.

    /// <summary>
    /// Semantic for mov.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">Argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_mov_ps(ZMM src,
                           Mask k,
                           ZMM a)
    {
        return mask_arith1<float>(src, k, a, [] (float x) { return x; });
    }
    
    /// <summary>
    /// Semantic for sqrt.
    /// </summary>
    /// <param name="a">Argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_sqrt_ps(ZMM a)
    {
        return arith1<float>(a, [] (float x) { return sqrt(x); });
    }

    // Arithmetic operations with 2 arguments.

    /// <summary>
    /// Semantic for add.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_add_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::plus<float>());
    }

    /// <summary>
    /// Semantic for add.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_add_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::plus<float>());
    }

    /// <summary>
    /// Semantic for add.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_add_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::plus<float>());
    }

    //

    /// <summary>
    /// Semantic for sub.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::minus<float>());
    }

    /// <summary>
    /// Semantic for sub.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_sub_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::minus<float>());
    }

    /// <summary>
    /// Semantic for sub.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_sub_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::minus<float>());
    }

    //

    /// <summary>
    /// Semantic for mul.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::multiplies<float>());
    }

    /// <summary>
    /// Semantic for mul.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_mul_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::multiplies<float>());
    }

    /// <summary>
    /// Semantic for mul.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_mul_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::multiplies<float>());
    }

    //

    /// <summary>
    /// Semantic for div.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_div_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::divides<float>());
    }

    /// <summary>
    /// Semantic for div.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_div_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::divides<float>());
    }

    /// <summary>
    /// Semantic for div.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_div_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::divides<float>());
    }

    //

    /// <summary>
    /// Semantic for min.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_min_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return std::min(x, y); });
    }

    /// <summary>
    /// Semantic for min.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_min_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, [] (float x, float y) { return std::min(x, y); });
    }

    /// <summary>
    /// Semantic for min.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_min_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, [] (float x, float y) { return std::min(x, y); });
    }

    //

    /// <summary>
    /// Semantic for max.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_max_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return std::max(x, y); });
    }

    /// <summary>
    /// Semantic for max.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_max_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, [] (float x, float y) { return std::max(x, y); });
    }

    /// <summary>
    /// Semantic for max.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_max_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, [] (float x, float y) { return std::max(x, y); });
    }

    //

    /// <summary>
    /// Semantic for pow.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_pow_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return pow(x, y); });
    }

    /// <summary>
    /// Semantic for pow.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_mask_pow_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, [] (float x, float y) { return pow(x, y); });
    }

    /// <summary>
    /// Semantic for pow.
    /// </summary>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_maskz_pow_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, [] (float x, float y) { return pow(x, y); });
    }

    // Arithmetic operations with 3 arguments.

    /// <summary>
    /// Semantic for fmadd.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="c">Third argument.</param>
    /// <returns>Result ZMM register.</returns>
    ZMM _mm512_fmadd_ps(ZMM a,
                        ZMM b,
                        ZMM c)
    {
        ZMM dst;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            dst.set<float>(i, a.get<float>(i) * b.get<float>(i) + c.get<float>(i));
        }

        return dst;
    }

    // Compare operations.

    /// <summary>
    /// Semantic for cmpeq.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpeq_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x == y; });
    }

    //

    /// <summary>
    /// Semantic for cmple.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmple_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x <= y; });
    }

    //

    /// <summary>
    /// Semantic for cmplt.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmplt_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x < y; });
    }

    /// <summary>
    /// Semantic for cmplt.
    /// </summary>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_mask_cmplt_ps_mask(Mask k,
                                   ZMM a,
                                   ZMM b)
    {
        return mask_compare<float>(k, a, b, [] (float x, float y) { return x < y; });
    }

    //

    /// <summary>
    /// Semantic for cmpneq.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpneq_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x != y; });
    }

    //

    /// <summary>
    /// Semantic for cmple.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpnle_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !(x <= y); });
    }

    //

    /// <summary>
    /// Semantic for cmpnlt.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpnlt_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !(x < y); });
    }

    //

    /// <summary>
    /// Semantic for cmpord.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpord_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !std::isnan(x) && !std::isnan(y); });
    }

    //

    /// <summary>
    /// Semantic for cmpunord.
    /// </summary>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_cmpunord_ps_mask(ZMM a,
                                 ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return std::isnan(x) || std::isnan(y); });
    }

    // Blend operations.

    /// <summary>
    /// Semantic for blend.
    /// </summary>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <returns>Result mask.</returns>
    ZMM _mm512_mask_blend_ps(Mask k,
                             ZMM a,
                             ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, k.is_set(i) ? b.get<float>(i) : a.get<float>(i));
        }

        return r;
    }

    // Operations with masks.

    /// <summary>
    /// Logical operation.
    /// </summary>
    /// <param name="count">Count of bits.</param>
    /// <param name="a">First mask.</param>
    /// <param name="b">Second mask.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result mask.</returns>
    Mask koperation(int count,
                    Mask a,
                    Mask b,
                    std::function<bool(bool, bool)> op)
    {
        Mask k;

        for (int i = 0; i < Mask::bits; i++)
        {
            if (i < count)
            {
                k.set(i, op(a.get(i), b.get(i)));
            }
            else
            {
                k.set(i, false);
            }
        }

        return k;
    }

    /// <summary>
    /// Semantic for kmov.
    /// </summary>
    /// <param name="a">Mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_kmov(Mask a)
    {
        return koperation(ZMM::count<float>(), a, a, [] (bool x, bool y) { return x; });
    }

    /// <summary>
    /// Semantic for knot.
    /// </summary>
    /// <param name="a">Mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_knot(Mask a)
    {
        return koperation(ZMM::count<float>(), a, a, [] (bool x, bool y) { return !x; });
    }

    /// <summary>
    /// Semantic for kand.
    /// </summary>
    /// <param name="a">First mask.</param>
    /// <param name="b">Second mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_kand(Mask a,
                     Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, std::bit_and<bool>());
    }

    /// <summary>
    /// Semantic for kandn.
    /// </summary>
    /// <param name="a">First mask.</param>
    /// <param name="b">Second mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_kandn(Mask a,
                      Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, [] (bool x, bool y) { return !x & y; });
    }

    /// <summary>
    /// Semantic for kor.
    /// </summary>
    /// <param name="a">First mask.</param>
    /// <param name="b">Second mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_kor(Mask a,
                    Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, std::bit_or<bool>());
    }

    /// <summary>
    /// Semantic for kxor.
    /// </summary>
    /// <param name="a">First mask.</param>
    /// <param name="b">Second mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_kxor(Mask a,
                     Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, std::bit_xor<bool>());
    }

    /// <summary>
    /// Semantic for kxnor.
    /// </summary>
    /// <param name="a">First mask.</param>
    /// <param name="b">Second mask.</param>
    /// <returns>Result mask.</returns>
    Mask _mm512_kxnor(Mask a,
                      Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, [] (bool x, bool y) { return !(x ^ y); });
    }
}
