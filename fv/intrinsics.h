#include "stdafx.h"

#ifndef INTRINSICS_H
#define INTRINSICS_H

#include "zmm.h"
#include "mask.h"
#include "global_stat.h"
#include "control_graph.h"

namespace fv
{
    /// <summary>
    /// Stub vector for vector operations.
    /// </summary>
    extern ZMM stub;

    /// <summary>
    /// Full mask for vector operations.
    /// </summary>
    extern Mask full;

    // Init operations.

    ZMM _mm512_set1_ps(float a);

    // Memory access.

    ZMM _mm512_load_ps(void const* mem_addr);

    void _mm512_store_ps(void* mem_addr,
                         ZMM& a);

    // Arithmetic operations with 1 argument. 

    /// <summary>
    /// Mask arithmetic operation with 1 argument.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">Argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result ZMM register.</returns>
    template <typename T>
    ZMM mask_arith1(ZMM& src,
                    Mask& k,
                    ZMM& a,
                    std::function<T(T)> op)
    {
        int w = ZMM::count<T>();
        ZMM dst;

        for (int i = 0; i < w; i++)
        {
            if (k.is_true(i))
            {
                dst.set<T>(i, op(a.get<T>(i)));
            }
            else
            {
                dst.set<T>(i, src.get<T>(i));
            }
        }

        GS.append_vector_oper(w, k.popcnt_tail(w));
        CG.reglink("arith1", src.get_id(), k.get_id(),
                   a.get_id(), dst.get_id());

        return dst;
    }

    /// <summary>
    /// Arithmetic operation with 1 argument.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="a">Argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result ZMM register.</returns>
    template <typename T>
    ZMM arith1(ZMM& a,
               std::function<T(T)> op)
    {
        return mask_arith1<T>(stub, full, a, op);
    }

    //

    ZMM _mm512_mask_mov_ps(ZMM& src,
                           Mask& k,
                           ZMM& a);

    //

    ZMM _mm512_abs_ps(ZMM& a);

    //

    ZMM _mm512_sqrt_ps(ZMM& a);

    ZMM _mm512_mask_sqrt_ps(ZMM& src,
                            Mask& k,
                            ZMM& a);

    // Arithmetic operations with 2 arguments.

    /// <summary>
    /// Mask arithmetic operation with 2 arguments.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="src">Source.</param>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="op">Operation.</param>
    /// <param name="is_z">Zero flag.</param>
    /// <returns>Result ZMM register.</returns>
    template <typename T>
    ZMM mask_arith2(ZMM& src,
                    Mask& k,
                    ZMM& a,
                    ZMM& b,
                    std::function<T(T, T)> op,
                    bool is_z = false)
    {
        int w = ZMM::count<T>();
        ZMM dst;

        for (int i = 0; i < w; i++)
        {
            if (k.is_true(i))
            {
                dst.set<T>(i, op(a.get<T>(i), b.get<T>(i)));
            }
            else
            {
                if (is_z)
                {
                    dst.set<T>(i, static_cast<T>(0));
                }
                else
                {
                    dst.set<T>(i, src.get<T>(i));
                }
            }
        }

        GS.append_vector_oper(w, k.popcnt_tail(w));
        CG.reglink("arith2", src.get_id(), k.get_id(),
                   a.get_id(), b.get_id(), dst.get_id());

        return dst;
    }

    /// <summary>
    /// Mask operation with 2 arguments.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="k">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result ZMM register.</returns>
    template <typename T>
    ZMM maskz_arith2(Mask& k,
                     ZMM& a,
                     ZMM& b,
                     std::function<T(T, T)> op)
    {
        return mask_arith2<T>(stub, k, a, b, op, true);
    }

    /// <summary>
    /// Arithmetic operation with 2 argument.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result ZMM register.</returns>
    template <typename T>
    ZMM arith2(ZMM& a,
               ZMM& b,
               std::function<T(T, T)> op)
    {
        return mask_arith2<T>(stub, full, a, b, op);
    }

    //

    ZMM _mm512_add_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_add_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_add_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    //

    ZMM _mm512_sub_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_sub_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_sub_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    //

    ZMM _mm512_mul_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_mul_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_mul_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    //

    ZMM _mm512_div_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_div_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_div_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    //

    ZMM _mm512_min_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_min_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_min_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    //

    ZMM _mm512_max_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_max_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_max_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    //

    ZMM _mm512_pow_ps(ZMM& a,
                      ZMM& b);

    ZMM _mm512_mask_pow_ps(ZMM& src,
                           Mask& k,
                           ZMM& a,
                           ZMM& b);

    ZMM _mm512_maskz_pow_ps(Mask& k,
                            ZMM& a,
                            ZMM& b);

    // Arithmetic operations with 3 arguments.

    /// <summary>
    /// Arithmetic operation with 3 arguments.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="c">Second argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result ZMM register.</returns>
    template <typename T>
    ZMM arith3(ZMM& a,
               ZMM& b,
               ZMM& c,
               std::function<T(T, T, T)> op)
    {
        int w = ZMM::count<T>();
        ZMM dst;

        for (int i = 0; i < w; i++)
        {
            dst.set<T>(i, op(a.get<T>(i), b.get<T>(i), c.get<T>(i)));
        }

        GS.append_vector_oper(w, w);
        CG.reglink("arith3", a.get_id(), b.get_id(), c.get_id(), dst.get_id());

        return dst;
    }

    //

    ZMM _mm512_fmadd_ps(ZMM& a,
                        ZMM& b,
                        ZMM& c);

    ZMM _mm512_fnmadd_ps(ZMM& a,
                         ZMM& b,
                         ZMM& c);

    //

    ZMM _mm512_fmsub_ps(ZMM& a,
                        ZMM& b,
                        ZMM& c);

    ZMM _mm512_fnmsub(ZMM& a,
                      ZMM& b,
                      ZMM& c);

    // Compare operations.

    /// <summary>
    /// Mask compare operation.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="k1">Mask.</param>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result mask.</returns>
    template <typename T>
    Mask mask_compare(Mask& k1,
                      ZMM& a,
                      ZMM& b,
                      std::function<bool(T, T)> op)
    {
        int w = ZMM::count<T>();
        Mask k;

        for (int i = 0; i < w; i++)
        {
            if (k1.is_true(i))
            {
                k.set(i, op(a.get<float>(i), b.get<float>(i)));
            }
            else
            {
                k.set(i, false);
            }
        }

        GS.append_vector_oper(w, k1.popcnt_tail(w));
        CG.reglink("cmp", k1.get_id(), a.get_id(), b.get_id(), k.get_id());

        return k;
    }

    /// <summary>
    /// Compare operation.
    /// </summary>
    /// <typeparam name="T">Type.</typeparam>
    /// <param name="a">First argument.</param>
    /// <param name="b">Second argument.</param>
    /// <param name="op">Operation.</param>
    /// <returns>Result mask.</returns>
    template <typename T>
    Mask compare(ZMM& a,
                 ZMM& b,
                 std::function<bool(T, T)> op)
    {
        return mask_compare<T>(full, a, b, op);
    }

    //

    Mask _mm512_cmpeq_ps_mask(ZMM& a,
                              ZMM& b);

    Mask _mm512_mask_cmpeq_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b);

    //

    Mask _mm512_cmple_ps_mask(ZMM& a,
                              ZMM& b);

    Mask _mm512_mask_cmple_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b);

    //

    Mask _mm512_cmplt_ps_mask(ZMM& a,
                              ZMM& b);

    Mask _mm512_mask_cmplt_ps_mask(Mask& k,
                                   ZMM& a,
                                   ZMM& b);

    //

    Mask _mm512_cmpneq_ps_mask(ZMM& a,
                               ZMM& b);

    Mask _mm512_mask_cmpneq_ps_mask(Mask& k,
                                    ZMM& a,
                                    ZMM& b);

    //

    Mask _mm512_cmpnle_ps_mask(ZMM& a,
                               ZMM& b);

    Mask _mm512_mask_cmpnle_ps_mask(Mask& k,
                                    ZMM& a,
                                    ZMM& b);

    //

    Mask _mm512_cmpnlt_ps_mask(ZMM& a,
                               ZMM& b);

    Mask _mm512_mask_cmpnlt_ps_mask(Mask& k,
                                    ZMM& a,
                                    ZMM& b);

    //

    Mask _mm512_cmpord_ps_mask(ZMM& a,
                               ZMM& b);

    //

    Mask _mm512_cmpunord_ps_mask(ZMM& a,
                                 ZMM& b);

    // Blend operations.

    ZMM _mm512_mask_blend_ps(Mask& k,
                             ZMM& a,
                             ZMM& b);

    // Operations with masks.

    Mask koperation(Mask& a,
                    Mask& b,
                    std::function<bool(bool, bool)> op);

    Mask _mm512_kmov(Mask& a);

    Mask _mm512_knot(Mask& a);

    Mask _mm512_kand(Mask& a,
                     Mask& b);

    Mask _mm512_kandn(Mask& a,
                      Mask& b);

    Mask _mm512_kor(Mask& a,
                    Mask& b);

    Mask _mm512_kxor(Mask& a,
                     Mask& b);

    Mask _mm512_kxnor(Mask& a,
                      Mask& b);
}

#endif
