#include "stdafx.h"

#ifndef INTRINSICS_H
#define INTRINSICS_H

#include "zmm.h"
#include "mask.h"

namespace fv
{
    // Memory access.

    ZMM _mm512_load_ps(void const* mem_addr);

    void _mm512_store_ps(void* mem_addr,
                         ZMM a);

    // Arithmetic operations with 1 argument.

    ZMM _mm512_mask_mov_ps(ZMM src,
                           Mask k,
                           ZMM a);

    ZMM _mm512_sqrt_ps(ZMM a);

    // Arithmetic operations with 2 arguments.

    template <typename T>
    ZMM mask_arith2(ZMM src,
                    Mask k,
                    ZMM a,
                    ZMM b,
                    std::function<T(T, T)> op,
                    bool is_z = false)
    {
        ZMM dst;

        for (int i = 0; i < ZMM::count<T>(); i++)
        {
            if (k.is_set(i))
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

        return dst;
    }

    template <typename T>
    ZMM maskz_arith2(ZMM src,
                     Mask k,
                     ZMM a,
                     ZMM b,
                     std::function<T(T, T)> op)
    {
        return mask_arith2(src, k, a, b, op, true);
    }

    template <typename T>
    ZMM arith2(ZMM a,
               ZMM b,
               std::function<T(T, T)> op)
    {
        return mask_arith2(ZMM(), Mask::full(), a, b, op);
    }

    //

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_add_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_add_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    //

    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_sub_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_sub_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    //

    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_mul_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_mul_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    //

    ZMM _mm512_div_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_div_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_div_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    //

    ZMM _mm512_min_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_min_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_min_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    //

    ZMM _mm512_max_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_max_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_max_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    //

    ZMM _mm512_pow_ps(ZMM a,
                      ZMM b);

    ZMM _mm512_mask_pow_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b);

    ZMM _mm512_maskz_pow_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b);

    // Arithmetic operations with 3 arguments.

    ZMM _mm512_fmadd_ps(ZMM a,
                        ZMM b,
                        ZMM c);

    // Compare operations.

    template <typename T>
    Mask compare(ZMM a,
                 ZMM b,
                 std::function<bool(T, T)> op)
    {
        Mask m;

        for (int i = 0; i < ZMM::count<T>(); i++)
        {
            m.set(i, op(a.get<T>(i), b.get<T>(i)));
        }

        return m;
    }

    //

    Mask _mm512_cmpeq_ps_mask(ZMM a,
                              ZMM b);

    //

    Mask _mm512_cmple_ps_mask(ZMM a,
                              ZMM b);

    //

    Mask _mm512_cmplt_ps_mask(ZMM a,
                              ZMM b);

    //

    Mask _mm512_cmpneq_ps_mask(ZMM a,
                               ZMM b);

    //

    Mask _mm512_cmpnle_ps_mask(ZMM a,
                               ZMM b);

    //

    Mask _mm512_cmpnlt_ps_mask(ZMM a,
                               ZMM b);

    Mask _mm512_mask_cmplt_ps_mask(Mask k,
                                   ZMM a,
                                   ZMM b);

    //

    Mask _mm512_cmpord_ps_mask(ZMM a,
                               ZMM b);

    //

    Mask _mm512_cmpunord_ps_mask(ZMM a,
                                 ZMM b);

    // Blend operations.

    ZMM _mm512_mask_blend_ps(Mask k,
                             ZMM a,
                             ZMM b);

    // Operations with masks.

    Mask koperation(int count,
                    Mask a,
                    Mask b,
                    std::function<bool(bool, bool)> op);

    Mask _mm512_kmov(Mask a);

    Mask _mm512_knot(Mask a);

    Mask _mm512_kand(Mask a,
                     Mask b);

    Mask _mm512_kandn(Mask a,
                      Mask b);

    Mask _mm512_kor(Mask a,
                    Mask b);

    Mask _mm512_kxor(Mask a,
                     Mask b);

    Mask _mm512_kxnor(Mask a,
                      Mask b);

    // Init operations.

    ZMM _mm512_set1_ps(float a);
}

#endif
