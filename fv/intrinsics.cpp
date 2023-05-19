#include "stdafx.h"

#include "intrinsics.h"

namespace fv
{
    // Memory access.

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

    ZMM _mm512_mask_mov_ps(ZMM src,
                           Mask k,
                           ZMM a)
    {
        std::cout << "_mm512_mask_mov_ps is not implemented" << std::endl;

        return ZMM();
    }

    ZMM _mm512_sqrt_ps(ZMM a)
    {
        std::cout << "_mm512_sqrt_ps is not implemented" << std::endl;

        return ZMM();
    }

    // Arithmetic operations with 2 arguments.

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::plus<float>());
    }

    ZMM _mm512_mask_add_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::plus<float>());
    }

    ZMM _mm512_maskz_add_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::plus<float>());
    }

    //

    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::minus<float>());
    }

    ZMM _mm512_mask_sub_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::minus<float>());
    }

    ZMM _mm512_maskz_sub_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::minus<float>());
    }

    //

    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::multiplies<float>());
    }

    ZMM _mm512_mask_mul_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::multiplies<float>());
    }

    ZMM _mm512_maskz_mul_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::multiplies<float>());
    }

    //

    ZMM _mm512_div_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, std::divides<float>());
    }

    ZMM _mm512_mask_div_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, std::divides<float>());
    }

    ZMM _mm512_maskz_div_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, std::divides<float>());
    }

    //

    ZMM _mm512_min_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return std::min(x, y); });
    }

    ZMM _mm512_mask_min_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, [] (float x, float y) { return std::min(x, y); });
    }

    ZMM _mm512_maskz_min_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, [] (float x, float y) { return std::min(x, y); });
    }

    //

    ZMM _mm512_max_ps(ZMM a,
                      ZMM b)
    {
        return arith2<float>(a, b, [] (float x, float y) { return std::max(x, y); });
    }

    ZMM _mm512_mask_max_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        return mask_arith2<float>(src, k, a, b, [] (float x, float y) { return std::max(x, y); });
    }

    ZMM _mm512_maskz_max_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        return maskz_arith2<float>(src, k, a, b, [] (float x, float y) { return std::max(x, y); });
    }

    //

    ZMM _mm512_pow_ps(ZMM a,
                      ZMM b)
    {
        std::cout << "_mm512_pow_ps is not implemented" << std::endl;

        return ZMM();
    }

    ZMM _mm512_mask_pow_ps(ZMM src,
                           Mask k,
                           ZMM a,
                           ZMM b)
    {
        std::cout << "_mm512_mask_pow_ps is not implemented" << std::endl;

        return ZMM();
    }

    ZMM _mm512_maskz_pow_ps(ZMM src,
                            Mask k,
                            ZMM a,
                            ZMM b)
    {
        std::cout << "_mm512_maskz_pow_ps is not implemented" << std::endl;

        return ZMM();
    }

    // Arithmetic operations with 3 arguments.

    ZMM _mm512_fmadd_ps(ZMM a,
                        ZMM b,
                        ZMM c)
    {
        std::cout << "_mm512_fmadd_ps is not implemented" << std::endl;

        return ZMM();
    }

    // Compare operations.

    Mask _mm512_cmpeq_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x == y; });
    }

    //

    Mask _mm512_cmple_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x <= y; });
    }

    //

    Mask _mm512_cmplt_ps_mask(ZMM a,
                              ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x < y; });
    }

    Mask _mm512_mask_cmplt_ps_mask(Mask k,
                                   ZMM a,
                                   ZMM b)
    {
        std::cout << "_mm512_mask_cmplt_ps_mask is not implemented" << std::endl;
    
        return Mask();
    }

    //

    Mask _mm512_cmpneq_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return x != y; });
    }

    //

    Mask _mm512_cmpnle_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !(x <= y); });
    }

    //

    Mask _mm512_cmpnlt_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !(x < y); });
    }

    //

    Mask _mm512_cmpord_ps_mask(ZMM a,
                               ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return !std::isnan(x) && !std::isnan(y); });
    }

    //

    Mask _mm512_cmpunord_ps_mask(ZMM a,
                                 ZMM b)
    {
        return compare<float>(a, b, [] (float x, float y) { return std::isnan(x) || std::isnan(y); });
    }

    // Blend operations.

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

    Mask _mm512_kmov(Mask a)
    {
        return koperation(ZMM::count<float>(), a, a, [] (bool x, bool y) { return x; });
    }

    Mask _mm512_knot(Mask a)
    {
        return koperation(ZMM::count<float>(), a, a, [] (bool x, bool y) { return !x; });
    }

    Mask _mm512_kand(Mask a,
                     Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, std::bit_and<bool>());
    }

    Mask _mm512_kandn(Mask a,
                      Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, [] (bool x, bool y) { return !x & y; });
    }

    Mask _mm512_kor(Mask a,
                    Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, std::bit_or<bool>());
    }

    Mask _mm512_kxor(Mask a,
                     Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, std::bit_xor<bool>());
    }

    Mask _mm512_kxnor(Mask a,
                      Mask b)
    {
        return koperation(ZMM::count<float>(), a, b, [] (bool x, bool y) { return !(x ^ y); });
    }

    // Init operations.

    ZMM _mm512_set1_ps(float a)
    {
        std::cout << "_mm512_set1_ps is not implemented" << std::endl;

        return ZMM();
    }
}
