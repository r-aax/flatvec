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

    // Arithmetic operations with 2 arguments.

    ZMM _mm512_add_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, a.get<float>(i) + b.get<float>(i));
        }

        return r;
    }

    ZMM _mm512_sub_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, a.get<float>(i) - b.get<float>(i));
        }

        return r;
    }

    ZMM _mm512_mul_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, a.get<float>(i) * b.get<float>(i));
        }

        return r;
    }

    ZMM _mm512_div_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, a.get<float>(i) / b.get<float>(i));
        }

        return r;
    }

    ZMM _mm512_min_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, std::min(a.get<float>(i), b.get<float>(i)));
        }

        return r;
    }

    ZMM _mm512_max_ps(ZMM a,
                      ZMM b)
    {
        ZMM r;

        for (int i = 0; i < ZMM::count<float>(); i++)
        {
            r.set<float>(i, std::max(a.get<float>(i), b.get<float>(i)));
        }

        return r;
    }
}
