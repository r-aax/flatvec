#include "stdafx.h"

#include "cases.h"

#include "array_manager.h"

namespace fv
{
    // arith_f32

    void scase_arith_f32(int n,
                         float* a,
                         float* b,
                         float* c)
    {
        for (int i = 0; i < n; i++)
        {
            c[i] = (a[i] + b[i]) + (a[i] - b[i])
                + (a[i] * b[i]) + (a[i] / b[i])
                + std::min(a[i], b[i]) + std::max(a[i], b[i]);
        }
    }

    void vcase_arith_f32(int n,
                         float* a,
                         float* b,
                         float* c)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM va = _mm512_load_ps(a + sh);
            ZMM vb = _mm512_load_ps(b + sh);

            ZMM vadd = _mm512_add_ps(va, vb);
            ZMM vsub = _mm512_sub_ps(va, vb);
            ZMM vmul = _mm512_mul_ps(va, vb);
            ZMM vdiv = _mm512_div_ps(va, vb);
            ZMM vmin = _mm512_min_ps(va, vb);
            ZMM vmax = _mm512_max_ps(va, vb);

            ZMM vc = _mm512_add_ps(vadd, vsub);
            vc = _mm512_add_ps(vc, vmul);
            vc = _mm512_add_ps(vc, vdiv);
            vc = _mm512_add_ps(vc, vmin);
            vc = _mm512_add_ps(vc, vmax);

            _mm512_store_ps(c + sh, vc);
        }
    }

    bool case_arith_f32(int len,
                        float random_lo,
                        float random_hi)
    {
        int n = len * ZMM::count<float>();

        ArrayManager<float> a(n);
        ArrayManager<float> b(n);
        ArrayManager<float> sc(n);
        ArrayManager<float> vc(n);

        a.generate_random(random_lo, random_hi);
        b.generate_random(random_lo, random_hi);

        scase_arith_f32(n, a.get_data(), b.get_data(), sc.get_data());
        vcase_arith_f32(n, a.get_data(), b.get_data(), vc.get_data());

        return vc.maxDiff(sc) == 0.0;
    }

    // blend_f32

    void scase_blend_f32(int n,
                         float* a,
                         float* b,
                         float* c)
    {
        for (int i = 0; i < n; i++)
        {
            if (a[i] > b[i])
            {
                c[i] = a[i] + b[i];
            }
            else
            {
                c[i] = a[i] * b[i];
            }
        }
    }

    void vcase_blend_f32(int n,
                         float* a,
                         float* b,
                         float* c)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM va = _mm512_load_ps(a + sh);
            ZMM vb = _mm512_load_ps(b + sh);
            Mask k = _mm512_cmpgt_ps_mask(va, vb);
            ZMM vadd, vmul, vc;

            vadd = _mm512_maskz_add_ps(vadd, k, va, vb);
            vmul = _mm512_maskz_mul_ps(vmul, _mm512_knot(k), va, vb);
            vc = _mm512_mask_blend_ps(k, vmul, vadd);

            _mm512_store_ps(c + sh, vc);
        }
    }

    bool case_blend_f32(int len,
                        float random_lo,
                        float random_hi)
    {
        int n = len * ZMM::count<float>();

        ArrayManager<float> a(n);
        ArrayManager<float> b(n);
        ArrayManager<float> sc(n);
        ArrayManager<float> vc(n);

        a.generate_random(random_lo, random_hi);
        b.generate_random(random_lo, random_hi);

        scase_blend_f32(n, a.get_data(), b.get_data(), sc.get_data());
        vcase_blend_f32(n, a.get_data(), b.get_data(), vc.get_data());

        return vc.maxDiff(sc) == 0.0;
    }

    // Riemann solver
    // guessp function
    
    void scase_guessp(int n,
                      float* dl,
                      float* ul,
                      float* pl,
                      float* cl,
                      float* dr,
                      float* ur,
                      float* pr,
                      float* cr,
                      float* pm)
    {
        float g {1.4f};
        float g1 = (g - 1.0f) / (2.0f * g);
        float g3 = 2.0f * g / (g - 1.0f);
        float g4 = 2.0f / (g - 1.0f);
        float g5 = 2.0f / (g + 1.0f);
        float g6 = (g - 1.0f) / (g + 1.0f);
        float g7 = (g - 1.0f) / 2.0f;

        for (int i = 0; i < n; i++)
        {
            float cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr, qmax, quser, um;

            quser = 2.0f;

            // Compute guess pressure from PVRS Riemann solver.
            cup = 0.25f * (dl[i] + dr[i]) * (cl[i] + cr[i]);
            ppv = 0.5f * (pl[i] + pr[i]) + 0.5f * (ul[i] - ur[i]) * cup;
            ppv = (ppv > 0.0f) ? ppv : 0.0f;
            pmin = (pl[i] < pr[i]) ? pl[i] : pr[i];
            pmax = (pl[i] > pr[i]) ? pl[i] : pr[i];
            qmax = pmax / pmin;

            if ((qmax <= quser) && (pmin <= ppv) && (ppv <= pmax))
            {
                // Select PVRS Riemann solver.
                pm[i] = ppv;
            }
            else
            {
                if (ppv < pmin)
                {
                    // Select Two-Rarefaction Riemann solver.
                    pq = pow(pl[i] / pr[i], g1);
                    um = (pq * ul[i] / cl[i] + ur[i] / cr[i] + g4 * (pq - 1.0f)) / (pq / cl[i] + 1.0f / cr[i]);
                    ptl = 1.0f + g7 * (ul[i] - um) / cl[i];
                    ptr = 1.0f + g7 * (um - ur[i]) / cr[i];
                    pm[i] = 0.5f * (pow(pl[i] * ptl, g3) + pow(pr[i] * ptr, g3));
                }
                else
                {
                    // Select Two-Shock Riemann solver with PVRS as estimate.
                    gel = sqrt((g5 / dl[i]) / (g6 * pl[i] + ppv));
                    ger = sqrt((g5 / dr[i]) / (g6 * pr[i] + ppv));
                    pm[i] = (gel * pl[i] + ger * pr[i] - (ur[i] - ul[i])) / (gel + ger);
                }
            }
        }
    }

    void vcase_guessp(int n,
                      float* dl,
                      float* ul,
                      float* pl,
                      float* cl,
                      float* dr,
                      float* ur,
                      float* pr,
                      float* cr,
                      float* pm)
    {
        ;
    }

    bool case_guessp(int len)
    {
        return false;
    }
}
