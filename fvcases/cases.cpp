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
                         float* a_p,
                         float* b_p,
                         float* c_p)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM a = _mm512_load_ps(a_p + sh);
            ZMM b = _mm512_load_ps(b_p + sh);

            ZMM add = _mm512_add_ps(a, b);
            ZMM sub = _mm512_sub_ps(a, b);
            ZMM mul = _mm512_mul_ps(a, b);
            ZMM div = _mm512_div_ps(a, b);
            ZMM min = _mm512_min_ps(a, b);
            ZMM max = _mm512_max_ps(a, b);

            ZMM c = _mm512_add_ps(add, sub);
            c = _mm512_add_ps(c, mul);
            c = _mm512_add_ps(c, div);
            c = _mm512_add_ps(c, min);
            c = _mm512_add_ps(c, max);

            _mm512_store_ps(c_p + sh, c);
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

        return vc.max_diff(sc) == 0.0;
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
                         float* a_p,
                         float* b_p,
                         float* c_p)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM a = _mm512_load_ps(a_p + sh);
            ZMM b = _mm512_load_ps(b_p + sh);
            Mask k = _mm512_cmpgt_ps_mask(a, b);
            ZMM add, mul, c;

            add = _mm512_maskz_add_ps(add, k, a, b);
            mul = _mm512_maskz_mul_ps(mul, _mm512_knot(k), a, b);
            c = _mm512_mask_blend_ps(k, mul, add);

            _mm512_store_ps(c_p + sh, c);
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

        return vc.max_diff(sc) == 0.0;
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
                      float* dl_p,
                      float* ul_p,
                      float* pl_p,
                      float* cl_p,
                      float* dr_p,
                      float* ur_p,
                      float* pr_p,
                      float* cr_p,
                      float* pm_p)
    {
        float sg {1.4f};
        float sg1 = (sg - 1.0f) / (2.0f * sg);
        float sg3 = 2.0f * sg / (sg - 1.0f);
        float sg4 = 2.0f / (sg - 1.0f);
        float sg5 = 2.0f / (sg + 1.0f);
        float sg6 = (sg - 1.0f) / (sg + 1.0f);
        float sg7 = (sg - 1.0f) / 2.0f;

        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        ZMM z = _mm512_set1_ps(0.0);
        ZMM one = _mm512_set1_ps(1.0);
        ZMM two = _mm512_set1_ps(2.0);
        ZMM half = _mm512_set1_ps(0.5);
        ZMM g1 = _mm512_set1_ps(sg1);
        ZMM g3 = _mm512_set1_ps(sg3);
        ZMM g4 = _mm512_set1_ps(sg4);
        ZMM g5 = _mm512_set1_ps(sg5);
        ZMM g6 = _mm512_set1_ps(sg6);
        ZMM g7 = _mm512_set1_ps(sg7);

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM dl = _mm512_load_ps(dl_p + sh);
            ZMM ul = _mm512_load_ps(ul_p + sh);
            ZMM pl = _mm512_load_ps(pl_p + sh);
            ZMM cl = _mm512_load_ps(cl_p + sh);
            ZMM dr = _mm512_load_ps(dr_p + sh);
            ZMM ur = _mm512_load_ps(ur_p + sh);
            ZMM pr = _mm512_load_ps(pr_p + sh);
            ZMM cr = _mm512_load_ps(cr_p + sh);
            ZMM pm;

            // Begin of calculation part.

            ZMM cup, ppv, pmin, pmax, qmax, pq, um, ptl, ptr, gel, ger, pqcr;
            Mask cond_pvrs, cond_ppv, ncond_ppv;

            cup = _mm512_mul_ps(_mm512_set1_ps(0.25), _mm512_mul_ps(_mm512_add_ps(dl, dr), _mm512_add_ps(cl, cr)));
            ppv = _mm512_mul_ps(half, _mm512_fmadd_ps(_mm512_sub_ps(ul, ur), cup, _mm512_add_ps(pl, pr)));
            ppv = _mm512_max_ps(ppv, z);
            pmin = _mm512_min_ps(pl, pr);
            pmax = _mm512_max_ps(pl, pr);
            qmax = _mm512_div_ps(pmax, pmin);

            // Conditions.
            cond_pvrs = _mm512_kand(_mm512_kand(_mm512_cmple_ps_mask(qmax, two),
                                                _mm512_cmple_ps_mask(pmin, ppv)),
                                    _mm512_cmple_ps_mask(ppv, pmax));
            cond_ppv = _mm512_mask_cmplt_ps_mask(_mm512_knot(cond_pvrs), ppv, pmin);
            ncond_ppv = _mm512_kand(_mm512_knot(cond_pvrs), _mm512_knot(cond_ppv));

            // The first branch.
            pm = _mm512_mask_mov_ps(pm, cond_pvrs, ppv);

            // The second branch.
            if (!cond_ppv.is_empty())
            {
                pq = _mm512_mask_pow_ps(z, cond_ppv, _mm512_mask_div_ps(z, cond_ppv, pl, pr), g1);
                pqcr = _mm512_mul_ps(pq, cr);
                um = _mm512_mask_div_ps(z, cond_ppv,
                                        _mm512_fmadd_ps(_mm512_fmadd_ps(_mm512_sub_ps(pqcr, cr), g4, ur), cl, _mm512_mul_ps(pqcr, ul)),
                                        _mm512_add_ps(pqcr, cl));
                ptl = _mm512_fmadd_ps(_mm512_mask_div_ps(z, cond_ppv, _mm512_sub_ps(ul, um), cl), g7, one);
                ptr = _mm512_fmadd_ps(_mm512_mask_div_ps(z, cond_ppv, _mm512_sub_ps(um, ur), cr), g7, one);
                pm = _mm512_mask_mul_ps(pm, cond_ppv, half,
                                        _mm512_add_ps(_mm512_mask_pow_ps(z, cond_ppv, _mm512_mul_ps(pl, ptl), g3),
                                        _mm512_mask_pow_ps(z, cond_ppv, _mm512_mul_ps(pr, ptr), g3)));
            }

            // The third branch.
            if (!ncond_ppv.is_empty())
            {
                gel = _mm512_sqrt_ps(_mm512_mask_div_ps(z, ncond_ppv, g5, _mm512_mul_ps(_mm512_fmadd_ps(g6, pl, ppv), dl)));
                ger = _mm512_sqrt_ps(_mm512_mask_div_ps(z, ncond_ppv, g5, _mm512_mul_ps(_mm512_fmadd_ps(g6, pr, ppv), dr)));
                pm = _mm512_mask_div_ps(pm, ncond_ppv,
                                        _mm512_fmadd_ps(gel, pl, _mm512_fmadd_ps(ger, pr, _mm512_sub_ps(ul, ur))),
                                        _mm512_add_ps(gel, ger));
            }

            // End of calculation part.

            _mm512_store_ps(pm_p + sh, pm);
        }
    }

    bool case_guessp(int len,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * ZMM::count<float>();

        ArrayManager<float> dl(n);
        ArrayManager<float> ul(n);
        ArrayManager<float> pl(n);
        ArrayManager<float> cl(n);
        ArrayManager<float> dr(n);
        ArrayManager<float> ur(n);
        ArrayManager<float> pr(n);
        ArrayManager<float> cr(n);
        ArrayManager<float> spm(n);
        ArrayManager<float> vpm(n);

        dl.generate_random(random_lo, random_hi);
        ul.generate_random(random_lo, random_hi);
        pl.generate_random(random_lo, random_hi);
        cl.generate_random(random_lo, random_hi);
        dr.generate_random(random_lo, random_hi);
        ur.generate_random(random_lo, random_hi);
        pr.generate_random(random_lo, random_hi);
        cr.generate_random(random_lo, random_hi);

        scase_guessp(n,
                     dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                     dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                     spm.get_data());
        vcase_guessp(n,
                     dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                     dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                     vpm.get_data());

        return vpm.max_diff(spm) < eps;
    }
}
