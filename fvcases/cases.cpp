#include "stdafx.h"

#include "cases.h"

#include "array_manager.h"

namespace fv
{
    // arith_f32

    /// <summary>
    /// Case arithmetic scalar.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="a">Input array.</param>
    /// <param name="b">Input array.</param>
    /// <param name="c">Output array.</param>
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

    /// <summary>
    /// Case arithmetic vector.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="a_p">Input array.</param>
    /// <param name="b_p">Input array.</param>
    /// <param name="c_p">Output array.</param>
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

    /// <summary>
    /// Case arithmetic.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
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

        bool res = (vc.max_diff(sc) == 0.0);

        if (!res)
        {
            std::cout << "max_diff : " << vc.max_diff(sc) << std::endl;
        }

        return res;
    }

    // blend_f32

    /// <summary>
    /// Case blend scalar.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="a">Input array.</param>
    /// <param name="b">Input array.</param>
    /// <param name="c">Output array.</param>
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

    /// <summary>
    /// Case blend vector.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="a_p">Input array.</param>
    /// <param name="b_p">Input array.</param>
    /// <param name="c_p">Output array.</param>
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

    /// <summary>
    /// Case blend.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
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

        bool res = (vc.max_diff(sc) == 0.0);

        if (!res)
        {
            std::cout << "max_diff : " << vc.max_diff(sc) << std::endl;
        }

        return res;
    }

    // Riemann solver
    // guessp function
    
    /// <summary>
    /// Case guessp scalar.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="dl">Input array.</param>
    /// <param name="ul">Input array.</param>
    /// <param name="pl">Input array.</param>
    /// <param name="cl">Input array.</param>
    /// <param name="dr">Input array.</param>
    /// <param name="ur">Input array.</param>
    /// <param name="pr">Input array.</param>
    /// <param name="cr">Input array.</param>
    /// <param name="pm">Output array.</param>
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

    /// <summary>
    /// Case guessp vector.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="dl_p">Input array.</param>
    /// <param name="ul_p">Input array.</param>
    /// <param name="pl_p">Input array.</param>
    /// <param name="cl_p">Input array.</param>
    /// <param name="dr_p">Input array.</param>
    /// <param name="ur_p">Input array.</param>
    /// <param name="pr_p">Input array.</param>
    /// <param name="cr_p">Input array.</param>
    /// <param name="pm_p">Output array.</param>
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

        ZMM z = _mm512_set1_ps(0.0f);
        ZMM one = _mm512_set1_ps(1.0f);
        ZMM two = _mm512_set1_ps(2.0f);
        ZMM half = _mm512_set1_ps(0.5f);
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

    /// <summary>
    /// Case guessp.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <param name="eps"></param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
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

        bool res = (vpm.max_diff(spm) < eps);

        if (!res)
        {
            std::cout << "max_diff : " << vpm.max_diff(spm) << std::endl;
        }

        return res;
    }

    // Riemann solver
    // prefun function.

    /// <summary>
    /// Case prefun scalar.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="f">Output array.</param>
    /// <param name="fd">Output array.</param>
    /// <param name="p">Input array.</param>
    /// <param name="dk">Input array.</param>
    /// <param name="pk">Input array.</param>
    /// <param name="ck">Input array.</param>
    void scase_prefun(int n,
                      float* f,
                      float* fd,
                      float* p,
                      float* dk,
                      float* pk,
                      float* ck)
    {
        float g {1.4f};
        float g1 = (g - 1.0f) / (2.0f * g);
        float g2 = (g + 1.0f) / (2.0f * g);
        float g4 = 2.0f / (g - 1.0f);
        float g5 = 2.0f / (g + 1.0f);
        float g6 = (g - 1.0f) / (g + 1.0f);

        for (int i = 0; i < n; i++)
        {
            float ak, bk, pratio, qrt;

            if (p[i] <= pk[i])
            {
                // Rarefaction wave.
                pratio = p[i] / pk[i];
                f[i] = g4 * ck[i] * (pow(pratio, g1) - 1.0f);
                fd[i] = (1.0f / (dk[i] * ck[i])) * pow(pratio, -g2);
            }
            else
            {
                // Shock wave.
                ak = g5 / dk[i];
                bk = g6 * pk[i];
                qrt = sqrt(ak / (bk + p[i]));
                f[i] = (p[i] - pk[i]) * qrt;
                fd[i] = (1.0f - 0.5f * (p[i] - pk[i]) / (bk + p[i])) * qrt;
            }
        }
    }

    /// <summary>
    /// Case prefun vector.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="f_p">Output array.</param>
    /// <param name="fd_p">Output array.</param>
    /// <param name="p_p">Input array.</param>
    /// <param name="dk_p">Input array.</param>
    /// <param name="pk_p">Input array.</param>
    /// <param name="ck_p">Input array.</param>
    void vcase_prefun(int n,
                      float* f_p,
                      float* fd_p,
                      float* p_p,
                      float* dk_p,
                      float* pk_p,
                      float* ck_p)
    {
        float sg {1.4f};
        float sg1 = (sg - 1.0f) / (2.0f * sg);
        float sg2 = (sg + 1.0f) / (2.0f * sg);
        float sg4 = 2.0f / (sg - 1.0f);
        float sg5 = 2.0f / (sg + 1.0f);
        float sg6 = (sg - 1.0f) / (sg + 1.0f);

        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        ZMM z = _mm512_set1_ps(0.0f);
        ZMM one = _mm512_set1_ps(1.0f);
        ZMM half = _mm512_set1_ps(0.5f);
        ZMM g1 = _mm512_set1_ps(sg1);
        ZMM g2 = _mm512_set1_ps(sg2);
        ZMM g4 = _mm512_set1_ps(sg4);
        ZMM g5 = _mm512_set1_ps(sg5);
        ZMM g6 = _mm512_set1_ps(sg6);

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();

            ZMM p = _mm512_load_ps(p_p + sh);
            ZMM dk = _mm512_load_ps(dk_p + sh);
            ZMM pk = _mm512_load_ps(pk_p + sh);
            ZMM ck = _mm512_load_ps(ck_p + sh);
            ZMM f, fd;

            ZMM pratio, ak, bkp, ppk, qrt;
            Mask cond, ncond;

            // Conditions.
            cond = _mm512_cmple_ps_mask(p, pk);
            ncond = _mm512_knot(cond);

            // The first branch.
            if (!cond.is_empty())
            {
                pratio = _mm512_mask_div_ps(z, cond, p, pk);
                f = _mm512_mask_mul_ps(f, cond,
                                       _mm512_mul_ps(g4, ck),
                                       _mm512_sub_ps(_mm512_mask_pow_ps(z, cond, pratio, g1), one));
                fd = _mm512_mask_div_ps(fd, cond,
                                        _mm512_mask_pow_ps(z, cond, pratio, _mm512_sub_ps(z, g2)),
                                        _mm512_mul_ps(dk, ck));
            }

            // The second branch.
            if (!ncond.is_empty())
            {
                ak = _mm512_mask_div_ps(z, ncond, g5, dk);
                bkp = _mm512_fmadd_ps(g6, pk, p);
                ppk = _mm512_sub_ps(p, pk);
                qrt = _mm512_mask_sqrt_ps(z, ncond,
                                          _mm512_mask_div_ps(z, ncond, ak, bkp));
                f = _mm512_mask_mul_ps(f, ncond, ppk, qrt);
                fd = _mm512_mask_mul_ps(fd, ncond, qrt,
                                        _mm512_fnmadd_ps(_mm512_mask_div_ps(z, ncond, ppk, bkp),
                                        _mm512_set1_ps(0.5), one));
            }

            _mm512_store_ps(f_p + sh, f);
            _mm512_store_ps(fd_p + sh, fd);
        }
    }

    /// <summary>
    /// Case prefun.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <param name="eps"></param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_prefun(int len,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * ZMM::count<float>();

        ArrayManager<float> sf(n);
        ArrayManager<float> vf(n);
        ArrayManager<float> sfd(n);
        ArrayManager<float> vfd(n);
        ArrayManager<float> p(n);
        ArrayManager<float> dk(n);
        ArrayManager<float> pk(n);
        ArrayManager<float> ck(n);

        p.generate_random(random_lo, random_hi);
        dk.generate_random(random_lo, random_hi);
        pk.generate_random(random_lo, random_hi);
        ck.generate_random(random_lo, random_hi);

        scase_prefun(n,
                     sf.get_data(), sfd.get_data(), p.get_data(),
                     dk.get_data(), pk.get_data(), ck.get_data());
        vcase_prefun(n,
                     vf.get_data(), vfd.get_data(), p.get_data(),
                     dk.get_data(), pk.get_data(), ck.get_data());

        bool res = (vf.max_diff(sf) + vfd.max_diff(sfd) < eps);

        if (!res)
        {
            std::cout << "max_diff : " << vf.max_diff(sf) << ", " << vfd.max_diff(sfd) << std::endl;
        }

        return res;
    }
}
