#include "stdafx.h"

#include "cases.h"

#include "array_manager.h"

namespace fv
{
    // Some constants

    ZMM zero = _mm512_set1_ps(0.0f);
    ZMM one = _mm512_set1_ps(1.0f);
    ZMM two = _mm512_set1_ps(2.0f);
    ZMM half = _mm512_set1_ps(0.5f);

    // arith_f32

    /// <summary>
    /// Case arithmetic scalar.
    /// </summary>
    /// <param name="a">Input.</param>
    /// <param name="b">Input.</param>
    /// <param name="c">Output.</param>
    void scase_arith_f32_1(float a,
                           float b,
                           float& c)
    {
        c = (a + b) + (a - b)
            + (a * b) + (a / b)
            + std::min(a, b) + std::max(a, b);
    }

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
            scase_arith_f32_1(a[i], b[i], c[i]);
        }
    }

    /// <summary>
    /// Case arithmetic vector.
    /// </summary>
    /// <param name="a">Input.</param>
    /// <param name="b">Input.</param>
    /// <param name="c">Output.</param>
    void vcase_arith_f32_1(ZMM a,
                           ZMM b,
                           ZMM& c)
    {
        ZMM add = _mm512_add_ps(a, b);
        ZMM sub = _mm512_sub_ps(a, b);
        ZMM mul = _mm512_mul_ps(a, b);
        ZMM div = _mm512_div_ps(a, b);
        ZMM min = _mm512_min_ps(a, b);
        ZMM max = _mm512_max_ps(a, b);

        c = _mm512_add_ps(add, sub);
        c = _mm512_add_ps(c, mul);
        c = _mm512_add_ps(c, div);
        c = _mm512_add_ps(c, min);
        c = _mm512_add_ps(c, max);
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
            ZMM c;

            vcase_arith_f32_1(_mm512_load_ps(a_p + sh),
                              _mm512_load_ps(b_p + sh),
                              c);

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
    /// <param name="a">Input.</param>
    /// <param name="b">Input.</param>
    /// <param name="c">Output.</param>
    void scase_blend_f32_1(float a,
                           float b,
                           float& c)
    {
        if (a > b)
        {
            c = a + b;
        }
        else
        {
            c = a * b;
        }
    }

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
            scase_blend_f32_1(a[i], b[i], c[i]);
        }
    }

    /// <summary>
    /// Case blend vector.
    /// </summary>
    /// <param name="a">Input.</param>
    /// <param name="b">Input.</param>
    /// <param name="c">Output.</param>
    void vcase_blend_f32_1(ZMM a,
                           ZMM b,
                           ZMM& c)
    {
        Mask k = _mm512_cmpgt_ps_mask(a, b);
        ZMM add, mul;

        add = _mm512_maskz_add_ps(add, k, a, b);
        mul = _mm512_maskz_mul_ps(mul, _mm512_knot(k), a, b);
        c = _mm512_mask_blend_ps(k, mul, add);
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
            ZMM c;

            vcase_blend_f32_1(_mm512_load_ps(a_p + sh),
                              _mm512_load_ps(b_p + sh),
                              c);

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

    // Constants for Riemann solver

    namespace riemann
    {
        float sg {1.4f};
        float sg1 = (sg - 1.0f) / (2.0f * sg);
        float sg2 = (sg + 1.0f) / (2.0f * sg);
        float sg3 = 2.0f * sg / (sg - 1.0f);
        float sg4 = 2.0f / (sg - 1.0f);
        float sg5 = 2.0f / (sg + 1.0f);
        float sg6 = (sg - 1.0f) / (sg + 1.0f);
        float sg7 = (sg - 1.0f) / 2.0f;
        float sg8 = (sg - 1.0f);

        ZMM g1 = _mm512_set1_ps(sg1);
        ZMM g2 = _mm512_set1_ps(sg2);
        ZMM g3 = _mm512_set1_ps(sg3);
        ZMM g4 = _mm512_set1_ps(sg4);
        ZMM g5 = _mm512_set1_ps(sg5);
        ZMM g6 = _mm512_set1_ps(sg6);
        ZMM g7 = _mm512_set1_ps(sg7);
        ZMM g8 = _mm512_set1_ps(sg8);
    }

    // Riemann solver
    // guessp function
    
    /// <summary>
    /// Case guessp scalar.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="cl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="cr">Input.</param>
    /// <param name="pm">Output.</param>
    void scase_guessp_1(float dl,
                        float ul,
                        float pl,
                        float cl,
                        float dr,
                        float ur,
                        float pr,
                        float cr,
                        float& pm)
    {
        float cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr, qmax, quser, um;

        quser = 2.0f;

        // Compute guess pressure from PVRS Riemann solver.
        cup = 0.25f * (dl + dr) * (cl + cr);
        ppv = 0.5f * (pl + pr) + 0.5f * (ul - ur) * cup;
        ppv = (ppv > 0.0f) ? ppv : 0.0f;
        pmin = (pl < pr) ? pl : pr;
        pmax = (pl > pr) ? pl : pr;
        qmax = pmax / pmin;

        if ((qmax <= quser) && (pmin <= ppv) && (ppv <= pmax))
        {
            // Select PVRS Riemann solver.
            pm = ppv;
        }
        else
        {
            if (ppv < pmin)
            {
                // Select Two-Rarefaction Riemann solver.
                pq = pow(pl / pr, riemann::sg1);
                um = (pq * ul / cl + ur / cr + riemann::sg4 * (pq - 1.0f)) / (pq / cl + 1.0f / cr);
                ptl = 1.0f + riemann::sg7 * (ul - um) / cl;
                ptr = 1.0f + riemann::sg7 * (um - ur) / cr;
                pm = 0.5f * (pow(pl * ptl, riemann::sg3) + pow(pr * ptr, riemann::sg3));
            }
            else
            {
                // Select Two-Shock Riemann solver with PVRS as estimate.
                gel = sqrt((riemann::sg5 / dl) / (riemann::sg6 * pl + ppv));
                ger = sqrt((riemann::sg5 / dr) / (riemann::sg6 * pr + ppv));
                pm = (gel * pl + ger * pr - (ur - ul)) / (gel + ger);
            }
        }
    }

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
        for (int i = 0; i < n; i++)
        {
            scase_guessp_1(dl[i], ul[i], pl[i], cl[i], dr[i], ur[i], pr[i], cr[i], pm[i]);
        }
    }

    /// <summary>
    /// Case guessp vector.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="cl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="cr">Input.</param>
    /// <param name="pm">Output.</param>
    void vcase_guessp_1(ZMM dl,
                        ZMM ul,
                        ZMM pl,
                        ZMM cl,
                        ZMM dr,
                        ZMM ur,
                        ZMM pr,
                        ZMM cr,
                        ZMM& pm)
    {
        // Begin of calculation part.

        ZMM cup, ppv, pmin, pmax, qmax, pq, um, ptl, ptr, gel, ger, pqcr;
        Mask cond_pvrs, cond_ppv, ncond_ppv;

        cup = _mm512_mul_ps(_mm512_set1_ps(0.25), _mm512_mul_ps(_mm512_add_ps(dl, dr), _mm512_add_ps(cl, cr)));
        ppv = _mm512_mul_ps(half, _mm512_fmadd_ps(_mm512_sub_ps(ul, ur), cup, _mm512_add_ps(pl, pr)));
        ppv = _mm512_max_ps(ppv, zero);
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
            pq = _mm512_mask_pow_ps(zero, cond_ppv, _mm512_mask_div_ps(zero, cond_ppv, pl, pr), riemann::g1);
            pqcr = _mm512_mul_ps(pq, cr);
            um = _mm512_mask_div_ps(zero, cond_ppv,
                                    _mm512_fmadd_ps(_mm512_fmadd_ps(_mm512_sub_ps(pqcr, cr), riemann::g4, ur), cl, _mm512_mul_ps(pqcr, ul)),
                                    _mm512_add_ps(pqcr, cl));
            ptl = _mm512_fmadd_ps(_mm512_mask_div_ps(zero, cond_ppv, _mm512_sub_ps(ul, um), cl), riemann::g7, one);
            ptr = _mm512_fmadd_ps(_mm512_mask_div_ps(zero, cond_ppv, _mm512_sub_ps(um, ur), cr), riemann::g7, one);
            pm = _mm512_mask_mul_ps(pm, cond_ppv, half,
                                    _mm512_add_ps(_mm512_mask_pow_ps(zero, cond_ppv, _mm512_mul_ps(pl, ptl), riemann::g3),
                                    _mm512_mask_pow_ps(zero, cond_ppv, _mm512_mul_ps(pr, ptr), riemann::g3)));
        }

        // The third branch.
        if (!ncond_ppv.is_empty())
        {
            gel = _mm512_sqrt_ps(_mm512_mask_div_ps(zero, ncond_ppv, riemann::g5,
                                                    _mm512_mul_ps(_mm512_fmadd_ps(riemann::g6, pl, ppv), dl)));
            ger = _mm512_sqrt_ps(_mm512_mask_div_ps(zero, ncond_ppv, riemann::g5,
                                                    _mm512_mul_ps(_mm512_fmadd_ps(riemann::g6, pr, ppv), dr)));
            pm = _mm512_mask_div_ps(pm, ncond_ppv,
                                    _mm512_fmadd_ps(gel, pl, _mm512_fmadd_ps(ger, pr, _mm512_sub_ps(ul, ur))),
                                    _mm512_add_ps(gel, ger));
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
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();
            ZMM pm;

            vcase_guessp_1(_mm512_load_ps(dl_p + sh),
                           _mm512_load_ps(ul_p + sh),
                           _mm512_load_ps(pl_p + sh),
                           _mm512_load_ps(cl_p + sh),
                           _mm512_load_ps(dr_p + sh),
                           _mm512_load_ps(ur_p + sh),
                           _mm512_load_ps(pr_p + sh),
                           _mm512_load_ps(cr_p + sh),
                           pm);

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
    /// <param name="f">Output.</param>
    /// <param name="fd">Output.</param>
    /// <param name="p">Input.</param>
    /// <param name="dk">Input.</param>
    /// <param name="pk">Input.</param>
    /// <param name="ck">Input.</param>
    void scase_prefun_1(float& f,
                        float& fd,
                        float p,
                        float dk,
                        float pk,
                        float ck)
    {
        float ak, bk, pratio, qrt;

        if (p <= pk)
        {
            // Rarefaction wave.
            pratio = p / pk;
            f = riemann::sg4 * ck * (pow(pratio, riemann::sg1) - 1.0f);
            fd = (1.0f / (dk * ck)) * pow(pratio, -riemann::sg2);
        }
        else
        {
            // Shock wave.
            ak = riemann::sg5 / dk;
            bk = riemann::sg6 * pk;
            qrt = sqrt(ak / (bk + p));
            f = (p - pk) * qrt;
            fd = (1.0f - 0.5f * (p - pk) / (bk + p)) * qrt;
        }
    }

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
        for (int i = 0; i < n; i++)
        {
            scase_prefun_1(f[i], fd[i], p[i], dk[i], pk[i], ck[i]);
        }
    }

    /// <summary>
    /// Case prefun vector.
    /// </summary>
    /// <param name="f">Output.</param>
    /// <param name="fd">Output.</param>
    /// <param name="p">Input.</param>
    /// <param name="dk">Input.</param>
    /// <param name="pk">Input.</param>
    /// <param name="ck">Input.</param>
    /// <param name="m">Mask.</param>
    void vcase_prefun_1(ZMM& f,
                        ZMM& fd,
                        ZMM p,
                        ZMM dk,
                        ZMM pk,
                        ZMM ck,
                        Mask m)
    {
        ZMM pratio, ak, bkp, ppk, qrt;
        Mask cond, ncond;

        // Conditions.
        cond = _mm512_kand(_mm512_cmple_ps_mask(p, pk), m);
        ncond = _mm512_kand(_mm512_knot(cond), m);

        // The first branch.
        if (!cond.is_empty())
        {
            pratio = _mm512_mask_div_ps(zero, cond, p, pk);
            f = _mm512_mask_mul_ps(f, cond,
                                   _mm512_mul_ps(riemann::g4, ck),
                                   _mm512_sub_ps(_mm512_mask_pow_ps(zero, cond, pratio, riemann::g1), one));
            fd = _mm512_mask_div_ps(fd, cond,
                                    _mm512_mask_pow_ps(zero, cond, pratio, _mm512_sub_ps(zero, riemann::g2)),
                                    _mm512_mul_ps(dk, ck));
        }

        // The second branch.
        if (!ncond.is_empty())
        {
            ak = _mm512_mask_div_ps(zero, ncond, riemann::g5, dk);
            bkp = _mm512_fmadd_ps(riemann::g6, pk, p);
            ppk = _mm512_sub_ps(p, pk);
            qrt = _mm512_mask_sqrt_ps(zero, ncond,
                                      _mm512_mask_div_ps(zero, ncond, ak, bkp));
            f = _mm512_mask_mul_ps(f, ncond, ppk, qrt);
            fd = _mm512_mask_mul_ps(fd, ncond, qrt,
                                    _mm512_fnmadd_ps(_mm512_mask_div_ps(zero, ncond, ppk, bkp),
                                    _mm512_set1_ps(0.5), one));
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
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();
            ZMM f, fd;

            vcase_prefun_1(f, fd,
                           _mm512_load_ps(p_p + sh),
                           _mm512_load_ps(dk_p + sh),
                           _mm512_load_ps(pk_p + sh),
                           _mm512_load_ps(ck_p + sh),
                           Mask::full());

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

    // Riemann solver
    // sample function

    /// <summary>
    /// Case sample scalar.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="vl">Input.</param>
    /// <param name="wl">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="cl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="vr">Input.</param>
    /// <param name="wr">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="cr">Input.</param>
    /// <param name="pm">Input.</param>
    /// <param name="um">Input.</param>
    /// <param name="d">Output.</param>
    /// <param name="u">Output.</param>
    /// <param name="v">Output.</param>
    /// <param name="w">Output.</param>
    /// <param name="p">Output.</param>
    void scase_sample_1(float dl,
                        float ul,
                        float vl,
                        float wl,
                        float pl,
                        float cl,
                        float dr,
                        float ur,
                        float vr,
                        float wr,
                        float pr,
                        float cr,
                        float pm,
                        float um,
                        float& d,
                        float& u,
                        float& v,
                        float& w,
                        float& p)
    {
        float c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

        if (0.0f <= um)
        {
            // Sampling point lies to the left of the contact discontinuity.
            v = vl;
            w = wl;

            if (pm <= pl)
            {
                // Left rarefaction.
                shl = ul - cl;

                if (0.0f <= shl)
                {
                    // Sampled point is left data state.
                    d = dl;
                    u = ul;
                    p = pl;
                }
                else
                {
                    cml = cl * pow(pm / pl, riemann::sg1);
                    stl = um - cml;

                    if (0.0f > stl)
                    {
                        // Sampled point is star left state.
                        d = dl * pow(pm / pl, 1.0f / riemann::sg);
                        u = um;
                        p = pm;
                    }
                    else
                    {
                        // Sampled point is inside left fan.
                        u = riemann::sg5 * (cl + riemann::sg7 * ul);
                        c = riemann::sg5 * (cl + riemann::sg7 * ul);
                        d = dl * pow(c / cl, riemann::sg4);
                        p = pl * pow(c / cl, riemann::sg3);
                    }
                }
            }
            else
            {
                // Left shock.
                pml = pm / pl;
                sl = ul - cl * sqrt(riemann::sg2 * pml + riemann::sg1);

                if (0.0 <= sl)
                {
                    // Sampled point is left data state.
                    d = dl;
                    u = ul;
                    p = pl;
                }
                else
                {
                    // Sampled point is star left state.
                    d = dl * (pml + riemann::sg6) / (pml * riemann::sg6 + 1.0f);
                    u = um;
                    p = pm;
                }
            }
        }
        else
        {
            // Sampling point lies to the right of the contact discontinuity.
            v = vr;
            w = wr;

            if (pm > pr)
            {
                // Right shock.
                pmr = pm / pr;
                sr = ur + cr * sqrt(riemann::sg2 * pmr + riemann::sg1);

                if (0.0f >= sr)
                {
                    // Sampled point is right data state.
                    d = dr;
                    u = ur;
                    p = pr;
                }
                else
                {
                    // Sampled point is star right state.
                    d = dr * (pmr + riemann::sg6) / (pmr * riemann::sg6 + 1.0f);
                    u = um;
                    p = pm;
                }
            }
            else
            {
                // Right rarefaction.
                shr = ur + cr;
                if (0.0f >= shr)
                {
                    // Sampled point is right data state.
                    d = dr;
                    u = ur;
                    p = pr;
                }
                else
                {
                    cmr = cr * pow(pm / pr, riemann::sg1);
                    str = um + cmr;

                    if (0.0f <= str)
                    {
                        // Sampled point is star right state.
                        d = dr * pow(pm / pr, 1.0f / riemann::sg);
                        u = um;
                        p = pm;
                    }
                    else
                    {
                        // Sampled point is inside left fan.
                        u = riemann::sg5 * (-cr + riemann::sg7 * ur);
                        c = riemann::sg5 * (cr - riemann::sg7 * ur);
                        d = dr * pow(c / cr, riemann::sg4);
                        p = pr * pow(c / cr, riemann::sg3);
                    }
                }
            }
        }
    }

    /// <summary>
    /// Case sample scalar.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="dl">Input array.</param>
    /// <param name="ul">Input array.</param>
    /// <param name="vl">Input array.</param>
    /// <param name="wl">Input array.</param>
    /// <param name="pl">Input array.</param>
    /// <param name="cl">Input array.</param>
    /// <param name="dr">Input array.</param>
    /// <param name="ur">Input array.</param>
    /// <param name="vr">Input array.</param>
    /// <param name="wr">Input array.</param>
    /// <param name="pr">Input array.</param>
    /// <param name="cr">Input array.</param>
    /// <param name="pm">Input array.</param>
    /// <param name="um">Input array.</param>
    /// <param name="d">Output array.</param>
    /// <param name="u">Output array.</param>
    /// <param name="v">Output array.</param>
    /// <param name="w">Output array.</param>
    /// <param name="p">Output array.</param>
    void scase_sample(int n,
                      float* dl,
                      float* ul,
                      float* vl,
                      float* wl,
                      float* pl,
                      float* cl,
                      float* dr,
                      float* ur,
                      float* vr,
                      float* wr,
                      float* pr,
                      float* cr,
                      float* pm,
                      float* um,
                      float* d,
                      float* u,
                      float* v,
                      float* w,
                      float* p)
    {
        for (int i = 0; i < n; i++)
        {
            scase_sample_1(dl[i], ul[i], vl[i], wl[i], pl[i], cl[i],
                           dr[i], ur[i], vr[i], wr[i], pr[i], cr[i],
                           pm[i], um[i], d[i], u[i], v[i], w[i], p[i]);
        }
    }

    /// <summary>
    /// Case sample vector.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="vl">Input.</param>
    /// <param name="wl">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="cl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="vr">Input.</param>
    /// <param name="wr">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="cr">Input.</param>
    /// <param name="pm">Input.</param>
    /// <param name="um">Input.</param>
    /// <param name="d">Output.</param>
    /// <param name="u">Output.</param>
    /// <param name="v">Output.</param>
    /// <param name="w">Output.</param>
    /// <param name="p">Output.</param>
    void vcase_sample_1(ZMM dl,
                        ZMM ul,
                        ZMM vl,
                        ZMM wl,
                        ZMM pl,
                        ZMM cl,
                        ZMM dr,
                        ZMM ur,
                        ZMM vr,
                        ZMM wr,
                        ZMM pr,
                        ZMM cr,
                        ZMM pm,
                        ZMM um,
                        ZMM& d,
                        ZMM& u,
                        ZMM& v,
                        ZMM& w,
                        ZMM& p)
    {
        ZMM c, ums, pms, sh, st, s, uc;
        Mask cond_um, cond_pm, cond_sh, cond_st, cond_s, cond_sh_st;

        // d/u/p/c/ums
        cond_um = _mm512_cmplt_ps_mask(um, zero);
        d = _mm512_mask_blend_ps(cond_um, dl, dr);
        u = _mm512_mask_blend_ps(cond_um, ul, ur);
        v = _mm512_mask_blend_ps(cond_um, vl, vr);
        w = _mm512_mask_blend_ps(cond_um, wl, wr);
        p = _mm512_mask_blend_ps(cond_um, pl, pr);
        c = _mm512_mask_blend_ps(cond_um, cl, cr);
        ums = um;
        u = _mm512_mask_sub_ps(u, cond_um, zero, u);
        ums = _mm512_mask_sub_ps(ums, cond_um, zero, ums);

        // Calculate main values.
        pms = _mm512_div_ps(pm, p);
        sh = _mm512_sub_ps(u, c);
        st = _mm512_fnmadd_ps(_mm512_pow_ps(pms, riemann::g1), c, ums);
        s = _mm512_fnmadd_ps(c, _mm512_sqrt_ps(_mm512_fmadd_ps(riemann::g2, pms, riemann::g1)), u);

        // Conditions.
        cond_pm = _mm512_cmple_ps_mask(pm, p);
        cond_sh = _mm512_mask_cmplt_ps_mask(cond_pm, sh, zero);
        cond_st = _mm512_mask_cmplt_ps_mask(cond_sh, st, zero);
        cond_s = _mm512_mask_cmplt_ps_mask(_mm512_knot(cond_pm), s, zero);

        // Store.
        d = _mm512_mask_mov_ps(d, cond_st, _mm512_mul_ps(d, _mm512_pow_ps(pms, _mm512_set1_ps(1.0f / riemann::sg))));
        d = _mm512_mask_mov_ps(d, cond_s, _mm512_mul_ps(d, _mm512_div_ps(_mm512_add_ps(pms, riemann::g6),
                                                                         _mm512_fmadd_ps(pms, riemann::g6, one))));
        u = _mm512_mask_mov_ps(u, _mm512_kor(cond_st, cond_s), ums);
        p = _mm512_mask_mov_ps(p, _mm512_kor(cond_st, cond_s), pm);

        // Low prob - ignnore it.
        cond_sh_st = _mm512_kand(cond_sh, _mm512_knot(cond_st));
        if (!cond_sh_st.is_empty())
        {
            u = _mm512_mask_mov_ps(u, cond_sh_st, _mm512_mul_ps(riemann::g5, _mm512_fmadd_ps(riemann::g7, u, c)));
            uc = _mm512_div_ps(u, c);
            d = _mm512_mask_mov_ps(d, cond_sh_st, _mm512_mul_ps(d, _mm512_pow_ps(uc, riemann::g4)));
            p = _mm512_mask_mov_ps(p, cond_sh_st, _mm512_mul_ps(p, _mm512_pow_ps(uc, riemann::g3)));
        }

        // Final store.
        u = _mm512_mask_sub_ps(u, cond_um, zero, u);
    }

    /// <summary>
    /// Case sample vector.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="dl_p">Input array.</param>
    /// <param name="ul_p">Input array.</param>
    /// <param name="vl_p">Input array.</param>
    /// <param name="wl_p">Input array.</param>
    /// <param name="pl_p">Input array.</param>
    /// <param name="cl_p">Input array.</param>
    /// <param name="dr_p">Input array.</param>
    /// <param name="ur_p">Input array.</param>
    /// <param name="vr_p">Input array.</param>
    /// <param name="wr_p">Input array.</param>
    /// <param name="pr_p">Input array.</param>
    /// <param name="cr_p">Input array.</param>
    /// <param name="pm_p">Input array.</param>
    /// <param name="um_p">Input array.</param>
    /// <param name="d_p">Output array.</param>
    /// <param name="u_p">Output array.</param>
    /// <param name="v_p">Output array.</param>
    /// <param name="w_p">Output array.</param>
    /// <param name="p_p">Output array.</param>
    void vcase_sample(int n,
                      float* dl_p,
                      float* ul_p,
                      float* vl_p,
                      float* wl_p,
                      float* pl_p,
                      float* cl_p,
                      float* dr_p,
                      float* ur_p,
                      float* vr_p,
                      float* wr_p,
                      float* pr_p,
                      float* cr_p,
                      float* pm_p,
                      float* um_p,
                      float* d_p,
                      float* u_p,
                      float* v_p,
                      float* w_p,
                      float* p_p)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();
            ZMM d, u, v, w, p;

            vcase_sample_1(_mm512_load_ps(dl_p + sh),
                           _mm512_load_ps(ul_p + sh),
                           _mm512_load_ps(vl_p + sh),
                           _mm512_load_ps(wl_p + sh),
                           _mm512_load_ps(pl_p + sh),
                           _mm512_load_ps(cl_p + sh),
                           _mm512_load_ps(dr_p + sh),
                           _mm512_load_ps(ur_p + sh),
                           _mm512_load_ps(vr_p + sh),
                           _mm512_load_ps(wr_p + sh),
                           _mm512_load_ps(pr_p + sh),
                           _mm512_load_ps(cr_p + sh),
                           _mm512_load_ps(pm_p + sh),
                           _mm512_load_ps(um_p + sh),
                           d, u, v, w, p);

            _mm512_store_ps(d_p + sh, d);
            _mm512_store_ps(u_p + sh, u);
            _mm512_store_ps(v_p + sh, v);
            _mm512_store_ps(w_p + sh, w);
            _mm512_store_ps(p_p + sh, p);
        }
    }

    bool case_sample(int len,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * ZMM::count<float>();

        ArrayManager<float> dl(n);
        ArrayManager<float> ul(n);
        ArrayManager<float> vl(n);
        ArrayManager<float> wl(n);
        ArrayManager<float> pl(n);
        ArrayManager<float> cl(n);
        ArrayManager<float> dr(n);
        ArrayManager<float> ur(n);
        ArrayManager<float> vr(n);
        ArrayManager<float> wr(n);
        ArrayManager<float> pr(n);
        ArrayManager<float> cr(n);
        ArrayManager<float> pm(n);
        ArrayManager<float> um(n);
        ArrayManager<float> sd(n);
        ArrayManager<float> vd(n);
        ArrayManager<float> su(n);
        ArrayManager<float> vu(n);
        ArrayManager<float> sv(n);
        ArrayManager<float> vv(n);
        ArrayManager<float> sw(n);
        ArrayManager<float> vw(n);
        ArrayManager<float> sp(n);
        ArrayManager<float> vp(n);

        dl.generate_random(random_lo, random_hi);
        ul.generate_random(random_lo, random_hi);
        vl.generate_random(random_lo, random_hi);
        wl.generate_random(random_lo, random_hi);
        pl.generate_random(random_lo, random_hi);
        cl.generate_random(random_lo, random_hi);
        dr.generate_random(random_lo, random_hi);
        ur.generate_random(random_lo, random_hi);
        vr.generate_random(random_lo, random_hi);
        wr.generate_random(random_lo, random_hi);
        pr.generate_random(random_lo, random_hi);
        cr.generate_random(random_lo, random_hi);
        pm.generate_random(random_lo, random_hi);
        um.generate_random(random_lo, random_hi);

        scase_sample(n,
                     dl.get_data(), ul.get_data(), vl.get_data(), wl.get_data(), pl.get_data(), cl.get_data(),
                     dr.get_data(), ur.get_data(), vr.get_data(), wr.get_data(), pr.get_data(), cr.get_data(),
                     pm.get_data(), um.get_data(),
                     sd.get_data(), su.get_data(), sv.get_data(), sw.get_data(), sp.get_data());
                     
        vcase_sample(n,
                     dl.get_data(), ul.get_data(), vl.get_data(), wl.get_data(), pl.get_data(), cl.get_data(),
                     dr.get_data(), ur.get_data(), vr.get_data(), wr.get_data(), pr.get_data(), cr.get_data(),
                     pm.get_data(), um.get_data(),
                     vd.get_data(), vu.get_data(), vv.get_data(), vw.get_data(), vp.get_data());

        bool res = (vd.max_diff(sd) + vu.max_diff(su) + vv.max_diff(sv) + vw.max_diff(sw) + vp.max_diff(sp) < eps);

        if (!res)
        {
            std::cout << "max_diff : "
                      << vd.max_diff(sd) << ", " << vu.max_diff(su) << ", " << vv.max_diff(sv) << ", "
                      << vw.max_diff(sw) << ", " << vp.max_diff(sp) << std::endl;
        }

        return res;
    }

    // Riemann solver
    // starpu function

    /// <summary>
    /// Case starpu scalar.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="cl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="cr">Input.</param>
    /// <param name="p">Output.</param>
    /// <param name="u">Output.</param>
    void scase_starpu_1(float dl,
                        float ul,
                        float pl,
                        float cl,
                        float dr,
                        float ur,
                        float pr,
                        float cr,
                        float& p,
                        float& u)
    {
        const int32_t nriter = 20;
        const float tolpre = 1.0e-6f;
        float change, fl, fld, fr, frd, pold, pstart, udiff;

        // Guessed value pstart is computed.
        scase_guessp_1(dl, ul, pl, cl, dr, ur, pr, cr, pstart);
        pold = pstart;
        udiff = ur - ul;

        int32_t i = 1;

        for (; i <= nriter; i++)
        {
            scase_prefun_1(fl, fld, pold, dl, pl, cl);
            scase_prefun_1(fr, frd, pold, dr, pr, cr);
            p = pold - (fl + fr + udiff) / (fld + frd);
            change = 2.0f * abs((p - pold) / (p + pold));

            if (change <= tolpre)
            {
                break;
            }

            if (p < 0.0)
            {
                p = tolpre;
            }

            pold = p;
        }

        if (i > nriter)
        {
            p = 99.0;
        }

        // compute velocity in star region
        u = 0.5f * (ul + ur + fr - fl);
    }

    /// <summary>
    /// Case starpu scalar.
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
    /// <param name="p">Output array.</param>
    /// <param name="u">Output array.</param>
    void scase_starpu(int n,
                      float* dl,
                      float* ul,
                      float* pl,
                      float* cl,
                      float* dr,
                      float* ur,
                      float* pr,
                      float* cr,
                      float* p,
                      float* u)
    {
        for (int i = 0; i < n; i++)
        {
            scase_starpu_1(dl[i], ul[i], pl[i], cl[i], dr[i], ur[i], pr[i], cr[i], p[i], u[i]);
        }
    }

    /// <summary>
    /// Case starpu vector.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="cl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="cr">Input.</param>
    /// <param name="p">Output.</param>
    /// <param name="u">Output.</param>
    void vcase_starpu_1(ZMM dl,
                        ZMM ul,
                        ZMM pl,
                        ZMM cl,
                        ZMM dr,
                        ZMM ur,
                        ZMM pr,
                        ZMM cr,
                        ZMM& p,
                        ZMM& u)
    {
        ZMM tolpre, tolpre2, udiff, pold, fl, fld, fr, frd, change;
        Mask cond_break, cond_neg, m;
        const int nriter = 20;
        int iter = 1;

        tolpre = _mm512_set1_ps(1.0e-6f);
        tolpre2 = _mm512_set1_ps(5.0e-7f);
        udiff = _mm512_sub_ps(ur, ul);

        vcase_guessp_1(dl, ul, pl, cl, dr, ur, pr, cr, pold);

        // Start with full mask.
        m = Mask::full_tail(ZMM::count<float>());

        for (; (iter <= nriter) && (!m.is_empty_tail(ZMM::count<float>())); iter++)
        {
            vcase_prefun_1(fl, fld, pold, dl, pl, cl, m);
            vcase_prefun_1(fr, frd, pold, dr, pr, cr, m);
            p = _mm512_mask_sub_ps(p, m, pold,
                                   _mm512_mask_div_ps(zero, m,
                                                      _mm512_add_ps(_mm512_add_ps(fl, fr), udiff),
                                                      _mm512_add_ps(fld, frd)));
            change = _mm512_abs_ps(_mm512_mask_div_ps(zero, m, _mm512_sub_ps(p, pold), _mm512_add_ps(p, pold)));
            cond_break = _mm512_mask_cmple_ps_mask(m, change, tolpre2);
            m = _mm512_kand(m, _mm512_knot(cond_break));
            cond_neg = _mm512_mask_cmplt_ps_mask(m, p, zero);
            p = _mm512_mask_mov_ps(p, cond_neg, tolpre);
            pold = _mm512_mask_mov_ps(pold, m, p);
        }

        // Check for divergence.
        if (iter > nriter)
        {
            std::cout << "divergence in Newton-Raphson iteration" << std::endl;
            exit(1);
        }

        u = _mm512_mul_ps(half, _mm512_add_ps(_mm512_add_ps(ul, ur), _mm512_sub_ps(fr, fl)));
    }

    /// <summary>
    /// Case starpu vector.
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
    /// <param name="p_p">Output array.</param>
    /// <param name="u_p">Output array.</param>
    void vcase_starpu(int n,
                      float* dl_p,
                      float* ul_p,
                      float* pl_p,
                      float* cl_p,
                      float* dr_p,
                      float* ur_p,
                      float* pr_p,
                      float* cr_p,
                      float* p_p,
                      float* u_p)
    {
        assert(n % ZMM::count<float>() == 0);
        int vn = n / ZMM::count<float>();

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * ZMM::count<float>();
            ZMM p, u;

            vcase_starpu_1(_mm512_load_ps(dl_p + sh),
                           _mm512_load_ps(ul_p + sh),
                           _mm512_load_ps(pl_p + sh),
                           _mm512_load_ps(cl_p + sh),
                           _mm512_load_ps(dr_p + sh),
                           _mm512_load_ps(ur_p + sh),
                           _mm512_load_ps(pr_p + sh),
                           _mm512_load_ps(cr_p + sh),
                           p, u);

            _mm512_store_ps(p_p + sh, p);
            _mm512_store_ps(u_p + sh, u);
        }
    }

    bool case_starpu(int len,
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
        ArrayManager<float> sp(n);
        ArrayManager<float> vp(n);
        ArrayManager<float> su(n);
        ArrayManager<float> vu(n);

        dl.generate_random(random_lo, random_hi);
        ul.generate_random(random_lo, random_hi);
        pl.generate_random(random_lo, random_hi);
        cl.generate_random(random_lo, random_hi);
        dr.generate_random(random_lo, random_hi);
        ur.generate_random(random_lo, random_hi);
        pr.generate_random(random_lo, random_hi);
        cr.generate_random(random_lo, random_hi);
        sp.generate_random(random_lo, random_hi);
        vp.generate_random(random_lo, random_hi);
        su.generate_random(random_lo, random_hi);
        vu.generate_random(random_lo, random_hi);

        scase_starpu(n,
                     dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                     dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                     sp.get_data(), su.get_data());

        vcase_starpu(n,
                     dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                     dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                     vp.get_data(), vu.get_data());

        bool res = (vp.max_diff(sp) + vu.max_diff(su) < eps);

        if (!res)
        {
            std::cout << "max_diff : " << vp.max_diff(sp) << ", " << vu.max_diff(su) << std::endl;
        }

        return res;
    }
}
