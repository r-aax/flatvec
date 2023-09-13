#include "stdafx.h"

#include "cases.h"

#include "array_manager.h"
#include "global_stat.h"

#ifndef LINUX_ICC_BUILD

#include "fv.h"

// Redefenition vector types with their implementation.
using _m512 = fv::ZMM;
using __mmask16 = fv::Mask;

#endif

#ifdef LINUX_ICC_BUILD
#define CNT_FLOAT 16
#else
#define CNT_FLOAT ZMM::count<float>()
#endif

namespace fv
{
    // Some constants

    /// <summary>
    /// Zero.
    /// </summary>
    _m512 zero = _mm512_set1_ps(0.0f);

    /// <summary>
    /// One.
    /// </summary>
    _m512 one = _mm512_set1_ps(1.0f);

    /// <summary>
    /// Two.
    /// </summary>
    _m512 two = _mm512_set1_ps(2.0f);

    /// <summary>
    /// Four.
    /// </summary>
    _m512 four = _mm512_set1_ps(4.0f);

    /// <summary>
    /// Half.
    /// </summary>
    _m512 half = _mm512_set1_ps(0.5f);

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
    void vcase_arith_f32_1(_m512 a,
                           _m512 b,
                           _m512& c)
    {
        _m512 add = _mm512_add_ps(a, b);
        _m512 sub = _mm512_sub_ps(a, b);
        _m512 mul = _mm512_mul_ps(a, b);
        _m512 div = _mm512_div_ps(a, b);
        _m512 min = _mm512_min_ps(a, b);
        _m512 max = _mm512_max_ps(a, b);

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
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 c;

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
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_arith_f32(int len,
                        int repeats,
                        float random_lo,
                        float random_hi)
    {
        int n = len * CNT_FLOAT;

        ArrayManager<float> a(n);
        ArrayManager<float> b(n);
        ArrayManager<float> sc(n);
        ArrayManager<float> vc(n);

        a.generate_random(random_lo, random_hi);
        b.generate_random(random_lo, random_hi);

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_arith_f32(n, a.get_data(), b.get_data(), sc.get_data());
        }
        
        GS.fix_time_middle();
        
        for (int i = 0; i < repeats; i++)
        {
            vcase_arith_f32(n, a.get_data(), b.get_data(), vc.get_data());
        }
        
        GS.fix_time_after();

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
    void vcase_blend_f32_1(_m512 a,
                           _m512 b,
                           _m512& c)
    {
        __mmask16 k = _mm512_cmple_ps_mask(b, a);
        _m512 add, mul;

        add = _mm512_add_ps(a, b);
        mul = _mm512_mul_ps(a, b);
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
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 c;

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
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_blend_f32(int len,
                        int repeats,
                        float random_lo,
                        float random_hi)
    {
        int n = len * CNT_FLOAT;

        ArrayManager<float> a(n);
        ArrayManager<float> b(n);
        ArrayManager<float> sc(n);
        ArrayManager<float> vc(n);

        a.generate_random(random_lo, random_hi);
        b.generate_random(random_lo, random_hi);

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_blend_f32(n, a.get_data(), b.get_data(), sc.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_blend_f32(n, a.get_data(), b.get_data(), vc.get_data());
        }

        GS.fix_time_after();

        bool res = (vc.max_diff(sc) == 0.0);

        if (!res)
        {
            std::cout << "max_diff : " << vc.max_diff(sc) << std::endl;
        }

        return res;
    }

    // Constants for Riemann solver

    /// <summary>
    /// Namespace for riemann cases.
    /// </summary>
    namespace riemann
    {
        /// <summary>
        /// Constant G.
        /// </summary>
        float sg {1.4f};

        /// <summary>
        /// Constant G1.
        /// </summary>
        float sg1 = (sg - 1.0f) / (2.0f * sg);

        /// <summary>
        /// Constant G2.
        /// </summary>
        float sg2 = (sg + 1.0f) / (2.0f * sg);

        /// <summary>
        /// Constant G3.
        /// </summary>
        float sg3 = 2.0f * sg / (sg - 1.0f);

        /// <summary>
        /// Constant G4.
        /// </summary>
        float sg4 = 2.0f / (sg - 1.0f);

        /// <summary>
        /// Constant G5.
        /// </summary>
        float sg5 = 2.0f / (sg + 1.0f);

        /// <summary>
        /// Constant G6.
        /// </summary>
        float sg6 = (sg - 1.0f) / (sg + 1.0f);

        /// <summary>
        /// Constant G7.
        /// </summary>
        float sg7 = (sg - 1.0f) / 2.0f;

        /// <summary>
        /// Constant G8.
        /// </summary>
        float sg8 = (sg - 1.0f);

        /// <summary>
        /// Vector constant G.
        /// </summary>
        _m512 g = _mm512_set1_ps(sg);

        /// <summary>
        /// Vector constant G1.
        /// </summary>
        _m512 g1 = _mm512_set1_ps(sg1);

        /// <summary>
        /// Vector constant G2.
        /// </summary>
        _m512 g2 = _mm512_set1_ps(sg2);

        /// <summary>
        /// Vector constant G3.
        /// </summary>
        _m512 g3 = _mm512_set1_ps(sg3);

        /// <summary>
        /// Vector constant G4.
        /// </summary>
        _m512 g4 = _mm512_set1_ps(sg4);

        /// <summary>
        /// Vector constant G5.
        /// </summary>
        _m512 g5 = _mm512_set1_ps(sg5);

        /// <summary>
        /// Vector constant G6.
        /// </summary>
        _m512 g6 = _mm512_set1_ps(sg6);

        /// <summary>
        /// Vector constant G7.
        /// </summary>
        _m512 g7 = _mm512_set1_ps(sg7);

        /// <summary>
        /// Vector constant G8.
        /// </summary>
        _m512 g8 = _mm512_set1_ps(sg8);
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
    void vcase_guessp_1(_m512 dl,
                        _m512 ul,
                        _m512 pl,
                        _m512 cl,
                        _m512 dr,
                        _m512 ur,
                        _m512 pr,
                        _m512 cr,
                        _m512& pm)
    {
        // Begin of calculation part.

        _m512 cup, ppv, pmin, pmax, qmax, pq, um, ptl, ptr, gel, ger, pqcr;
        __mmask16 cond_pvrs, cond_ppv, ncond_ppv;

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
	    if (cond_ppv != 0x0)
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
	    if (ncond_ppv != 0x0)
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
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 pm;

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
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <param name="eps">Max deviation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_guessp(int len,
                     int repeats,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * CNT_FLOAT;

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

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_guessp(n,
                         dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                         dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                         spm.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_guessp(n,
                         dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                         dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                         vpm.get_data());
        }

        GS.fix_time_after();

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
    void vcase_prefun_1(_m512& f,
                        _m512& fd,
                        _m512 p,
                        _m512 dk,
                        _m512 pk,
                        _m512 ck,
                        __mmask16 m)
    {
        _m512 pratio, ak, bkp, ppk, qrt;
        __mmask16 cond, ncond;

        // Conditions.
        cond = _mm512_kand(_mm512_cmple_ps_mask(p, pk), m);
        ncond = _mm512_kand(_mm512_knot(cond), m);

        // The first branch.
	    if (cond != 0x0)
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
	    if (ncond != 0x0)
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
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 f, fd;

            vcase_prefun_1(f, fd,
                           _mm512_load_ps(p_p + sh),
                           _mm512_load_ps(dk_p + sh),
                           _mm512_load_ps(pk_p + sh),
                           _mm512_load_ps(ck_p + sh),
                           0xffff);

            _mm512_store_ps(f_p + sh, f);
            _mm512_store_ps(fd_p + sh, fd);
        }
    }

    /// <summary>
    /// Case prefun.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Lo value for random generation.</param>
    /// <param name="random_hi">Hi value for random generation.</param>
    /// <param name="eps">Max deviation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_prefun(int len,
                     int repeats,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * CNT_FLOAT;

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

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_prefun(n,
                         sf.get_data(), sfd.get_data(), p.get_data(),
                         dk.get_data(), pk.get_data(), ck.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_prefun(n,
                         vf.get_data(), vfd.get_data(), p.get_data(),
                         dk.get_data(), pk.get_data(), ck.get_data());
        }

        GS.fix_time_after();

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
    void vcase_sample_1(_m512 dl,
                        _m512 ul,
                        _m512 vl,
                        _m512 wl,
                        _m512 pl,
                        _m512 cl,
                        _m512 dr,
                        _m512 ur,
                        _m512 vr,
                        _m512 wr,
                        _m512 pr,
                        _m512 cr,
                        _m512 pm,
                        _m512 um,
                        _m512& d,
                        _m512& u,
                        _m512& v,
                        _m512& w,
                        _m512& p)
    {
        _m512 c, ums, pms, sh, st, s, uc;
        __mmask16 cond_um, cond_pm, cond_sh, cond_st, cond_s, cond_sh_st;

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
	    if (cond_sh_st != 0x0)
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
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 d, u, v, w, p;

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

    /// <summary>
    /// Case sample.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Lo bound for random.</param>
    /// <param name="random_hi">Hi bound for random.</param>
    /// <param name="eps">Max deviation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_sample(int len,
                     int repeats,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * CNT_FLOAT;

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

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_sample(n,
                         dl.get_data(), ul.get_data(), vl.get_data(), wl.get_data(), pl.get_data(), cl.get_data(),
                         dr.get_data(), ur.get_data(), vr.get_data(), wr.get_data(), pr.get_data(), cr.get_data(),
                         pm.get_data(), um.get_data(),
                         sd.get_data(), su.get_data(), sv.get_data(), sw.get_data(), sp.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_sample(n,
                         dl.get_data(), ul.get_data(), vl.get_data(), wl.get_data(), pl.get_data(), cl.get_data(),
                         dr.get_data(), ur.get_data(), vr.get_data(), wr.get_data(), pr.get_data(), cr.get_data(),
                         pm.get_data(), um.get_data(),
                         vd.get_data(), vu.get_data(), vv.get_data(), vw.get_data(), vp.get_data());
        }

        GS.fix_time_after();

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
    void vcase_starpu_1(_m512 dl,
                        _m512 ul,
                        _m512 pl,
                        _m512 cl,
                        _m512 dr,
                        _m512 ur,
                        _m512 pr,
                        _m512 cr,
                        _m512& p,
                        _m512& u)
    {
        _m512 tolpre, tolpre2, udiff, pold, fl, fld, fr, frd, change;
        __mmask16 cond_break, cond_neg, m;
        const int nriter = 20;
        int iter = 1;

        tolpre = _mm512_set1_ps(1.0e-6f);
        tolpre2 = _mm512_set1_ps(5.0e-7f);
        udiff = _mm512_sub_ps(ur, ul);

        vcase_guessp_1(dl, ul, pl, cl, dr, ur, pr, cr, pold);

        // Start with full mask.
	    m = 0xffff;

	    for (; (iter <= nriter) && (m != 0x0); iter++)
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
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 p, u;

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

    /// <summary>
    /// Case starpu.
    /// </summary>
    /// <param name="len">Count of vectors.</param>
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Lo bound for random.</param>
    /// <param name="random_hi">Hi bound for random.</param>
    /// <param name="eps">Max deviation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_starpu(int len,
                     int repeats,
                     float random_lo,
                     float random_hi,
                     float eps)
    {
        int n = len * CNT_FLOAT;

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

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_starpu(n,
                         dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                         dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                         sp.get_data(), su.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_starpu(n,
                         dl.get_data(), ul.get_data(), pl.get_data(), cl.get_data(),
                         dr.get_data(), ur.get_data(), pr.get_data(), cr.get_data(),
                         vp.get_data(), vu.get_data());
        }

        GS.fix_time_after();

        bool res = (vp.max_diff(sp) + vu.max_diff(su) < eps);

        if (!res)
        {
            std::cout << "max_diff : " << vp.max_diff(sp) << ", " << vu.max_diff(su) << std::endl;
        }

        return res;
    }

    /// <summary>
    /// Scalar riemann case.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="vl">Input.</param>
    /// <param name="wl">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="vr">Input.</param>
    /// <param name="wr">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="d">Output.</param>
    /// <param name="u">Output.</param>
    /// <param name="v">Output.</param>
    /// <param name="w">Output.</param>
    /// <param name="p">Output.</param>
    void scase_riemann_1(float dl,
                         float ul,
                         float vl,
                         float wl,
                         float pl,
                         float dr,
                         float ur,
                         float vr,
                         float wr,
                         float pr,
                         float& d,
                         float& u,
                         float& v,
                         float& w,
                         float& p)
    {
        float pm, um, cl, cr;

        pm = 0.0;

        // Sound speeds.
        cl = sqrt(riemann::sg * pl / dl);
        cr = sqrt(riemann::sg * pr / dr);

        // Check for vacuum.
        if (riemann::sg4 * (cl + cr) <= (ur - ul))
        {
            std::cout << "VACUUM" << std::endl;
            exit(1);
        }

        // Exact solution.
        scase_starpu_1(dl, ul, pl, cl, dr, ur, pr, cr, pm, um);
        scase_sample_1(dl, ul, vl, wl, pl, cl,
                       dr, ur, vr, wr, pr, cr,
                       pm, um,
                       d, u, v, w, p);
    }

    /// <summary>
    /// Scalar case riemann.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="dl">Input array.</param>
    /// <param name="ul">Input array.</param>
    /// <param name="vl">Input array.</param>
    /// <param name="wl">Input array.</param>
    /// <param name="pl">Input array.</param>
    /// <param name="dr">Input array.</param>
    /// <param name="ur">Input array.</param>
    /// <param name="vr">Input array.</param>
    /// <param name="wr">Input array.</param>
    /// <param name="pr">Input array.</param>
    /// <param name="d">Output array.</param>
    /// <param name="u">Output array.</param>
    /// <param name="v">Output array.</param>
    /// <param name="w">Output array.</param>
    /// <param name="p">Output array.</param>
    void scase_riemann(int n,
                       float* dl,
                       float* ul,
                       float* vl,
                       float* wl,
                       float* pl,
                       float* dr,
                       float* ur,
                       float* vr,
                       float* wr,
                       float* pr,
                       float* d,
                       float* u,
                       float* v,
                       float* w,
                       float* p)
    {
        for (int i = 0; i < n; i++)
        {
            scase_riemann_1(dl[i], ul[i], vl[i], wl[i], pl[i],
                            dr[i], ur[i], vr[i], wr[i], pr[i],
                            d[i], u[i], v[i], w[i], p[i]);
        }
    }

    /// <summary>
    /// Vector riemann case.
    /// </summary>
    /// <param name="dl">Input.</param>
    /// <param name="ul">Input.</param>
    /// <param name="vl">Input.</param>
    /// <param name="wl">Input.</param>
    /// <param name="pl">Input.</param>
    /// <param name="dr">Input.</param>
    /// <param name="ur">Input.</param>
    /// <param name="vr">Input.</param>
    /// <param name="wr">Input.</param>
    /// <param name="pr">Input.</param>
    /// <param name="d">Output.</param>
    /// <param name="u">Output.</param>
    /// <param name="v">Output.</param>
    /// <param name="w">Output.</param>
    /// <param name="p">Output.</param>
    void vcase_riemann_1(_m512 dl,
                         _m512 ul,
                         _m512 vl,
                         _m512 wl,
                         _m512 pl,
                         _m512 dr,
                         _m512 ur,
                         _m512 vr,
                         _m512 wr,
                         _m512 pr,
                         _m512& d,
                         _m512& u,
                         _m512& v,
                         _m512& w,
                         _m512& p)
    {
        _m512 cl, cr, pm, um;
        __mmask16 vacuum_mask;

        cl = _mm512_sqrt_ps(_mm512_div_ps(_mm512_mul_ps(riemann::g, pl), dl));
        cr = _mm512_sqrt_ps(_mm512_div_ps(_mm512_mul_ps(riemann::g, pr), dr));
        vacuum_mask = _mm512_cmple_ps_mask(_mm512_mul_ps(riemann::g4, _mm512_add_ps(cl, cr)),
                                           _mm512_sub_ps(ur, ul));

        // Vacuum check.
	    if (vacuum_mask != 0x0)
        {
            std::cout << "VACUUM" << std::endl;
            exit(1);
        }
        vcase_starpu_1(dl, ul, pl, cl, dr, ur, pr, cr, pm, um);
        vcase_sample_1(dl, ul, vl, wl, pl, cl,
                       dr, ur, vr, wr, pr, cr,
                       pm, um,
                       d, u, v, w, p);
    }

    /// <summary>
    /// Vectorized riemann case.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="dl_p">Input array.</param>
    /// <param name="ul_p">Input array.</param>
    /// <param name="vl_p">Input array.</param>
    /// <param name="wl_p">Input array.</param>
    /// <param name="pl_p">Input array.</param>
    /// <param name="dr_p">Input array.</param>
    /// <param name="ur_p">Input array.</param>
    /// <param name="vr_p">Input array.</param>
    /// <param name="wr_p">Input array.</param>
    /// <param name="pr_p">Input array.</param>
    /// <param name="d_p">Output array.</param>
    /// <param name="u_p">Output array.</param>
    /// <param name="v_p">Output array.</param>
    /// <param name="w_p">Output array.</param>
    /// <param name="p_p">Output array.</param>
    void vcase_riemann(int n,
                       float* dl_p,
                       float* ul_p,
                       float* vl_p,
                       float* wl_p,
                       float* pl_p,
                       float* dr_p,
                       float* ur_p,
                       float* vr_p,
                       float* wr_p,
                       float* pr_p,
                       float* d_p,
                       float* u_p,
                       float* v_p,
                       float* w_p,
                       float* p_p)
    {
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 d, u, v, w, p;

            vcase_riemann_1(_mm512_load_ps(dl_p + sh),
                            _mm512_load_ps(ul_p + sh),
                            _mm512_load_ps(vl_p + sh),
                            _mm512_load_ps(wl_p + sh),
                            _mm512_load_ps(pl_p + sh),
                            _mm512_load_ps(dr_p + sh),
                            _mm512_load_ps(ur_p + sh),
                            _mm512_load_ps(vr_p + sh),
                            _mm512_load_ps(wr_p + sh),
                            _mm512_load_ps(pr_p + sh),
                            d, u, v, w, p);

            _mm512_store_ps(d_p + sh, d);
            _mm512_store_ps(u_p + sh, u);
            _mm512_store_ps(v_p + sh, v);
            _mm512_store_ps(w_p + sh, w);
            _mm512_store_ps(p_p + sh, p);
        }
    }

    /// <summary>
    /// Case riemann.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="repeats">Repeats count.</param>
    /// <param name="eps">Max deviation.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_riemann(int len,
                      int repeats,
                      float eps)
    {
        int n = len * CNT_FLOAT;

        // Files ../fvcases/data/case_riemann_{dl,ul,pl,dr,ur,pr}.txt
        // contain 419997 float elements.
        // So it has no sense to pass len more than 26249.

        ArrayManager<float> dl(n, "../fvcases/data/case_riemann_dl.txt");
        ArrayManager<float> ul(n, "../fvcases/data/case_riemann_ul.txt");
        ArrayManager<float> vl(n);
        ArrayManager<float> wl(n);
        ArrayManager<float> pl(n, "../fvcases/data/case_riemann_pl.txt");
        ArrayManager<float> dr(n, "../fvcases/data/case_riemann_dr.txt");
        ArrayManager<float> ur(n, "../fvcases/data/case_riemann_ur.txt");
        ArrayManager<float> vr(n);
        ArrayManager<float> wr(n);
        ArrayManager<float> pr(n, "../fvcases/data/case_riemann_pr.txt");
        ArrayManager<float> sd(n);
        ArrayManager<float> su(n);
        ArrayManager<float> sv(n);
        ArrayManager<float> sw(n);
        ArrayManager<float> sp(n);
        ArrayManager<float> vd(n);
        ArrayManager<float> vu(n);
        ArrayManager<float> vv(n);
        ArrayManager<float> vw(n);
        ArrayManager<float> vp(n);

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_riemann(n,
                          dl.get_data(), ul.get_data(), vl.get_data(), wl.get_data(), pl.get_data(),
                          dr.get_data(), ur.get_data(), vr.get_data(), wr.get_data(), pr.get_data(),
                          sd.get_data(), su.get_data(), sv.get_data(), sw.get_data(), sp.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_riemann(n,
                          dl.get_data(), ul.get_data(), vl.get_data(), wl.get_data(), pl.get_data(),
                          dr.get_data(), ur.get_data(), vr.get_data(), wr.get_data(), pr.get_data(),
                          vd.get_data(), vu.get_data(), vv.get_data(), vw.get_data(), vp.get_data());
        }

        GS.fix_time_after();

        bool res = (vd.max_diff(sd) + vu.max_diff(su) + vv.max_diff(sv) + vw.max_diff(sw) + vp.max_diff(sp) < eps);

        if (!res)
        {
            std::cout << "max_diff : "
                      << vd.max_diff(sd) << ", " << vu.max_diff(su) << ", " << vv.max_diff(sv) << ", "
                      << vw.max_diff(sw) << ", " << vp.max_diff(sp)
                      << std::endl;
        }

        return res;
    }

    /// <summary>
    /// Case square equation scalar.
    /// </summary>
    /// <param name="a">Input.</param>
    /// <param name="b">Input.</param>
    /// <param name="c">Input.</param>
    /// <param name="h">Output.</param>
    void scase_square_equation_1(float a,
                                 float b,
                                 float c,
                                 float& h)
    {
        if (a == 0.0f)
        {
            if (b == 0.0f)
            {
                h = 0.0f;
            }
            else
            {
                h = -c / b;

                if (h < 0.0f);
                {
                    h = 0.0f;
                }
            }
        }
        else
        {
            b /= a;
            c /= a;

            float d = b * b - 4.0f * c;

            if (d < 0.0f)
            {
                h = 0.0f;
            }
            else
            {
                float sd = sqrt(d);
                float h1 = 0.5 * (-b - sd);

                if (h1 > 0.0f)
                {
                    h = h1;
                }
                else
                {
                    float h2 = 0.5 * (-b + sd);

                    if (h2 > 0.0f)
                    {
                        h = h2;
                    }
                    else
                    {
                        h = 0.0f;
                    }
                }
            }
        }
    }

    /// <summary>
    /// Case square equation scalar.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="a">Input array.</param>
    /// <param name="b">Input array.</param>
    /// <param name="c">Input array.</param>
    /// <param name="h">Output array.</param>
    void scase_square_equation(int n,
                               float* a,
                               float* b,
                               float* c,
                               float* h)
    {
        for (int i = 0; i < n; i++)
        {
            scase_square_equation_1(a[i], b[i], c[i], h[i]);
        }
    }

    /// <summary>
    /// Case square equation vector.
    /// </summary>
    /// <param name="a">Input.</param>
    /// <param name="b">Input.</param>
    /// <param name="c">Input.</param>
    /// <param name="h">Output.</param>
    void vcase_square_equation_1(_m512 a,
                                 _m512 b,
                                 _m512 c,
                                 _m512& h)
    {
        h = zero;
        __mmask16 p1 = _mm512_cmpeq_ps_mask(a, zero);
        __mmask16 p2 = _mm512_mask_cmpneq_ps_mask(p1, b, zero);
        h = _mm512_maskz_max_ps(p2, _mm512_maskz_div_ps(p2, _mm512_sub_ps(zero, c), b), zero);
        __mmask16 np1 = _mm512_knot(p1);
        b = _mm512_mask_div_ps(b, np1, b, a);
        c = _mm512_mask_div_ps(c, np1, c, a);
        _m512 d = _mm512_sub_ps(_mm512_mul_ps(b, b), _mm512_mul_ps(four, c));
        __mmask16 p4 = _mm512_mask_cmpge_ps_mask(np1, d, zero);
        _m512 sd = _mm512_mask_sqrt_ps(sd, p4, d);
        _m512 h1 = _mm512_mul_ps(half, _mm512_sub_ps(zero, _mm512_add_ps(b, sd)));
        h = _mm512_mask_mov_ps(h, p4, _mm512_mask_blend_ps(_mm512_cmpgt_ps_mask(h1, zero),
                                                           _mm512_max_ps(_mm512_mul_ps(half, _mm512_sub_ps(sd, b)), zero), h1));
    }

    /// <summary>
    /// Case for square equation.
    /// </summary>
    /// <param name="n">Count.</param>
    /// <param name="a_p">Input array.</param>
    /// <param name="b_p">Input array.</param>
    /// <param name="c_p">Input array.</param>
    /// <param name="h_p">Output array.</param>
    void vcase_square_equation(int n,
                               float* a_p,
                               float* b_p,
                               float* c_p,
                               float* h_p)
    {
        assert(n % CNT_FLOAT == 0);
        int vn = n / CNT_FLOAT;

        for (int vi = 0; vi < vn; vi++)
        {
            int sh = vi * CNT_FLOAT;
            _m512 h;

            vcase_square_equation_1(_mm512_load_ps(a_p + sh),
                                    _mm512_load_ps(b_p + sh),
                                    _mm512_load_ps(c_p + sh),
                                    h);

            _mm512_store_ps(h_p + sh, h);
        }
    }

    /// <summary>
    /// Find minimal positive root of square equation.
    /// </summary>
    /// <param name="len">Vectors count.</param>
    /// <param name="repeats">Repeats count.</param>
    /// <param name="random_lo">Low limit for random.</param>
    /// <param name="random_hi">High limit for random.</param>
    /// <returns>
    /// true - OK result,
    /// false - ERROR result.
    /// </returns>
    bool case_square_equation(int len,
                              int repeats,
                              float random_lo,
                              float random_hi)
    {
        int n = len * CNT_FLOAT;

        ArrayManager<float> a(n);
        ArrayManager<float> b(n);
        ArrayManager<float> c(n);
        ArrayManager<float> sh(n);
        ArrayManager<float> vh(n);

        a.generate_random(random_lo, random_hi);
        b.generate_random(random_lo, random_hi);
        c.generate_random(random_lo, random_hi);

        GS.fix_time_before();

        for (int i = 0; i < repeats; i++)
        {
            scase_square_equation(n, a.get_data(), b.get_data(), c.get_data(), sh.get_data());
        }

        GS.fix_time_middle();

        for (int i = 0; i < repeats; i++)
        {
            vcase_square_equation(n, a.get_data(), b.get_data(), c.get_data(), vh.get_data());
        }

        GS.fix_time_after();

        bool res = (vh.max_diff(sh) == 0.0);

        if (!res)
        {
            std::cout << "max_diff : " << vh.max_diff(sh) << std::endl;
        }

        return res;
    }

    //

    /// <summary>
    /// General case.
    /// </summary>
    /// <param name="name">Name.</param>
    /// <param name="fun">Function.</param>
    /// <param name="count">Count.</param>
    /// <param name="repeats">Repeats count.</param>
    void test_case(std::string name,
                   std::function<bool(int, int)> fun,
                   int count,
                   int repeats)
    {
        GS.clean();
        std::cout << name << " : " << (fun(count, repeats) ? "OK" : "ERROR") << std::endl;
        GS.print();
    }
}
