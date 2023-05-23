#include "stdafx.h"

#ifndef CASES_H
#define CASES_H

#include "fv.h"

namespace fv
{
    // arith_f32

    void scase_arith_f32(int n,
                         float* a,
                         float* b,
                         float* c);

    void vcase_arith_f32(int n,
                         float* a,
                         float* b,
                         float* c);

    bool case_arith_f32(int len,
                        float random_lo = 10.0,
                        float random_hi = 100.0);

    // blend_f32

    void scase_blend_f32(int n,
                         float* a,
                         float* b,
                         float* c);

    void vcase_blend_f32(int n,
                         float* a,
                         float* b,
                         float* c);

    bool case_blend_f32(int len,
                        float random_lo = 10.0,
                        float random_hi = 100.0);

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
                      float* pm);

    void vcase_guessp(int n,
                      float* dl,
                      float* ul,
                      float* pl,
                      float* cl,
                      float* dr,
                      float* ur,
                      float* pr,
                      float* cr,
                      float* pm);

    bool case_guessp(int len,
                     float random_lo = 1.0,
                     float random_hi = 2.0,
                     float eps = 1.0e-4);

    // Riemann solver
    // prefun function

    void scase_prefun(int n,
                      float* f,
                      float* fd,
                      float* p,
                      float* dk,
                      float* pk,
                      float* ck);

    void vcase_prefun(int n,
                      float* f,
                      float* fd,
                      float* p,
                      float* dk,
                      float* pk,
                      float* ck);

    bool case_prefun(int len,
                     float random_lo = 1.0,
                     float random_hi = 2.0,
                     float eps = 1.0e-6);

    // Riemann solver
    // sample function

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
                      float* p);

    void vcase_sample(int n,
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
                      float* p);

    bool case_sample(int len,
                     float random_lo = 1.0,
                     float random_hi = 2.0,
                     float eps = 1.0e-6);

    // Riemann solver
    // starpu function

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
                      float* u);

    void vcase_starpu(int n,
                      float* dl,
                      float* ul,
                      float* pl,
                      float* cl,
                      float* dr,
                      float* ur,
                      float* pr,
                      float* cr,
                      float* p,
                      float* u);

    bool case_starpu(int len,
                     float random_lo = 1.0,
                     float random_hi = 2.0,
                     float eps = 1.0e-6);
}

#endif
