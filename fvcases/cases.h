#ifndef CASES_H
#define CASES_H

#include "stdafx.h"

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

    bool case_guessp(int len);
}

#endif
