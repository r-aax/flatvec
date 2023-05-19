#ifndef CASES_H
#define CASES_H

#include "stdafx.h"

#include "fv.h"

namespace fv
{
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
}

#endif
