#include "stdafx.h"

#ifndef CASES_H
#define CASES_H

namespace fv
{
    // arith_f32
    bool case_arith_f32(int len,
                        int repeats,
                        float random_lo = 10.0f,
                        float random_hi = 100.0f);

    // blend_f32
    bool case_blend_f32(int len,
                        int repeats,
                        float random_lo = 10.0f,
                        float random_hi = 100.0f);

    // Riemann solver
    // guessp function
    bool case_guessp(int len,
                     int repeats,
                     float random_lo = 1.0f,
                     float random_hi = 2.0f,
                     float eps = 1.0e-4f);

    // Riemann solver
    // prefun function
    bool case_prefun(int len,
                     int repeats,
                     float random_lo = 1.0f,
                     float random_hi = 2.0f,
                     float eps = 1.0e-6f);

    // Riemann solver
    // sample function
    bool case_sample(int len,
                     int repeats,
                     float random_lo = 2.0f,
                     float random_hi = 3.0f,
                     float eps = 1.0e-6f);

    // Riemann solver
    // starpu function
    bool case_starpu(int len,
                     int repeats,
                     float random_lo = 1.0f,
                     float random_hi = 1.1f,
                     float eps = 1.0e-6f);

    // Riemann solver
    // riemann function
    bool case_riemann(int len,
                      int repeats,
                      float eps = 1.0e-3f);

    //

    // Check general case.
    void test_case(std::string name,
                   std::function<bool(int, int)> fun,
                   int count,
                   int repeats);
}

#endif
