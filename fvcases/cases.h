#ifndef CASES_H
#define CASES_H

#include "stdafx.h"

#include "fv.h"

namespace fv
{
    void scase_add_f32(int n,
                       float* a,
                       float* b,
                       float* c);

    void vcase_add_f32(int n,
                       float* a,
                       float* b,
                       float* c);
}

#endif
