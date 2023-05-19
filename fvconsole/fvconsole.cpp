// fvconsole.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "fvcases.h"

using namespace fv;

int main()
{
    int n = 1 * ZMM::f32_count;
    ArrayManager<float> a(n);
    ArrayManager<float> b(n);
    ArrayManager<float> sc(n);
    ArrayManager<float> vc(n);

    a.generate_random(0.0, 100.0);
    b.generate_random(0.0, 100.0);
    a.print();
    b.print();

    scase_add_f32(n, a.getData(), b.getData(), sc.getData());
    vcase_add_f32(n, a.getData(), b.getData(), vc.getData());

    sc.print();
    vc.print();
    std::cout << vc.maxDiff(sc);

    getchar();
    return 0;
}
