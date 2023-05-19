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

    a.generate_random(10.0, 100.0);
    b.generate_random(10.0, 100.0);

    a.print();
    b.print();

    scase_arith_f32(n, a.getData(), b.getData(), sc.getData());
    vcase_arith_f32(n, a.getData(), b.getData(), vc.getData());
    sc.print();
    vc.print();
    std::cout << vc.maxDiff(sc);
    std::cout << " " << std::endl;

    if (vc.maxDiff(sc) > 0.0)
    {
        std::cout << "ERROR!" << std::endl;
    }

    getchar();
    return 0;
}
