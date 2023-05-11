// fvconsole.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>

#include "fv.h"

int main()
{
    __m512 a, b, r;

    a.setF32(3, 5.0);
    a.setF32(6, 3.0);
    b.setF32(6, 4.0);
    b.setF32(10, 10.0);

    r = _mm512_add_ps(a, b);

    std::cout << a.getF32Representation() << std::endl;
    std::cout << b.getF32Representation() << std::endl;
    std::cout << r.getF32Representation() << std::endl;

    getchar();
    return 0;
}
