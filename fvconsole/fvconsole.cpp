// fvconsole.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>

#include "fv.h"

int main()
{
    fv::ZMM zmm;
    zmm.setI32(3, 100000);

    std::cout << zmm.getI32Representation() << std::endl;
    std::cout << zmm.getF32Representation() << std::endl;
    std::cout << zmm.getF64Representation() << std::endl;

    getchar();
    return 0;
}

