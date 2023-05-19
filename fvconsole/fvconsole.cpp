// fvconsole.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "fvcases.h"

using namespace fv;

int main()
{
    if (case_arith_f32(10))
    {
        std::cout << "OK" << std::endl;
    }
    else
    {
        std::cout << "ERR" << std::endl;
    }

    getchar();
    return 0;
}
