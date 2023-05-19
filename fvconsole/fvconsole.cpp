// fvconsole.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "fvcases.h"

using namespace fv;

int main()
{
    std::cout << (case_arith_f32(10) ? "arith_f32 : OK" : "arith_f32 : ERROR") << std::endl;
    std::cout << (case_blend_f32(10) ? "blend_f32 : OK" : "blend_f32 : ERROR") << std::endl;
    std::cout << (case_guessp(10)    ? "guessp    : OK" : "guessp    : ERROR") << std::endl;

    getchar();
    return 0;
}
