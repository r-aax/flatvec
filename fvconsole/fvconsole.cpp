// fvconsole.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "fvcases.h"

using namespace fv;

/// <summary>
/// Main function.
/// </summary>
/// <returns>Return code.</returns>
int main()
{
    // Cases for detecting some inefficiency in the vector code.

    // Trivial cases.
    test_case("arith_f32", [] (int n, int r) { return case_arith_f32(n, r); },       1, 1, false);
    test_case("blend_f32", [] (int n, int r) { return case_blend_f32(n, r); },       1, 1, false);
    test_case("square_eq", [] (int n, int r) { return case_square_equation(n, r); }, 1, 1, false);

    // Cases from real applications.
    test_case("guessp   ", [] (int n, int r) { return case_guessp(n, r); },          1, 1, false);
    test_case("prefun   ", [] (int n, int r) { return case_prefun(n, r); },          1, 1, false);
    test_case("sample   ", [] (int n, int r) { return case_sample(n, r); },          1, 1, false);
    test_case("starpu   ", [] (int n, int r) { return case_starpu(n, r); },          1, 1, false);
#ifndef LINUX_GCC_BUILD
    test_case("riemann  ", [] (int n, int r) { return case_riemann(n, r); },         1, 1, false);
#endif

    return 0;
}
