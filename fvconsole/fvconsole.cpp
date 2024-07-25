// fvconsole.cpp : Defines the entry point for the console application.

#include "stdafx.h"

#include "fvcases.h"

using namespace fv;

/// <summary>
/// Main function.
/// </summary>
/// <returns>Return code.</returns>
int main()
{

#if 0
    // Cases for detecting some inefficiency in the vector code.
    test_case("warn_hang",
              [] (int n, int r) { return case_warn_hang(n, r); },
              1, 1, false);
    test_case("warn_multi_rewrite",
              [] (int n, int r) { return case_warn_multi_rewrite(n, r); },
              1, 1, false);
#endif

#if 0
    // Trivial cases.
    test_case("arith_f32",
              [] (int n, int r) { return case_arith_f32(n, r); },
              1, 1, false);
    test_case("blend_f32",
              [] (int n, int r) { return case_blend_f32(n, r); },
              1, 1, false);
    test_case("square_eq",
              [] (int n, int r) { return case_square_equation(n, r); },
              1, 1, false);
#endif

    //
    // Cases from real applications.
    //

#if 0
    test_case("guessp",
              [] (int n, int r) { return case_guessp(n, r); },
              1, 1, false);
#endif

#if 1
    test_case("prefun",
              [] (int n, int r) { return case_prefun(n, r); },
              /*141186*/ 70000, 1, false);
#endif

#if 0
    test_case("sample",
              [] (int n, int r) { return case_sample(n, r); },
              1, 1, false);
#endif

#if 0
    test_case("starpu",
              [] (int n, int r) { return case_starpu(n, r); },
              1, 1, false);
#endif

#if 0
#ifndef LINUX_GCC_BUILD
    test_case("riemann",
              [] (int n, int r) { return case_riemann(n, r); },
              26249, 1, false);
#endif
#endif

    return 0;
}
