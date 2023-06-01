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
    test_case("arith_f32", [] (int n) { return case_arith_f32(n); }, 10);
    test_case("blend_f32", [] (int n) { return case_blend_f32(n); }, 10);
    test_case("guessp   ", [] (int n) { return case_guessp(n); }, 10);
    test_case("prefun   ", [] (int n) { return case_prefun(n); }, 10);
    test_case("sample   ", [] (int n) { return case_sample(n); }, 10);
    test_case("starpu   ", [] (int n) { return case_starpu(n); }, 10);
    test_case("riemann  ", [] (int n) { return case_riemann(n); }, 10);

    return 0;
}
