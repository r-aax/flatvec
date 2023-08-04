#include "stdafx.h"

#include "array_manager.h"

namespace fv
{
    /// <summary>
    /// Set elements with random numbers.
    /// </summary>
    /// <param name="lo">Lo value.</param>
    /// <param name="hi">Hi value.</param>
#ifdef LINUX_GCC_BUILD
    template<>
#endif
    void ArrayManager<float>::generate_random(float lo,
                                              float hi)
    {
        // See examples in:
        // https://www.techiedelight.com/generate-random-float-value-in-cpp/

        std::random_device d;
        std::default_random_engine gen(d());
        std::uniform_real_distribution<float> dist(lo, hi);

        for (int i = 0; i < size; i++)
        {
            data[i] = dist(gen);
        }
    }
}
